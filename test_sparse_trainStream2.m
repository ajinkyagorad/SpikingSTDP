%% Solve for a neural network
% modify to display activity in graph
clear all; close all; clc;
set(0,'DefaultFigureWindowStyle','docked'); % Figures are docked when generated
SHOW_PROCESS = 0;
USE_GPU = 0; % 1/0 use GPU for computation, gpuArrays will be used
USE_FULL_ARRAY_GPU = 1; %1/0 use full matrix instead of sparse for GPU
RECORD_VIDEO = 0;
LOG_DATA.DATA = 1;
LOG_DATA.SINGLE_NEURON = 1;
STDP_RULE = 1;

%% load spike data
tic();
%load('rec2_spike.mat'); % loads variable spike_data
%spike_data = audiodata;
%fprintf('playing audio to be fed\r\n');
%soundsc(audiodata.sig);
%Nin  =length(spike_data.SP(:,1));
Fs =  1E3; % must be
%% Give neural network input

Nout =0; % output clxassification neurons %resSize,w,r0,k0,f_inhibit,tau,
[X,Xn,Tau,W,R,E]=createNetwork([7 7 7],[0.5 0.5;-0.2 -0.2],2,[0.5,0.5; 0.5 0.5],0.2,0);
Nres = length(E); % Neurons in the reservoir
%[X,Xn,Tau,W,R,E] = full_network(X,Xn,Tau,W,R,E,Nout);
%display_network(X,Xn,Tau,W,R,E,Nout);

%% Make all the network information in sparse matrices
N = Nres+Nout;
W = sparse(X,Xn,W,N,N);
A = sparse(X,Xn,1,N,N);
Tau = sparse(X,Xn,Tau,N,N);
TauK = ceil(Tau*1E3);
Nsyn = length(X);
%% Find strongly connected components of the graph
SCC = tarjan(A')
% singly connected components could be taken as output
% input nodes for one of the node in multiple connected nodes
input_nodes = [];
output_nodes = [];
for i = 1:length(SCC)
    if(length(SCC{i})>1)
        input_nodes(end+1:end+length(SCC{i})) = SCC{i};
    else
        output_nodes(end+1) = SCC{i};
    end
end


%% Parameters of the LIF neuron
% parameters ending with K (format *K) have their  values for milliseconds
% timestamp simulation
% EL = 0 ; always
VT = 1;  % parameter to be fixed to 1, so other variables could be changed
RpK = 5;
tau_mK = 50; % LIF time decay time constant of membrane integration
tau1K = 50;   % second order PSP time constant 1 (PSP : post synaptic potential)
tau2K = 12;   % second order PSP time constant 2
tauLK = 30;
alpha = tau1K/tau2K;
K0 = 1/(exp(-log(alpha)/(alpha-1))-exp(-log(alpha)/(alpha-1))^alpha); % PSP norm constant
%% Weight update parameters
Wmin = 0;
Wmax = 2;
%% Variables of simulation (not logged)
V = zeros(N,1);
Iin = zeros(N,1);
Iapp = zeros(N,1);
v1 = sparse(X,Xn,0,N,N);
v2 = sparse(X,Xn,0,N,N);
v3 = sparse(X,Xn,0,N,N);
v4 = sparse(N,1);
if(USE_FULL_ARRAY_GPU && USE_GPU)
    v1 = full(v1);
    v2 = full(v2);
    v3 = full(v3);
    v4 = full(v4);
    W = full(W);
    Wmax = full(Wmax);
    A = full(A);
end
%spikeArrival= sparse(N,N); % history matrix
lastSpike = -1000*ones(N,1);
spikes = zeros(N,1);
x = zeros(0,1);
xn = zeros(0,1);
ta = zeros(0,1); % because matrix concat error for two empty matrices during execution
%C = sparse(applied_current_matrix(Nin,Nres,1));
% make all gpuArray
if(USE_GPU)
    V = gpuArray(V);
    Iin = gpuArray(Iin);
    Iapp = gpuArray(Iapp);
    v1 = gpuArray(v1);
    v2 = gpuArray(v2);
    lastSpike=gpuArray(lastSpike);
    spikes = gpuArray(spikes);
    x = gpuArray(x);
    xn = gpuArray(xn);
    ta = gpuArray(ta);
    W = gpuArray(W);
    Wmax = gpuArray(Wmax);
    TauK = gpuArray(TauK);
    A = gpuArray(A);
end
%% Print information about network
if(SHOW_PROCESS)
     % variables for display
    ColorDef = zeros(N,3); % color coding of neurons
    ColorDef(find(E>0),:) = ones(length(find(E>0)),3).*repmat([0 0 1],length(find(E>0)),1); % excitatory neuronsosn coded as blue
    ColorDef(find(E<0),:) = ones(length(find(E<0)),3).*repmat([1 0 0],length(find(E<0)),1); % inhibitory neurons coded red
    
    G = digraph(sparse(X,Xn,1,N,N));
    %h =plot(G,'XData',R(:,1),'YData',R(:,2),'ZData',R(:,3),'EdgeAlpha',0.2); hold on;
    h =plot(G,'Layout','layered','Sinks',20,'EdgeAlpha',0.2); hold on;
    %https://in.mathworks.com/help/matlab/ref/graph.plot.html
    %h =plot(G,'Layout','subspace3','Dimension',length(SCC),'EdgeAlpha',0.2,'ButtonDownFcn',@cbHighlightNode); % autoposition plot
    bgColor = [0.01 0.01 0.01];
    set(gca,'Color',bgColor);set(gcf,'color',bgColor);
    set(0,'DefaultFigureWindowStyle','normal');
    %set(gcf, 'Position', get(0, 'Screensize'));
    %set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    
    hold on;
    %light('Position',[-1 0 0],'Style','local')
    %axis equal
    camproj('perspective')
    
end
%% Simulation time parameters
% Fs must be defined above
dt = 1/Fs;
Tsim = 2; % simulation time in seconds
K = Tsim*Fs;
%% Movie Parameters
% ref : http://kawahara.ca/matlab-make-a-movievideo/
if(RECORD_VIDEO)
    
    writerObj = VideoWriter([int2str(now*1E6) '.avi'],'MPEG-4'); % Name it.
    fprintf('Video writing to file %s\r\n',[int2str(now) '.avi']);
    open(writerObj);
end
%% Log parameters
tSpikeLog = sparse(zeros(N,0));
LOG.meanW = [];
if(LOG_DATA.SINGLE_NEURON)
    LOG.preneuron = find(E>0); % XXX : could put preneurons for this neurons also ( this be 'neuron', pre->neuron->post)
    LOG.preneuron = LOG.preneuron(1);
    LOG.postneuron = Xn(X==LOG.preneuron);
    LOG.postneuron = LOG.postneuron(1);
    LOG.v3 = [];
    LOG.v4 = [];
    LOG.Vpre = [];
    LOG.Vpost = [];
    LOG.W = [];
end
%% Simloop
display.az = 0;
t = 0;
tic();
for k = 1:K
    t = t+dt;
    if(mod(k,200)==0); fprintf('real2sim time ratio %2.2f\r\n',toc/t); end
    if(mod(k,200)==0)
        Iin(:)=0;
        V(:) = 0;
        v1(:) = 0;
        v2(:) = 0;
        v3(:) = 0;
        v4(:) = 0;
        x = []; xn=[]; ta=[];
    end
    if(1)
        Iapp(input_nodes) =10*rand(length(input_nodes),1)<0.5;%*rand(16,1);% give impulse current if there are no spikes arriving
    else
        Iapp(:) = 0;
    end
    %if(1);Iapp(1:5) =  rand(5,1);end
    V = V*(1-1/tau_mK)+Iapp+Iin;
    V(k-lastSpike<RpK) = 0;
    spikedNeurons = find(V>VT);
    spikes(:) = 0;spikes(spikedNeurons) = 1;
    lastSpike(spikedNeurons) = k; % for refractory period
    V(V>VT) = 0;
    V(V<0) = 0;
    [x_,xn_,ta_] = find(diag(spikes)*TauK);
    if(isempty(ta)) % Could it be done in loop faster(few spikes)?, can input spikes be given here!?
        x = x_; xn = xn_; ta = k+ta_;
    else
        x = [x;x_];  xn = [xn;xn_]; ta = [ta;k+ta_];
    end
    spikeArrived = full(sparse(x(ta==k),xn(ta==k),1,N,N));
    x(ta<=k) = []; xn(ta<=k) = []; ta(ta<=k) = [];
    
    % v1,v2 determine the PSP shape
    % v3,v4 determine the shape of the STDP
    
    v1 = v1*(1-1/tau1K)+spikeArrived;
    v2 = v2*(1-1/tau2K)+spikeArrived;
    v3(find(spikeArrived)) = 0;
    v3 = v3*(1-2/tauLK)+spikeArrived;
    v4(find(spikes)) = 0; % reset values due to earlier spikes
    v4 = v4*(1-1/tauLK)+spikes;
    %     [I, J, v12] = find(v1-v2);
    %     [~, ~, w] = find(W);
    %     C = sparse(I, J, v12.*w);
    Isyn = K0*(v1-v2); % these must be dense(full) matrices for gpu
    Iin = 1E-2*sum(Isyn.*W,1)';
    %Iin = Iin./(max(Iin(:))+10);
    %STDP rule
    if(STDP_RULE)
        % d = (-diag(v4)*spikeArrived+v3*diag(spikes));
        % d(d==0) = [];
        W = W.*(1+0.1.*(Wmax-W).*(-spikeArrived*diag(v4)+v3*diag(spikes)));
        %W = W.*(1+0.001.*(Wmax.^2-W.^2).*(-spikeArrived*diag(v4)+0.1));
        % diag(spikes)*v3 gets for pre spike arrival after post spike &  v3*diag(spikes) gets for pre spike arrival before post spike
        %figure(9); hist(d(:),100); drawnow;
    end
    %% Log Data
    if(LOG_DATA.DATA)
        tSpikeLog(spikedNeurons,end+1) = 1;
        LOG.meanW(end+1) = sum(sum(gather(W),2),1)/N^2;
        if(LOG_DATA.SINGLE_NEURON)
            LOG.v3(:,end+1) = gather(v3(LOG.preneuron,LOG.postneuron));
            LOG.v4(:,end+1) = gather(v4(LOG.postneuron));
            LOG.Vpre(:,end+1) = gather(V(LOG.preneuron));
            LOG.Vpost(:,end+1) = gather(V(LOG.postneuron));
            LOG.W(:,end+1) = gather(W(LOG.preneuron,LOG.postneuron));
        end
    end
    %% Display INFO
    if(SHOW_PROCESS )
        %Display it on graph
        %h.LineWidth = 10*abs(W)/Wmax;
         %Spikes(k,spikedNeurons) = 1;
            Color = ColorDef;
            r = gather(0.001+abs(V/VT));          
            if(~isempty(find(spikedNeurons)))
                Color(spikedNeurons(E(spikedNeurons)>0),2) = 1; % excitatory neurons
                Color(spikedNeurons(E(spikedNeurons)>0),3) = 0;
                Color(spikedNeurons(E(spikedNeurons)<0),2) = 1; % inhibitory neurons
                r(spikedNeurons) = 1;
            end
            h.NodeCData = v4;
            h.EdgeCData = (v1(A~=0));
            h.LineWidth = (5*abs(W(A~=0)));
            h.MarkerSize = (5*r);
            title(['t = ' num2str(t)],'Color','w');
        if(RECORD_VIDEO)
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
        end
        display.az = display.az+0.2;
        drawnow;
    end
end
toc();
%% MOVIE END
if(RECORD_VIDEO)
    close(writerObj);
end
%% Display logged data
if(LOG_DATA.SINGLE_NEURON && LOG_DATA.DATA)
    figure('name','loggedData');
    hold on;
    %plot(-[diff(LOG.Vpre)'+0.5 diff(LOG.Vpost)']);    plot(diff(LOG.v3)');
    plot(1+4*(LOG.W'-min(LOG.W))/(max(LOG.W)-min(LOG.W)));
    plot(LOG.v3');
    plot(LOG.v4');
    TauK_temp = full(TauK);
    tspikepre = find(tSpikeLog(LOG.preneuron,:));
    tspikepre_arrive = tspikepre+TauK_temp(LOG.preneuron,LOG.postneuron);
    tspikepost = find(tSpikeLog(LOG.postneuron,:));
    plot(tspikepre,0.99*ones(1,length(tspikepre)),'x','Color','b');
    plot(tspikepre_arrive,ones(1,length(tspikepre_arrive)),'^','Color','r');
    plot(tspikepost,ones(1,length(tspikepost)),'*','Color','g');
    legend('Wij','prespike','preNeuronPercievedactivity','postNeuronActivity','spikearrive','postspike'); ylim([0 5]);
    plot_STDP; % plot the STDP of the logged neuron
end
%%
if(LOG_DATA.DATA);
    figure('name','raster_plot');
    raster_plot(tSpikeLog([1:N] ,:),(1:length(tSpikeLog(1,:)))/Fs,E);
end




