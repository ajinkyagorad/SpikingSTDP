%% Solve for a neural network
clear all; close all; clc;
set(0,'DefaultFigureWindowStyle','docked'); % Figures are docked when generated
SHOW_PROCESS = 1;
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
[X,Xn,Tau,W,R,E]=createNetworkD([5 5 5],[0 1 1],[0.5 0.5;-0.2 -0.2],2,[0.5,0.5; 0.5 0.5],0.2,0);
Nres = length(E); % Neurons in the reservoir
%[X,Xn,Tau,W,R,E] = full_network(X,Xn,Tau,W,R,E,Nout);
display_network(X,Xn,Tau,W,R,E,Nout);

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
Wmax = 1;
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
    
    % Figures
    sim_fig = figure('name','Simulation');
    fprintf('Number of Synapses/Neurons = %i/%i(%i)\r\n',Nsyn,N,floor(Nsyn/N));
    visual_subplot = subplot(4,4,[1 2 5 6]);
    h_neurons=scatter3(R(:,1),R(:,2),R(:,3),20*0.001,ColorDef); hold on;
    % uncomment following to display synapses
    %h_synapse=scatter3((R(Xn,1)+ R(X,1))/2,(R(Xn,2)+R(X,2))/2,(R(Xn,3)+ R(X,3))/2,'x');
    %h_synapse2=plot3([R(Xn,1)'; R(X,1)'],[R(Xn,2)'; R(X,2)'],[R(Xn,3)'; R(X,3)']);
    %h_synapse.SizeData= 0.1;
    raster_subplot = subplot(4,4,[9 10 13 14]); hold on;
    sim_vis = subplot(4,4,[11 12 15 16]);
    sim_vis2 = subplot(4,4,[3 4 7 8]);
    
end
%% Simulation time parameters
% Fs must be defined above
dt = 1/Fs;
Tsim = 2.5; % simulation time in seconds
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
%% Generate  patterns for input
pattern1=rand(100,length(input_nodes),1)<0.1;
pattern2=rand(100,length(input_nodes),1)<0.1;
figure('name','InputPatterns');
imagesc(cat(3,pattern1,pattern2,zeros(size(pattern1))));
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
    
    
if(mod(k,200)<=100)
        Iapp(input_nodes) =pattern1(mod(k,100)+1,:);%*rand(16,1);% give impulse current if there are no spikes arriving
        if(k<1500) % learn for 1 second (2.5 epochs)
            Iapp(output_nodes(end-3:end-2)) = 2;
            Iapp(output_nodes(end-1:end)) = -2;
        end
    elseif(mod(k,200)>100 && mod(k,200)<=200)
        Iapp(input_nodes) =pattern2(mod(k,100)+1,:);%higher spikes
        if(k<1500)
            Iapp(output_nodes(end-3:end-2)) = -2;
            Iapp(output_nodes(end-1:end)) = +2;
        end
    end
    if(k>1500)
        Iapp(output_nodes(end-3:end-2)) = 0;
            Iapp(output_nodes(end-1:end)) = 0;
    end
    %}
    %if(0);Iapp(1:N) =  1*(rand(N,1)<0.1);end
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
    v1(find(spikeArrived)) = 0;
    v2(find(spikeArrived)) = 0;
    v1 = v1*(1-1/tau1K)+spikeArrived;
    v2 = v2*(1-1/tau2K)+spikeArrived;
    v3(find(spikeArrived)) = 0;
    v3 = v3*(1-1/tauLK)+spikeArrived;
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
        %Wold = W;
        W = W.*(1+0.05.*(Wmax-W).*(-spikeArrived*diag(v4)+v3*diag(spikes)));
        %max(abs(Wold(:)-W(:)))
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
        figure(sim_fig);
        if(1)
            subplot(visual_subplot);hold off;
            title(['t =' num2str(t)]);
            
            %Spikes(k,spikedNeurons) = 1;
            Color = ColorDef;
            r = gather(0.001+abs(V/VT));
            
            if(~isempty(find(spikedNeurons)))
                Color(spikedNeurons(E(spikedNeurons)>0),2) = 1; % excitatory neurons
                Color(spikedNeurons(E(spikedNeurons)>0),3) = 0;
                Color(spikedNeurons(E(spikedNeurons)<0),2) = 1; % inhibitory neurons
                r(spikedNeurons) = 1;
            end
            
            %Color_W=squeeze(ind2rgb(255*abs(W-Wmin)/(Wmax-Wmin),'colorcube'));
            view(display.az,30);
            set(h_neurons,'SizeData',20*r,'CData',Color);
            %set(h_synapse,'SizeData',4*abs(W-Wmin)/(Wmax-Wmin)+0.001);%,'CData',Color_W);
        end
        if(1)
            subplot(raster_subplot); hold on;
            plot(t*ones(1,length(spikedNeurons(E(spikedNeurons)>0))),spikedNeurons(E(spikedNeurons)>0),'.','MarkerSize',0.01,'Color','b');
            plot(t*ones(1,length(spikedNeurons(E(spikedNeurons)<0))),spikedNeurons(E(spikedNeurons)<0),'.','MarkerSize',0.01,'Color','r');
        end
        if(1)
            subplot(sim_vis);
            imagesc(Isyn); hold on;set(gca,'YDir','normal'); plot(N*V,'w'); hold off;
            subplot(sim_vis2);
            imagesc(W);set(gca,'YDir','normal');
            hold on;
            histogram(10*W(W~=0)./Wmax+N/2,200,'FaceAlpha',0.4,'EdgeAlpha',0.1,'EdgeColor','r');
            plot(min(N,length(LOG.meanW))*(1:length(LOG.meanW))/min(N,length(LOG.meanW)),5000*LOG.meanW,'r');hold off;
        end
        drawnow;
        if(RECORD_VIDEO)
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
        end
        display.az = display.az+0.2;
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
    plot(2+4*(LOG.W'-min(LOG.W))/(max(LOG.W)-min(LOG.W)));
    plot(LOG.v3');
    plot(LOG.v4');
    TauK_temp = full(TauK);
    tspikepre = find(tSpikeLog(LOG.preneuron,:));
    tspikepre_arrive = tspikepre+TauK_temp(LOG.preneuron,LOG.postneuron);
    tspikepost = find(tSpikeLog(LOG.postneuron,:));
    plot(tspikepre,0.99*ones(1,length(tspikepre)),'x','Color','b');
    plot(tspikepre_arrive,ones(1,length(tspikepre_arrive)),'^','Color','r');
    plot(tspikepost,ones(1,length(tspikepost)),'*','Color','g');
    legend('Wij','preNeuronPercievedactivity','postNeuronActivity','prespike','spikearrive','postspike'); ylim([0 5]);
    plot_STDP; % plot the STDP of the logged neuron
end
%%
if(LOG_DATA.DATA);
    figure('name','raster_plot');
    raster_plot(tSpikeLog([[1:N] output_nodes],:),(1:length(tSpikeLog(1,:)))/Fs,E([[1:N] output_nodes]));
end




