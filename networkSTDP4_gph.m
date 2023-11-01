%% Solve for a neural network
clear all; close all; clc;
set(0,'DefaultFigureWindowStyle','docked'); % Figures are docked when generated
SHOW_PROCESS = 1;
RECORD_VIDEO = 0;
tSpikeLogX= [];
%% Give neural network input
%  K = 50; L = 50;
%  X = [];Xn =[]; W=[]; Tau = [];
%  for i = 1:K*L
%      if(mod(i,K)~=0 && i<(L-1)*K)
%          X(end+1) = i; Xn(end+1) = i+1; W(end+1) =2; Tau(end+1) = 5E-3;
%          X(end+1) = i; Xn(end+1) = i+K; W(end+1) =2; Tau(end+1) = 5E-3;
%      end
%      if(mod(i,K)==0 && i<(L-1)*K)
%          X(end+1) = i; Xn(end+1) = i+K; W(end+1) =2; Tau(end+1) = 5E-3;
%      end
%      if(i>(L-1)*K && i<L*K)
%          X(end+1) = i; Xn(end+1) = i+1; W(end+1) =2; Tau(end+1) = 5E-3;
%      end
%  end
% N = K*L;
% X = reshape(X,length(X),1);
% Xn = reshape(Xn,length(Xn),1);
% W = reshape(W,length(W),1);
% Tau = reshape(Tau,length(Tau),1);
% [Xn,sort_id] = sort(Xn);
% X = X(sort_id);
% W = W(sort_id);
% Tau = Tau(sort_id);
% E = rand(N,1); E(E>0.2) = 1;  E(E<=0.2) = -1;
% R = [ mod(0:N-1,K); ((0:N-1)-mod(0:N-1,K))/K;zeros(1,N)]';
N = 100;
R(:,1) = (randn(N,1)+100); % radius
R(:,2) = rand(N,1)*2*pi;% theta
R(:,3) = randn(N,1)/1000; % phi
[R(:,1) R(:,2) R(:,3)]  = sph2cart(R(:,2),R(:,3),R(:,1));
%addpath specialNetworks
%[X,Xn,Tau,W,R,E] = networkDelayDetector();
[X,Xn,Tau,W,R,E]=createNetworkD([4,4,4],[1 1 1],[2 1; -2 -1],2,[1 1; 1 1],0.2,2E-3);
%[X,Xn,Tau,W,R,E]=createNetworkM([3 3 3 3],[2 1; -2 -1],3,[0.4 0.3; 0.2 0],0.2,2E-3);%(resSize,w,r0,k0,f_inhibit,tau,show)
%[X,Xn,Tau,W,R,E] = p2Network(25,0.2,3,-6,20);%spatialNetwork2(200,0.2,3,1,120,0.2,1)
Nres = length(E); % Neurons in the reservoir
%doing full network out neuron pool
% hopefully will find characteristics in the input
Nout = 0;
%[X,Xn,Tau,W,R,E] = full_network(X,Xn,Tau,W,R,E,Nout);
%display_network(X,Xn,Tau,W,R,E,Nout);
N = length(E);    % Total neurons
Nsyn = length(X);
[Yn,inv_map] = sort(X); % just get inv_map
Y = Xn(inv_map);
SynEx = find(E(X)>0); % Synapse Id's coming from excitatory neurons
SynIn = find(E(X)<0); % Synapse Id's coming from inhibitory neurons
%% INFO : till now must have connectivity information
% X,Xn,Tau,W,R,E ,N,Nsyn , Y,Yn
% [X;Xn] pair describes the unidirectional input synaptic connections(fanin)
% [Y;Yn] pair describes outgoing synaptic connections
% Tau,W are delays and initial weights for the given synaptic connections
% R x,y,z coordinates of N neurons
% E is an binary (+1/-1) array indicating exitatory(+1)/inhibitatory(-1) neurons
% N : number of neurons , Nsyn = number of synapses
% X,Xn,Tau,W:(1xNsyn arrays),R = Nx3 matrix, E= 1xN;N = integer, Nsyn = integer
%% Parameters for LIF model of Neural network
C = 300E-12; % Farad
gL = 30E-9; %Siemens
EL = 0; % V
Rp = 2E-3; %Sec
I0 = 1*3500E-7; % Synaptic current unit
Iapp0 = 50E-3; % Input Synapse current unit
VT = 20E-3; %  Threshold voltage(V)
tau_m = 5E-3; % LIF time constant
tauL = 50E-3; % STDP time constant
tau = 4E-3; % decay time constant of membrane integration

alpha = 1.2; % ratio between tau & tauS
tauS = tau/alpha;% decay time constant of synaptic current
K0 = 1/(exp(-log(alpha)/(alpha-1))-exp(-log(alpha)/(alpha-1))^alpha); % PSP norm constant
K = @(t) K0*(exp(-t/tau)-exp(-t/tauS)).*heaviside(t); % Normalised PSP
% dummy weight  information
Wmin = 0;
Wmax = 20;
%% Initial values, variable initialization
V =EL*ones(1,N); % These variables are function of time for every neuron
Iapp = zeros(1,N); % External applied current for each neuron % defined for all neurons (Nres+Nout)
Isyn = zeros(1,N); % Input synaptic current for each neuron % defined  for all synapses except input synapses from BSA
tSpike = -10*ones(1,N); % recent spike timings for each neuron (-Inf) %XXX : initial value of tSpikes (must be ideally -Inf) % all neurons (Nres+Nout)
PSP = zeros(1,Nsyn); % post synaptic potential due to Weights
clear v1 v2
v1(:,1) = zeros(Nsyn,1); % finding PSP,decaying term at rate tau in exp(-t/tau)
v2(:,1) = zeros(Nsyn,1); % finding PSP,decaying term at rate tauS in exp(-t/tauS)
v3(:,1) = zeros(Nsyn,1);
v4(:,1) = zeros(N,1);
spikeArrived = logical(zeros(1,Nsyn));
tSa = [];   % stores future spike time arrival due to activity till time t
Syn = [];   % synapse id corresponding to tSa
tSaR = -10*ones(size(X)); % stores the recent spike time arrival  entries
Spikes = zeros(1,N); % storing latest time of spike for each neuron for all time <t
ColorDef = zeros(N,3); % color coding of neurons
ColorDef(find(E>0),:) = ones(length(find(E>0)),3).*repmat([0 0 1],length(find(E>0)),1); % excitatory neuronsosn coded as blue
ColorDef(find(E<0),:) = ones(length(find(E<0)),3).*repmat([1 0 0],length(find(E<0)),1); % inhibitory neurons coded red
A = synaptic_connection_matrix(Xn,N);


%% Print information about network
if(SHOW_PROCESS)
     % variables for display
    ColorDef = zeros(N,3); % color coding of neurons
    ColorDef(find(E>0),:) = ones(length(find(E>0)),3).*repmat([0 0 1],length(find(E>0)),1); % excitatory neuronsosn coded as blue
    ColorDef(find(E<0),:) = ones(length(find(E<0)),3).*repmat([1 0 0],length(find(E<0)),1); % inhibitory neurons coded red
    
    G = digraph(sparse(X,Xn,1,N,N));
    %h =plot(G,'XData',R(:,1),'YData',R(:,2),'ZData',R(:,3),'EdgeAlpha',0.2); hold on;
    h =plot(G,'Layout','layered','Source',10,'EdgeAlpha',0.3); hold on;
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
    %h2 =plot3(0,0,0);
end
%% Time simulation parameters
Fs = 1000;
dt = 1/Fs;
Tper = 0.5; % Time of each period for repeated input
Tsim = 100*Tper; % T seconds

K = Tsim*Fs;
k = 0; % time counter (goes till k_max)
%% Movie Parameters
% ref : http://kawahara.ca/matlab-make-a-movievideo/
if(RECORD_VIDEO)
    Vid_filename = [int2str(now*1E6)];
     writerObj = VideoWriter(Vid_filename,'MPEG-4'); % Name it.
    fprintf('Video writing to file %s\r\n',[Vid_filename '.mp4']);
    open(writerObj);
end
%% Simulate
az=0;
tSpikeLog = sparse(logical(zeros(N,0)));
cLog =[];
IteachLog = [];
wLog = [];
%tSpikeLog = logical(zeros(N,1)); % XXX: use only if required
%VLog = zeros(N,length(tt));
%IsynLog = zeros(N,length(tt));
%tSpikeALog = logical(zeros(Nsyn,length(tt)));
tic();
t=0;
Delay = -Tper/2*Fs;
for k = 1:K %time in
    t = t+dt;
    k_ = mod(k,Tper*Fs);
    Iapp(1:2) = 0;
    if(k_-Tper/2*Fs==1)
        Iapp(1:5) = 1;
    end
    
    if(k_-Tper/2*Fs==Delay+1)
        Iapp(6:10) = 1;
    end
    
    if(k_==0) % increase delay after a period
        Delay=Delay+5;
    end
    spikeArrived(:) = 0;
    spikeArrived(Syn(tSa>t-dt & tSa<=t))=logical(1);%all synapse spikes arriving in time (t-dt,t]
    tSaR(Syn(tSa>t-dt & tSa<=t)) =  t; % update recent time arrival entries
    Syn(tSa<=t) = []; % first  remove entries for t<=t from synapse ID
    tSa(tSa<=t) = []; % remove older processed spiked entries <=t
    v1(spikeArrived) = 0;
    v2(spikeArrived) = 0;
    v1 = v1*(1-dt/tau)+spikeArrived';
    v2 = v2*(1-dt/tauS)+spikeArrived';
    v3(spikeArrived) = 1;
    v3 = v3*(1-dt/tau);
    v4 = v4*(1-dt/tau);
    PSP = K0*(v1-v2);
    Isyn = I0*(A*(W.*PSP))';               % XXX: input synaptic current (correct) X-(W)->Xn
    
    %V = V *(1-dt*gL/C)+gL*EL*dt/C+(Iapp+Isyn)*dt/C;%XXX : LIF model (correct)
    V = V*(1-dt/tau_m)+(Iapp+Isyn);
    V(find(t-tSpike<=Rp))=EL; % Refractory period      % XXX : correct (refractory period condition) placement in program?
    spikedNeurons = find(V>VT);                 % XXX : get spiked neurons (correct)
    
    tSpike(spikedNeurons) = t;
    V(spikedNeurons) = EL;                      % XXX : % correct
    V(V<EL) = EL;    % XXX : caused by inhibitory neurons
        
    if(~isempty(spikedNeurons))    % change weights for spiked neurons
        v4(spikedNeurons) = 1;
        tSa = [tSa t+Tau(ismember(X,spikedNeurons))']; %tSa add  time of arrival only of spiked Neurons
        Syn = [Syn find(ismember(X,spikedNeurons))']; % Syn add information of newly spiked neurons
    end
        
    
        % following is weight change rules (comment out if not required)
        id_p = find(ismember(Xn,spikedNeurons)); % incoming connection
        id_d = find(ismember(Yn,spikedNeurons)); % outgoing connections
        id_p_ex = intersect(id_p, SynEx);    % only excitatory synapse are potentiated
        id_d_ex = intersect(id_d, SynEx);    % only excitatory synapse are depreciated
        id_p_in = intersect(id_p, SynIn);    % only inhibitory synapse are potentiated
        id_d_in = intersect(id_d, SynIn);    % only inhibitory synapse are depreciated
        % Change explicitly by pre-post spike time interval
        W(id_p_ex) = W(id_p_ex).*(1-0.1.*exp(-(t-tSaR(id_p_ex))/tauL));
        W(id_d_ex) = W(id_d_ex).*(1+0.1.*exp(-(t+Tau(id_d_ex)-tSpike(X(id_d_ex))')/tauL));
        W(id_p_in) = W(id_p_in).*(1+0.1.*exp(-(t-tSaR(id_p_in))/tauL));
        W(id_d_in) = W(id_d_in).*(1-0.1.*exp(-(t+Tau(id_d_in)-tSpike(X(id_d_in))')/tauL));
        
        
        % Find deltaW vs delta t
        %limit weigths for both excitatory and inhibitory from both sides
        W(W>Wmax) = Wmax;
        W(W>0 & W< Wmin) = Wmin;
        W(W<0 & W>-Wmin) = -Wmin;
        W(W<-Wmax) = -Wmax;
    %end
    %}
    %% Log  data
    tSpikeLog(spikedNeurons,end+1) = 1;
    %VLog(:,k) = V;
    
    
    %% Display information
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
            %h.NodeCData = v4;
            h.EdgeCData = v1;
            h.LineWidth = 4*abs(W/Wmax);
            h.MarkerSize = 5*r;
            title(['t = ' num2str(t)],'Color','w');
%             h2.XData = 1:N;
%             h2.YData = 10*V;
%             h2.ZData = -1*ones(1,N);
        if(RECORD_VIDEO)
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
        end
        
        %   az =az+0.2;
        drawnow;
    end
    
end
if(RECORD_VIDEO)
    close(writerObj);
end
toc();
%% Final results (if any)
tt =(0:length(tSpikeLog(1,:))-1)/Fs;
figure('name','Raster')
raster_plot(tSpikeLog,tt,E);
title('raster plot');

%% Post process
% process on raster plot to order neuron postspikes according to second
% spike time

%% Plot neuron sensitivity
% color coded with number of spikes as measure
figure('name','NeuronSensitivity');
sensitivity = sum(tSpikeLog,2);
sensitivity = (sensitivity-min(sensitivity))/(max(sensitivity)-min(sensitivity));
h_neurons=scatter3(R(:,1),R(:,2),R(:,3),sensitivity*40+1,sensitivity,'filled'); colormap copper;



