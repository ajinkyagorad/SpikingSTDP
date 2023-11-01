%% Solve for a neural network
clear all; close all; clc;  
set(0,'DefaultFigureWindowStyle','docked'); % Figures are docked when generated
SHOW_PROCESS = 0;
SHOW_PROCESS2 = 1;
SAVE_VIDEO = 0;
tSpikeLogX= [];
%% Give neural network input
N = 300;
R = rand(N,3);
[X,Xn,Tau,W,R,E]=createNetworkD([5 5 5],[1 1 1],[2 1; -2 -1],0.2,[1 1; 1 1],0.2,2E-3);
%[X,Xn,Tau,W,R,E]=createNetworkM([10 10 1],[2 1; -2 -1],2,[0.4 0.3; 0.2 0],0.2,0);%(resSize,w,r0,k0,f_inhibit,tau,show)
%[X,Xn,Tau,W,R,E] = p2Network(25,0.2,3 ,-6,20);%spatialNetwork2(200,0.2,3,1,120,0.2,1)
Nres = length(E); % Neurons in the rese rvoir
%doing full network out neuron pool
% hopefully will find characteristics in the input
Nout = 0;
%[X,Xn,Tau,W,R,E] = full_network(X,Xn,Tau,W,R,E,Nout);
display_network(X,Xn,Tau,W,R,E,Nout);
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
Rp = 10E-3; %Sec
I0 = 1*3500E-7; % Synaptic current unit
Iapp0 = 50E-3; % Input Synapse current unit
VT = 20E-3; %  Threshold voltage(V)
tau_m = 10E-3; % LIF time constant
tauL = 50E-3; % STDP time constant
tau = 32E-3; % decay time constant of membrane integration
tauc = 64E-3; % decay time constant of calcium concentration
alpha = 2; % ratio between tau & tauS
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

fprintf('Number of Synapses/Neurons = %i/%i(%i)\r\n',Nsyn,N,floor(Nsyn/N));
if(SHOW_PROCESS)
    sim_fig = figure('name','Simulation');
    
    visual_subplot = subplot(4,4,[1 2 5 6]);
    h_neurons=scatter3(R(:,1),R(:,2),R(:,3),20*0.001,ColorDef); hold on;
    %h_synapse=scatter3((R(Xn,1)+ R(X,1))/2,(R(Xn,2)+R(X,2))/2,(R(Xn,3)+ R(X,3))/2,'x');
    %h_synapse2=plot3([R(Xn,1)'; R(X,1)'],[R(Xn,2)'; R(X,2)'],[R(Xn,3)'; R(X,3)']);
    %h_synapse.SizeData= 0.1;
    raster_subplot = subplot(4,4,[9 10 13 14]); hold on;
    W_hist_subplot = subplot(4,4,[11 12 15 16]);
    Vex_hist_subplot = subplot(4,4,[3 4]);
    Vin_hist_subplot = subplot(4,4,[7 8]);
end
if(SHOW_PROCESS2)
    sim_fig2 = figure('name','Simulation2');
    sim2.r = 4;
    sim2.c = 6;
    %sim2.s = reshape([1:sim2.r*sim2.c],sim2.c,sim2.r)';
    sim2_raster_plot = subplot(sim2.r,sim2.c,[1 2 3 4 7 8 9 10]);title('rasterPlot');
    sim2_raster_sum1 = subplot(sim2.r,sim2.c,[13 14 15 16 19 20 21 22]);title('Sum1');
    sim2_raster_sum2 = subplot(sim2.r,sim2.c,[5 6 11 12]);title('sum2');
    sim2_param_plot = subplot(sim2.r,sim2.c,[17 18 23 24]); title('Weights');
end
%figure('name','Input');
% display_input(tt,Iapp_in/max(Iapp_in(:)));
% title('Iapp_in');

%% Time simulation parameters
Fs = 1000;
dt = 1/Fs;
Tper = 1; % Time of each period for repeated input
Tsim = 20*Tper; % T seconds

K = Tsim*Fs;
k = 0; % time counter (goes till k_max)
%% Movie Parameters
% ref : http://kawahara.ca/matlab-make-a-movievideo/
if(SAVE_VIDEO)
    writerObj = VideoWriter('media\Sim2.avi'); % Name it.
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

for k = 1:K %time in
    t = t+dt;
    k_ = mod(k,Tper*Fs);
            Iapp = zeros(1,N);
    if(k_<=0.1*Fs)%||k_==0.06*Fs||k_==0.20*Fs) % give impulse to only excitatory neurons
        Iapp(1:10) =Iapp0*rand(10,1);
    end
    Iapp(E<0) = 0;
    spikeArrived(:) = 0;
    spikeArrived(Syn(tSa>t-dt & tSa<=t))=1;%all synapse spikes arriving in time (t-dt,t]
    tSaR(Syn(tSa>t-dt & tSa<=t)) =  t; % update recent time arrival entries
    Syn(tSa<=t) = []; % first  remove entries for t<=t from synapse ID
    tSa(tSa<=t) = []; % remove older processed spiked entries <=t
    v1 = v1*(1-dt/tau)+spikeArrived';
    v2 = v2*(1-dt/tauS)+spikeArrived';
    PSP = K0*(v1-v2);
    Isyn = I0*(A*(W.*PSP))';               % XXX: input synaptic current (correct) X-(W)->Xn
    
    %V = V *(1-dt*gL/C)+gL*EL*dt/C+(Iapp+Isyn)*dt/C;%XXX : LIF model (correct)
    V = V*(1-dt/tau_m)+(Iapp+Isyn);
    V(find(t-tSpike<=Rp))=EL; % Refractory period      % XXX : correct (refractory period condition) placement in program?
    spikedNeurons = find(V>VT);                 % XXX : get spiked neurons (correct)
    not_spikedNeurons = find(V<=VT);
    tSpike(spikedNeurons) = t;
    V(spikedNeurons) = EL;                      % XXX : % correct
    V(V<EL) = EL;    % XXX : caused by inhibitory neurons

    
    
    if(~isempty(spikedNeurons))    % change weights for spiked neurons
        tSa = [tSa t+Tau(ismember(X,spikedNeurons))']; %tSa add  time of arrival only of spiked Neurons
        Syn = [Syn find(ismember(X,spikedNeurons))']; % Syn add information of newly spiked neurons
        
        % following is weight change rules (comment out if not required)
        id_p = find(ismember(Xn,spikedNeurons)); % incoming connection
        id_d = find(ismember(Yn,spikedNeurons)); % outgoing connections
        id_p_ex = intersect(id_p, SynEx);    % only excitatory synapse are potentiated
        id_d_ex = intersect(id_d, SynEx);    % only excitatory synapse are depreciated
        id_p_in = intersect(id_p, SynIn);    % only inhibitory synapse are potentiated
        id_d_in = intersect(id_d, SynIn);    % only inhibitory synapse are depreciated
        % Change explicitly by pre-post spike time interval
        W(id_p_ex) = W(id_p_ex).*(1-0.01*exp(-(t-tSaR(id_p_ex))/tauL));
        W(id_d_ex) = W(id_d_ex).*(1+0.01*exp(-(t+Tau(id_d_ex)-tSpike(X(id_d_ex))')/tauL));
        W(id_p_in) = W(id_p_in).*(1+0.01*exp(-(t-tSaR(id_p_in))/tauL));
        W(id_d_in) = W(id_d_in).*(1-0.01*exp(-(t+Tau(id_d_in)-tSpike(X(id_d_in))')/tauL));
        
        
        % Find deltaW vs delta t
        %limit weigths for both excitatory and inhibitory from both sides
        W(W>Wmax) = Wmax;
        W(W>0 & W< Wmin) = Wmin;
        W(W<0 & W>-Wmin) = -Wmin;
        W(W<-Wmax) = -Wmax;
        % Weights increment
        % increase the input weights of the neurons which did not fired
        
        W(find(ismember(Xn,not_spikedNeurons))) = W(find(ismember(Xn,not_spikedNeurons)))+0.01;
        % update delay also
        %Tau = (1+20-abs(W)*20/Wmax)*1E-3;
    end
    
    
    %% Log  data
    tSpikeLog(spikedNeurons,end+1) = 1;
    %VLog(:,k) = V;
    
    
    %% Display information
    if(mod(k,10*Tper*Fs)==0 && SHOW_PROCESS2 == 1)
%         if(mod(k,10*Tper*Fs)>=Tper*Fs)        
%             tSpikeLog = logical([]);
%             continue; 
%         end
        %if(k==0) continue; end
        figure(sim_fig2);
        subplot(sim2_raster_plot); 
        raster_plot(tSpikeLog,[],E);ylim([1 N]);
        subplot(sim2_param_plot);
        W_sparse=display_network_param(X,Xn,W,N);
        subplot(sim2_raster_sum1);
        plot(1:length(tSpikeLog(1,:)),sum(tSpikeLog,1));
        subplot(sim2_raster_sum2);
        plot(sum(tSpikeLog,2),1:length(tSpikeLog(:,1))); ylim([1 N]);
        drawnow;
        
        tSpikeLog = sparse(logical(zeros(N,0)));
        if(SAVE_VIDEO)
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end
    end
    
    if(SHOW_PROCESS)
        figure(sim_fig);
        subplot(visual_subplot); title(mat2str(t));hold off;
        if(1)
            %Spikes(k,spikedNeurons) = 1;
            Color = ColorDef;
            r = 0.001+abs(V-EL)/(VT-EL);
            
            if(~isempty(spikedNeurons))
                Color(spikedNeurons(E(spikedNeurons)>0),2) = 1; % excitatory neurons
                Color(spikedNeurons(E(spikedNeurons)>0),3) = 0;
                Color(spikedNeurons(E(spikedNeurons)<0),2) = 1; % inhibitory neurons
                r(spikedNeurons) = 1;
            end
            
            Color_W=squeeze(ind2rgb(255*abs(W-Wmin)/(Wmax-Wmin),'colorcube'));
            view(az,30);
            set(h_neurons,'SizeData',20*r,'CData',Color);
            %set(h_synapse,'SizeData',4*abs(W-Wmin)/(Wmax-Wmin)+0.001);%,'CData',Color_W);
            subplot(raster_subplot);
            plot(t*ones(1,length(spikedNeurons(E(spikedNeurons)>0))),spikedNeurons(E(spikedNeurons)>0),'.','MarkerSize',0.01,'Color','b');
            plot(t*ones(1,length(spikedNeurons(E(spikedNeurons)<0))),spikedNeurons(E(spikedNeurons)<0),'.','MarkerSize',0.01,'Color','r');
        end
        if(1)
            subplot(W_hist_subplot);
            [N_W,W_edge] = histcounts(W,500);
            bar(W_edge(2:end),log10(N_W));% xlim([-Wmax Wmax]);
        end
        if(1)
            subplot(Vex_hist_subplot);
            [N_W,W_edge] = histcounts(V(E>0),100);
            bar(W_edge(2:end),log10(N_W)); xlim([EL VT]);
            subplot(Vin_hist_subplot);
            [N_W,W_edge] = histcounts(V(E<0),100);
            bar(W_edge(2:end),log10(N_W)); xlim([EL VT]);
        end
        % To show directional graph
        %         figure(1)
        %         G = digraph(X,Xn,W);
        %         plot(G,'XData',R(:,1),'YData',R(:,2),'ZData',R(:,3),'EdgeLabel',G.Edges.Weight)
        drawnow;
        az = az+0.2;
    end
    
end
%% MOVIE END 
if(SAVE_VIDEO)
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
%% 


