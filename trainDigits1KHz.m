%% Solve for a neural network
clear all; close all; clc;
set(0,'DefaultFigureWindowStyle','docked'); % Figures are docked when generated
SHOW_PROCESS = 1;
SHOW_PROCESS2 = 0;
%% load spike data
tic();

load('digit_spikes1KHz.mat'); % loads variable spike_data
[~,Nutt,Nsubjects] = size(spike_data);

Nin  =length(spike_data(1,1,1).SP(:,1));
Fs =  1E3;
%% Give neural network input
Ndgt =4;
Nout =Ndgt; % output classification neurons %resSize,w,r0,k0,f_inhibit,tau,
[X,Xn,Tau,W,R,E]=createNetwork([20 4 4],[3 6; -2 -2],2,[0.45 0.3;0.6 0.3],0.2,2E-3);

Nres = length(E); % Neurons in the reservoir
[X,Xn,Tau,W,R,E] = full_network(X,Xn,Tau,W,R,E,Nout);
%display_network(X,Xn,Tau,W,R,E,Nout);

N = length(E);    % Total neurons
Nsyn = length(X);
[Yn,inv_map] = sort(X); % just get inv_map
Y = Xn(inv_map);
SynEx = find(E(X)>0); % Synapse Id's coming from excitatory neurons
SynIn = find(E(X)<0); % Synapse Id's coming from inhibitory neurons
%must have connectivity information
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
%ws = 3000;
I0 = 100E-7; % Synaptic current unit
VT = 20E-3; %  Threshold voltage(V)
tau_m = 32E-3; % LIF time constant
tauL = 20E-3; % STDP time constant
tau = 4E-3; % decay time constant of membrane integration
tauc = 64E-3; % decay time constant of calcium concentration
alpha = 2; % ratio between tau & tauS
tauS = tau/alpha;% decay time constant of synaptic current
K0 = 1/(exp(-log(alpha)/(alpha-1))-exp(-log(alpha)/(alpha-1))^alpha); % PSP norm constant
K = @(t) K0*(exp(-t/tau)-exp(-t/tauS)).*heaviside(t); % Normalised PSP

%% Initial values, variable initialization
V =EL*ones(1,N); % These variables are function of time for every neuron % defined for all neurons (Nres+Nout)
Iapp = zeros(1,N); % External applied current for each neuron % defined for all neurons (Nres+Nout)
Isyn = zeros(1,N); % Input synaptic current for each neuron % defined  for all synapses except input synapses from BSA
tSpike = -10*ones(1,N); % recent spike timings for each neuron (-Inf) %XXX : initial value of tSpikes (must be ideally -Inf) % all neurons (Nres+Nout)
PSP = zeros(1,Nsyn); % post synaptic potential due to Weights
c = zeros(1,Nout); % calcium concentration of output neurons;
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
A = (synaptic_connection_matrix(Xn,N));
C = sparse(applied_current_matrix(Nin,Nres,4));
Iteach = zeros(Nout,1);
%% Print information about network
if(SHOW_PROCESS);
    sim_fig = figure('name','Simulation');
    fprintf('Number of Synapses/Neurons = %i/%i(%i)\r\n',Nsyn,N,floor(Nsyn/N));
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
    sim_fig2 = figure('name','Sim2');
end
%figure('name','Input');
% display_input(tt,Iapp_in/max(Iapp_in(:)));
% title('Iapp_in');

%% create teaching current and its parameters
c_theta = 5;
c_delta = 1;
It0 = 1000000; % teaching current unit
%% Weight Parameters
C_th = c_theta;
del_C = 4.5;
Wmin = -8;
Wmax = 8;

%assumed params
del_W = 0.05;
%p_inc = 0.9;
%p_dec = 0.9;
%% Time simulation parameters
dt = 1/Fs;
T = 1; % T seconds
k = 0; % time counter (goes till k_max)
%% Simulate
az=0;

tSpikeLog = logical([]);
cLog =[];
IteachLog = [];
wLog = [];
vLog = [];
%tSpikeLog = logical(zeros(N,1)); % XXX: use only if required
%VLog = zeros(N,length(tt));
%IsynLog = zeros(N,length(tt));
%tSpikeALog = logical(zeros(Nsyn,length(tt)));
t = 0;
%%


training_epochs = 1 ;

training_Nutt = 4;
training_Nsub = 15;
testing_Nsub = 16;
testing_Nutt =  4;
result =    zeros(Ndgt, Nout, testing_Nutt,testing_Nsub);
spike_sum = zeros(Ndgt,Nres,testing_Nutt,testing_Nsub);
toc();
for epoch = 1:training_epochs+1
    if(epoch<=training_epochs)
        nutt_ = training_Nutt;
        nsub_ = training_Nsub;
    else
        nutt_ = testing_Nutt;
        nsub_ = testing_Nsub;
    end
    for digit = 1:Nout
        for sub = 1:nsub_
            for utt = 1:nutt_
                fprintf('%i/%i/%i/%i\r\n',epoch,sub,digit,utt);
                % Create/Get input training data spikes and get applied current
                spike_pattern = spike_data(digit,utt,sub).SP;
                [Nin,k_max] = size(spike_pattern);
                Iin  = double(spike_pattern);
                Iin_res = 100E-4*C*Iin; % reservoir synaptic current unit
                c_d = zeros(1,Nout); % set for desired output neuron for training input
                c_d(digit) = 1;
                %figure(sim_fig2); subplot(211);imagesc(Iin_res); title('Input Reservoir current');
                for k = 1:k_max %time in
                    t = t+dt;
                    
                    Iteach(:) = 0;
                    if (epoch<=training_epochs) % train and test simultaneously
                        %Iteach(:) = -It0;
                        Iteach(find((c_d==1).*(c<c_theta+c_delta))) = It0;
                        Iteach(find((c_d==0).*(c>c_theta-c_delta))) = -It0;
                        
                    end
                    
                    Iapp = [Iin_res(:,k);Iteach]';
                    spikeArrived(:) = 0;
                    spikeArrived(Syn(tSa>t-dt & tSa<=t))=1;%all synapse spikes arriving in time (t-dt,t]
                    tSaR(Syn(tSa>t-dt & tSa<=t)) =  t; % update recent time arrival entries
                    Syn(tSa<=t) = []; % first  remove entries for t<=t from synapse ID
                    tSa(tSa<=t) = []; % remove older processed spiked entries <=t
                    v1 = v1*(1-dt/tau)+spikeArrived';
                    v2 = v2*(1-dt/tauS)+spikeArrived';
                    PSP = K0*(v1-v2);
                    Isyn = I0*(A*(W.*PSP))';               % XXX: input synaptic current (correct) X-(W)->Xn
                    
                    %   V = V *(1-dt*gL/C)+ gL*EL*dt/C+(Iapp+Isyn)*dt/C;%XXX : LIF model (correct)
                    V = V*(1-dt/tau_m)+(Iapp+Isyn);
                    V((t-tSpike<=Rp))=EL; % Refractory period      % XXX : correct (refractory period condition) placement in program?
                    %spikedNeurons = find(V>VT);                 % XXX : get spiked neurons (correct)
                    reservoirSpiked = (V(1:Nres)>VT);
                    outputSpiked = (V(Nres+1:Nres+Nout)>VT); % get output  spiked neurons
                    spikedNeurons = find([reservoirSpiked outputSpiked]);
                    tSpike(spikedNeurons) = t;
                    V(spikedNeurons) = EL;                      % XXX : % correct
                    V(V<EL) = EL;
                    % update calcium concentration due to spiked neurons
                    c = c*(1-dt/tauc)+outputSpiked; % define output spiked
                    
                    if(~isempty(spikedNeurons))    % change weights for spiked neurons
                        tSa = [tSa t+Tau(ismember(X,spikedNeurons))']; %tSa add  time of arrival only of spiked Neurons
                        Syn = [Syn find(ismember(X,spikedNeurons))']; % Syn add information of newly spiked neurons
                    end
                    if(epoch<=training_epochs)
                        spike_sum(digit,:,utt,sub)=spike_sum(digit,:,utt,sub)+reservoirSpiked;
                    else
                        result(digit,:,utt,sub) = result(digit,:,utt,sub)+outputSpiked;
                    end
                    %{
                    if (epoch<=training_epochs) % train
                      
                        pos_upd_id = find(c>C_th & c<C_th+del_C)+Nres;
                        neg_upd_id = find(c<C_th & c>C_th-del_C)+Nres;
                        syn_id_r = ismember(X,spikedNeurons); % replace this function by a faster function
                        syn_id_u = ismember(Xn,pos_upd_id);
                        syn_id_d = ismember(Xn,neg_upd_id);
                        inc_ids = syn_id_r.*syn_id_u;
                        dec_ids = syn_id_r.*syn_id_d;
                        %generating the probability arrays
                        %p_inc_arr = (rand(1,length(X))<p_inc)';
                        %p_dec_arr = (rand(1,length(X))<p_dec)';
                        %W = W + del_W*(inc_ids.*p_inc_arr - dec_ids.*p_dec_arr);
                        W = W + del_W*(inc_ids - dec_ids);
                        W(W>Wmax) = Wmax;
                        W(W<Wmin) = Wmin;
                        
                        end
                    %}
                    % XXX : caused by inhibitory neurons
                    
                    %% Log  data
                    %if(epoch>training_epochs )
                    %   figure(9); plot(V(Nres+1:end));
                    %Wout = W(end-Nout*Nres+1:end);
                    %Wout = reshape(Wout,Nres,Nout);
                    %imagesc(Wout);
                    %  drawnow;
                    %end
                    if(k==1)
                        Wout = W(end-Nout*Nres+1:end);
                        Wout = reshape(Wout,Nres,Nout);
                        wLog(end+1,:) = mean(Wout,1);
                    end
                    if(epoch>=training_epochs-1)
                        tSpikeLog(spikedNeurons,end+1) = 1;
                        cLog(:,end+1)=c;
                        IteachLog(:,end+1) = Iteach;
                        %vLog(end+1,:) = V;
                        
                        
                    end
                    %tSpikeALog(spikeArrived,k) = 1;
                    %IsynLog(:,k) = Isyn;
                    %VLog(:,k) = V;
                    %
                    %% Display information
                    
                    if(k==k_max && SHOW_PROCESS2)
                        figure(sim_fig2);subplot(212); hold on;
                        plot(sort(W(end-Nout*Nres+1:end))); drawnow;
                    end
                    if(SHOW_PROCESS)
                        figure(sim_fig);
                        subplot(visual_subplot);hold off;
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
                        if(k==1)
                            subplot(Vex_hist_subplot);
                            %[N_W,W_edge] = histcounts(V(E>0),100);
                            %bar(W_edge(2:end),log10(N_W)); xlim([EL VT]);
                            plot(Iteach);
                            subplot(Vin_hist_subplot);
                            %[N_W,W_edge] = histcounts(V(E<0),100);
                            %bar(W_edge(2:end),log10(N_W)); xlim([EL VT]);
                            plot(c);
                        end
                        drawnow;
                        az = az+0.2;
                    end
                    
                end
                if(epoch<=training_epochs)
                    % calculate weights  directly
                    W(end-Nout*Nres+1:end) = (6E4*reshape(pinv(spike_sum(:,:,1,1)),Nout*Nres,1)+W(end-Nout*Nres+1:end))/2;
                end
            end
        end
        
    end
    
    
end

toc();
%% Final results (if any)
tt =(0:length(tSpikeLog(1,:))-1)/Fs;
figure('name','Results')
raster_plot(tSpikeLog,tt,E);
title('raster plot');
figure('name','Calcium Concentration')
plot(tt,cLog'); hold on;
%figure('name','Teach Current');
plot(tt,IteachLog'.*max(cLog(:))/It0/5);
line([0 max(tt)],[c_theta c_theta]);
line([0 max(tt)],[c_theta+del_C c_theta+del_C],'LineStyle','--');
line([0 max(tt)],[c_theta-del_C c_theta-del_C],'LineStyle','--');
line([0 max(tt)],[c_theta+c_delta c_theta+c_delta],'LineStyle','--');
line([0 max(tt)],[c_theta-c_delta c_theta-c_delta],'LineStyle','--');
figure('name','Mean Weights');
plot(wLog);
figure('name','Potential');
imagesc(vLog');
%%
% figure('name','result');
% [x,y,z] = ndgrid(1:size(result,1), 1:size(result,2), 1:size(result,3));
% pointsize = result(:,:,:,2);
% scatter3(x(:), y(:), z(:),5*pointsize(:)+1,'o','filled');
[~,id_recognition]=max(result,[],2);
%% figure
toc();
%% Play sound
load train
sound(y,Fs)


