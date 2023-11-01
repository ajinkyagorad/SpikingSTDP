% Locate neurons in spatial positions with their connections
% No self connection
%
% N: Number of neurons
% D: Neuron effective distance (fraction)
% w0: default weight
% L: dimensions of neuron network cube in m
% V: conduction velocity between neurons (m/s)
%
% displays :
%    -3D scatter plto of network with connections as lines
%    -Histogram of delay for Interneuron synapses
%    -Synapse distribution for neurons
%    -Mean input weight distribution based on excitatory & inhibitory
%       connections
% returns :
%   X :  array of connections for each neuron, can be interpreted as
%      either fan-in or fan-out
%   T : contains delays for each of the connections for each neuron in 
%       array 
%   W : contains weights for each of the connections for each neuron in
%       array
%   R : gives the 3D x,y,z coordinates of the neurons
% 
%   E : array which lists excitatory neurons as +1 & inhibitory as -1


function [X,Xn,T,W,R,E]= spatialNetwork2(N,Deff,w0,l,v,inhibit_f,fig)
if(nargin<1)
    N = 1000;
end
if(nargin<2)
    Deff = 0.05;
end
if(nargin<3)
    w0 = 3000;
end
if(nargin<4)
    l = 1;
end
if(nargin<5)
    v = 120;
end
if(nargin<6)
    inhibit_f = 0.2;
end
R = l*rand(N,3); % random locations of neurons (in between 0 & 1)
% if(1) % construct radial  distribution instead
%     R(:,1) = (randn(N,1)+10)/10; % radius
%     R(:,2) = rand(N,1)*2*pi;% theta
%     R(:,3) = randn(N,1)/100; % phi
%     [R(:,1) R(:,2) R(:,3)]  = sph2cart(R(:,2),R(:,3),R(:,1));
% end
E = ones(N,1); % E : excitatory neurons
E(rand(N,1)<inhibit_f) = -1;
%% Generate network
% X,Xn containing the node mapping
% Tau containing time delays from X->Xn
% W containing  weights from X->Xn
n = 1:N;
Dx = R(:,1)-R(:,1)';
Dy = R(:,2)-R(:,2)';
Dz = R(:,3)-R(:,3)';
X = (Dx.^2+Dy.^2+Dz.^2);
clear Dx Dy Dz
%X(X>Deff.^2) = 0;
X(rand(N,N)>exp(-X/2/Deff^2)) = 0; % make random connections locally
P = (E+1)+(E'+1); % P is the possibility of connection (1 or 0)
X = X.*P;
clear P;
%X = sparse(X);
T = X*l/v;
T = reshape(T,[1 numel(T)])';
T = T(T~=0);

X = logical(X);
fan_in = sum(X,1);
fan_out = sum(X,2);
W = w0*logical(X)*diag(E);
w_in = sum(W,1)/N;
w_out = sum(W,2)/N;
W = reshape(W,[1 numel(W)])';
W = W(W~=0);


Xn = X*diag(n);
Xn = reshape(Xn,[1 numel(X)])';
Xn = Xn(Xn~=0);

X = diag(n)*X;
X = reshape(X,[1 numel(X)])';
X = X(X~=0);


if(nargout==0 || nargin>6)
    figure('name','Neuron Positions');
    subplot(3,3,[1 2 4 5]);
    %whitebg('k')
    scatter3(R(:,1),R(:,2),R(:,3),10,R); hold on;% R is for color
    h=plot3([R(Xn,1)'; R(X,1)'],[R(Xn,2)' ;R(X,2)'],[R(Xn,3)'; R(X,3)'],'-');
    
%     for i = 1:N % plot synapses individually
%         RXi = R(X{i},:);
%         Ri = R(i,:);
%         Ri = ones(length(RXi(:,1)),1)*Ri;
%         C = [240 180 200]/255;
%         plot3([Ri(:,1) RXi(:,1)],[Ri(:,2) RXi(:,2)],[Ri(:,3) RXi(:,3)],'Color',C);
%     end
    subplot(3,3,3);
    hist(T,1000);
    title('Delays'); ylabel('#Synapses'); xlabel('{\tau_{delay}}(s)');
    subplot(3,3,7);     hist(fan_in,1000); 
    ylabel('#Neurons'); xlabel('Fan-in');
    subplot(3,3,8);    hist(fan_out,1000); 
    ylabel('#Neurons'); xlabel('Fan-out');
    subplot(3,3,9);
    [N_W,W_edge] = histcounts(w_in,500);
    bar(W_edge(2:end),log10(N_W));
    title('Mean W_{in}'); ylabel('#Neurons'); xlabel('Weight');
    subplot(3,3,6);
    [N_W,W_edge] = histcounts(w_out,500);
    bar(W_edge(2:end),log10(N_W));
    title('Mean W_{out}'); ylabel('#Neurons'); xlabel('Weight');
    drawnow;
    hold off;
end
end


