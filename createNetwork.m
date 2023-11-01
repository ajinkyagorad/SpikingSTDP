%% Create reservoir network for Liquid State Machine (LSM)
% pure random/spatial random/ pure spatial
% resSize (reservoir size) : 1x3 dimensional array with entries as k l m
% for k x l x m sized reservoir (multidimentional networks?)
% w (weight) : is 2x2 matrix with fixed weight for entries {E,I}->{E,I} [E->E E->I;I->E I->I];
function [X,Xn,T,W,R,E]= createNetwork(resSize,w,r0,k0,f_inhibit,tau,show)
if(nargin<1) resSize = [3 3 15]; end
if(nargin<2) w = [3 6 ;-2 -2]; end
if(nargin<3) r0 = 2; end
if(nargin<4) k0 = [0.45 0.3;0.6 0.15]; end
if(nargin<5) f_inhibit = 0.2; end
if(nargin<6) tau = 1E-3; end
if(nargin<7) show = 0; end
%% align input information
K = resSize(1);
L = resSize(2);
M = resSize(3);
N = K*L*M; % number of neurons
%% generate position of neurons
R = zeros(K,L,M,3);
for k = 1:K ;   R(k,:,:,1) = k;end
for l = 1:L ;   R(:,l,:,2) = l;end
for m = 1:M ;   R(:,:,m,3) = m;end
R = reshape(R,[N 3]);
%% Assign excitatory/inhibitory behaviour
E = ones(N,1); % E : excitatory (+1)
E(rand(N,1)<f_inhibit) = -1; %  inhibitory(-1)
%% Assign synapses to network
% Pure random/ spatial random/ pure spatial
%XXX : doesn't work in old matlab versions
[~,ver] = version;
if(posixtime(datetime(ver))<posixtime(datetime('March 27, 2017')))
    Dx = repmat(R(:,1),1,N)-repmat(R(:,1),1,N)';
    Dy = repmat(R(:,2),1,N)-repmat(R(:,2),1,N)';
    Dz = repmat(R(:,3),1,N)-repmat(R(:,3),1,N)';
else
    Dx = R(:,1)-R(:,1)';
    Dy = R(:,2)-R(:,2)';
    Dz = R(:,3)-R(:,3)';
end

X = sqrt(Dx.^2+Dy.^2+Dz.^2);
clear Dx Dy Dz
nE = sum(E>0);
nI = sum(E<0);
% assign on basis of E/I information
% sort distances based on E
[~,sort_E] =sort(E);
[~,sort_back] = sort(sort_E);
X = X(sort_E,sort_E);
% XXX : (optimise following for performance & memory)
ConnProb = [k0(2,2)*exp(-X(E<0,E<0)/r0^2) k0(2,1)*exp(-X(E<0,E>0)/r0^2);...
    k0(1,2)*exp(-X(E>0,E<0)/r0^2) k0(1,1)*exp(-X(E>0,E>0)/r0^2)];
X(ConnProb<rand(N,N)) =0;
clear ConnProb; % XXX : delete all NxN matrices if not needed

X = X(sort_back,sort_back);
%% Assign weights
X = logical(X);
fan_in = sum(X,1); % fan in of neuron
fan_out = sum(X,2); % fan out of neuron
W = sparse(size(X));
W(E>0,E>0) = w(1,1);
W(E>0,E<0) = w(1,2);
W(E<0,E>0) = w(2,1);
W(E<0,E<0) = w(2,2);
W(X==0) = 0;
w_in = sum(W,1); % net input weight to neuron
w_out = sum(W,2); % net output weight of neuron
%% Create arrays out of 2D matrices and reduce them
W = reshape(W,[1 numel(W)])';
W = W(W~=0);
n = 1:N; % id of neuron with  order (((1...K)...L)...M)
% Connections formed are of the form X->Xn
Xn = X*diag(n);
Xn = reshape(Xn,[1 numel(X)])';
Xn = Xn(Xn~=0); % destination neuron

X = diag(n)*X;
X = reshape(X,[1 numel(X)])';
X = X(X~=0);  % input neurons for destination neuron

%% Assign synaptic time delay
if(tau==0)
    T = 1E-3*randi([1 20],size(X));
else
T = tau*ones(size(X));
end
% Xij = i->j
%% Display out information
fprintf('Reservoir with size %ix%ix%i created \r\n',K,L,M);
if(nargout<1 || show >=1)
    % Color coding
    Color = zeros(N,3); % color coding of neurons
    if(posixtime(datetime(ver))<posixtime(datetime('March 27, 2017')))
        Color(find(E>0),:) =repmat([0 0 1],length(find(E>0)),1);% ones(length(find(E>0)),1).*[0 0 1]; % excitatory neuronsosn coded as blue
        Color(find(E<0),:) = repmat([1 0 0],length(find(E<0)),1); % inhibitory neurons coded red
    else
        
        Color(find(E>0),:) = ones(length(find(E>0)),1).*[0 0 1]; % excitatory neuronsosn coded as blue
        Color(find(E<0),:) = ones(length(find(E<0)),1).*[1 0 0]; % inhibitory neurons coded red
    end
    
    % Display Network
    figure('name','Neuron Positions');
    subplot(3,3,[1 2 4 5]);
    scatter3(R(:,1),R(:,2),R(:,3),10,Color); hold on;% R is for color
    if(length(X)<10000)
        h=plot3([R(Xn,1)'; R(X,1)'],[R(Xn,2)' ;R(X,2)'],[R(Xn,3)'; R(X,3)'],'-');
    else
        fprintf('Too many synapses to display');
    end
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