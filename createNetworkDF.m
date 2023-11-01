%% Create multidimensional directional reservoir network for Liquid State Machine (LSM)
% pure random/spatial random/ pure spatial
% resSize (reservoir size) : 1x3 dimensional array with entries as k l m
% create a multidimensional network
% w (weight) : is 2x2 matrix with fixed weight for entries {E,I}->{E,I} [E->E E->I;I->E I->I];
function [X,Xn,T,W,R,E]= createNetworkDF(resSize,direction,w,r0,k0,f_inhibit,tau,show)
if(nargin<1) resSize = [5 5 5]; end
if(nargin<2) direction = [1 0 0]; end
if(nargin<3) w = [3 6 ;-2 -2]; end
if(nargin<4) r0 = 2; end
if(nargin<5) k0 = [0.45 0.3;0.6 0.15]; end
if(nargin<6) f_inhibit = 0.2; end
if(nargin<7) tau = 1E-3; end
if(nargin<8) show = 0; end
%% generate Nd coordinates of neurons
Ndim = length(resSize);
N = prod(resSize); % number of neurons
R = [1:resSize(1)]';
for i_dim = 2:Ndim
    dim = resSize(i_dim);
    R = repmat(R,dim,1);
    r =[1:dim]'*ones(1,prod(resSize(1:i_dim-1)));
    r = reshape(r',length(r(:)),1);
    R(:,i_dim) = r;
end
%% Assign excitatory/inhibitory behaviour
E = ones(N,1); % E : excitatory (+1)
E(rand(N,1)<f_inhibit) = -1; %  inhibitory(-1)
%% Assign synapses to network
% Pure random/ spatial random/ pure spatial
%XXX : doesn't work in old matlab versions
[~,ver] = version;
if(posixtime(datetime(ver))<posixtime(datetime('March 27, 2017')))
    D = direction(1)*(repmat(R(:,1),1,N)-repmat(R(:,1),1,N)');
    for i = 2:Ndim
        D = D + direction(i)*(repmat(R(:,i),1,N)-repmat(R(:,i),1,N)');
    end
else
    
    D = direction(1)*(R(:,1)-R(:,1)');
    for i = 2:Ndim
        D = D + direction(i)*(R(:,i)-R(:,i)');
    end
end
D = -D; % correct the orientation
D(D<0) = Inf;
%%
nE = sum(E>0);
nI = sum(E<0);
% assign on basis of E/I information
% sort distances based on E
[~,sort_E] =sort(E);
[~,sort_back] = sort(sort_E);
D = D(sort_E,sort_E);
% XXX : (optimise following for performance & memory)
ConnProb = [k0(2,2)*exp(-D(E<0,E<0)/r0^2) k0(2,1)*exp(-D(E<0,E>0)/r0^2);...
    k0(1,2)*exp(-D(E>0,E<0)/r0^2) k0(1,1)*exp(-D(E>0,E>0)/r0^2)];
D(D>r0) =0;
clear ConnProb; % XXX : delete all NxN matrices if not needed
D = D(sort_back,sort_back);
%% remove reflected coefficients by choosing one random connection
DL = logical(D);
%DL = (DL);
[r,c] = find(triu(DL.*DL'));
reflected_loops = length(r); % reflected loops in the network
if(reflected_loops>0)
    loop_id = randperm(reflected_loops);
    partition_id = ceil(reflected_loops/2); % create halfway partition for lower and upper triangular matrix
    upper_loops = loop_id(1:partition_id);
    lower_loops = loop_id(partition_id+1:end);
    %{
subplot(121); imagesc(D); hold on;
plot(r(upper_loops),c(upper_loops),'o');
plot(c(lower_loops),r(lower_loops),'o');
subplot(122);
imagesc(logical(D.*D')); hold on;
plot(r(upper_loops),c(upper_loops),'o');
%plot(r(lower_loops),c(lower_loops),'o');
%plot(c(upper_loops),r(upper_loops),'o');
plot(c(lower_loops),r(lower_loops),'o');
    %}
    %remove_ids = sparse([r(upper_loops);c(lower_loops)],[c(upper_loops);r(lower_loops)],1);
    %D(remove_ids==1)= 0;
    for i=upper_loops
        D(r(i),c(i)) = 0;
    end
    for i = lower_loops
        D(c(i),r(i)) = 0;
    end
end
if(find(D.*D'))
    warning('Network has reflected loops');
end
%D(logical(tril(D.*D'))) = 0; % remove all remaining reflected coefficients
%% Assign weights
X = logical(D);

%% Assign synaptic time delay
if(tau~=0)
    T = tau*logical(X);
else
    T = D*1E-3;
end
T = reshape(T,[1 numel(T)])';
T = T(T~=0);
clear D;
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

% X = full(X);
% Xn = full(Xn);
% W = full(W);
% T = full(T);

% Xij = i->j
%% remove connection which form closed loop of 2 neurons
%{
remove_ids= [];
scramble = randperm(length(X));
X = X(scramble); % scramble them so removal is random
Xn = Xn(scramble);
for i = 1:length(X)
    if(Xn(X==Xn(i))==X(i))
        remove_ids(end+1) = i;
        X(i) = 0; % mark the processed node as zero
        Xn(i) = 0;
    end
end
%%
X(remove_ids) = [];
Xn(remove_ids) = [];
[Xn,sort_id] = sort(Xn); % sort them back according to order
X = X(sort_id);
W(remove_ids) = [];
T(remove_ids) = [];
%}
%% Display out information
fprintf('Reservoir with size %s created %i/%i \r\n',mat2str(resSize),length(X),length(E));

if(nargout<1 || show >=1)
    % Color coding
    Color = zeros(N,3); % color coding of neurons
    Color(find(E>0),:) = ones(length(find(E>0)),1).*[0 0 1]; % excitatory neuronsosn coded as blue
    Color(find(E<0),:) = ones(length(find(E<0)),1).*[1 0 0]; % inhibitory neurons coded red
    % Display Network
    figure('name','Neuron Positions');
    subplot(3,3,[1 2 4 5]);
    scatter3(R(:,1),R(:,2),R(:,3),10,Color); hold on;% R is for color
    if(length(X)<10000)
        h=plot3([R(Xn,1)'; R(X,1)'],[R(Xn,2)' ;R(X,2)'],[R(Xn,3)'; R(X,3)'],'-');
    else
        fprintf('Too many synapses to display\r\n');
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
% figure('name','Directional graph');
% display_network_dir_graph(X,Xn,R,W);

end