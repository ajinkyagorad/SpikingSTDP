% modifyNetwork.m
% modify the given input network to a new configuration with gaussian
% probability, by updating positions of neuron

function [X,Xn,T,W,R,E]= modifyNetwork(v,X,Xn,T,W,R,E,w,r0,k0,tau)
N = length(E);
Ndim = length(v);
dR = v.*randn(N,3);
R = R+dR;

% network is modified only interms of connection in this case
% weights (and other) are kept default parameters (not modified)


%% Assign synapses to network
% Pure random/ spatial random/ pure spatial
%XXX : doesn't work in old matlab versions
[~,ver] = version;
if(posixtime(datetime(ver))<posixtime(datetime('March 27, 2017')))
    D = (repmat(R(:,1),1,N)-repmat(R(:,1),1,N)').^2;
    for i = 2:Ndim
        D = D + (repmat(R(:,i),1,N)-repmat(R(:,i),1,N)').^2;
    end
else
   
    D = (R(:,1)-R(:,1)').^2;
    for i = 2:Ndim
        D = D + (R(:,i)-R(:,i)').^2;
    end
end

D = sqrt(D);
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
D(ConnProb<rand(N,N)) =0;  
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

end
