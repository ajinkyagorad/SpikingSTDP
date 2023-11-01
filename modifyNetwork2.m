% modifyNetwork.m
% modify the given input network with new positions according to binding
% and breaking thresholds
% Network is represented in compact source-target(X->Xn) notation
%

function [X,Xn,Tau]= modifyNetwork(X,Xn,R,breakBindDist)
if(nargin<4)
    breakingDistance = 1;
    bindingDistance= 0.5;
else
    breakingDistance = breakBindDist(2);
    bindingDistance = breakBindDist(1);
end
[N,Ndim] = size(R);
%% Assign synapses to network
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
activeSyn =logical(sparse(X,Xn,1));
ind = sub2ind(size(D),X,Xn);
brokenSynapses = D(ind)>breakingDistance;
X(brokenSynapses) = [];
Xn(brokenSynapses) = [];
[x,xn,~] = find((D<bindingDistance).*(~activeSyn));
X = [X;x]; Xn = [Xn;xn];
Tau = sparse(X,Xn,1,N,N).*D;
[~,~,Tau] = find(Tau);
end
