%% Generate network for problem 2
function [X,Xn,Tau,W,R,E,N] = p2Network(N,f_inhibit,we,wi,num_syn)
if(nargin<1)
    N = 500;
end
if(nargin<2)
    f_inhibit = 0.2;
end
if(nargin<3)
    we = 3000;
end
if(nargin<4)
    wi = -we;
end
if(nargin<5)
    num_syn = N/10;
end
Nsyn =N*num_syn; % each neuron connects to N/10 synapses
neu_ex = floor((1-f_inhibit)*N);
neu_in = N-neu_ex;
%first synapses are excitatory, can connect to any neuron with id = 1 to N
X = [];
for i = 1:N
    if(i<=neu_ex)
        X = [X ;randperm(N,num_syn)'];
    else
        X = [X ;randperm(neu_ex,num_syn)'];
    end
end
Xn = reshape(ones(num_syn,N)*diag([1:N]),1,N*num_syn)';
Tau =[randi([1 20],1,neu_ex*num_syn) ones(1,neu_in*num_syn)]'*1E-3;
%W = [ we*ones(1,neu_ex*num_syn) wi*ones(1,neu_in*num_syn)]';
W = we*ones(1,Nsyn)';
W(ismember(X,[neu_ex+1:N])) = wi;
R = rand(N,3);
E = [ones(1,neu_ex) -1*ones(1,neu_in)]';

end
