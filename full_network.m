% add the output layer neurons to the network created by createNetwork
% each output neuron is fully connected to the reservoir (input network)
% Pass the paramters created, with additional parameters 
% X,Xn,T,W,R,E are generated by createNetwork
% Nout is the number of added neurons
% final size of each variable increases by N*Nout along its length
function [X_f,Xn_f,T_f,W_f,R_f,E_f]= full_network(X,Xn,T,W,R,E,Nout)
    if(nargin<7)
    Nout = 10; %no. of neurons in the output layer
    end
        
    N = length(E);
    %index of new added neurons = 136 to 135+Nout
    %generating matrices to be appended to X and Xn
    a = 1:N;
    a = reshape((repmat(a',[1,Nout])'),[Nout*N,1]);
    b = (N+1:N+Nout)';
    b = repmat(b,[N,1]);

    %extending T (Nsynx1)
    T_f = [T;T(1)*ones(Nout*N,1)]; % extend by assuming time of first synapse

    %extending E (Nx1)
    E_f = [E; 2*ones(Nout,1)]; % all neurons are excitatory, added neurons are with id 2

    %extending R (Nx3)
    R_f = [R;[-1*ones(Nout,1) -1*ones(Nout,1) (1:Nout)']]; % extend on the vertical line in space

    %extending weight ( Nsynx1)
    W0 =1; % parameters of the weights to be extended
    Wo = W0*ones(Nout*N,1);
    W_f = [W;Wo];

    %extending Xn and X (Nsynx1 each)
    X_f = [X;a];
    Xn_f = [Xn;b];

%     [X_f,sort_index] = sort(X_f); % make them proper
%     Xn_f = Xn_f(sort_index);
%     W_f = W_f(sort_index);
    %XXX order T if variable!
    [Xn_f,sort_index] = sort(Xn_f); % make them proper
    X_f = X_f(sort_index);
    W_f = W_f(sort_index);
    %XXX order T if variable!
end