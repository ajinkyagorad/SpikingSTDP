%% Engineer the parameters to be used for the network
% input is the network with connections X & Xn
% return weights such that the given node neuron fires with a threshold
% fraction
function W = engineerWeights(X,Xn)
% steps
% look for effective fan in and put the weights after knowing the threshold
% and post synaptic potential function
%
% weights will be at least such that every neuron has the possibility of
% firing

% for now put a simple rule