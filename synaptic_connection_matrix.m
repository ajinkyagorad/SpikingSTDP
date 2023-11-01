% create synaptic connection matrix
function A = synaptic_connection_matrix(Xn,N,Nsyn)
Nsyn = length(Xn);
A =logical(sparse(N,Nsyn)); % for finding synaptic current use matrix operation
%construct matrix where corresponding columns(input connections)of a row(single neuron)
for i = 1:N
    A(i,Xn==i) = 1;
end
end