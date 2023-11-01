% count loops in a directed network
%function loops = get_loops(X,Xn)
%clear all;
%[X,Xn,Tau,W,R,E]=createNetwork([2 2 2]);
function Loop3 =get_loops(X,Xn,N)
%[X,sortid_X] = sort(X);
%Xn = Xn(sortid_X);
%N = length(E);
A = sparse(X,Xn,1,N,N);
Loop3 = zeros(3,0);
for i = 1:N
    %out = Xn(ismember(X,out))';
    a = sparse(zeros(N,N));
    a(Xn(X==i),X(Xn==i)) = 1;
    [node2,node3, ~] = find(a.*A);
    Loop3 = [Loop3 [i*ones(1,length(node2));node2';node3']];
end
%size(Loop3)
%hist(Loop3(1,:),500)
