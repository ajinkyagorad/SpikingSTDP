% display directional network
% s (X) as source nodes
% t (Xn) as target nodes
% e as type of nodes (excitatory/ inhibitory)
function display_network_dir(X,Xn,E)
if(nargin<3||isempty(E))
    N = max(max(X,Xn));
    E = ones(1,N);
else
    N = length(E);
end
Nsyn = length(X);
nodes_pos_y = 1:N;
nodes_pos_x = ones(N,1);
%ColorDef = zeros(N,3); % color coding of neurons
%ColorDef(find(E>0),:) = ones(length(find(E>0)),3).*repmat([0 0 1],length(find(E>0)),1); % excitatory neuronsosn coded as blue
%ColorDef(find(E<0),:) = ones(length(find(E<0)),3).*repmat([1 0 0],length(find(E<0)),1); % inhibitory neurons coded re
figure('name','digraph-oneway');
hold on;
plot(nodes_pos_x(E>0),nodes_pos_y(E>0),'o','MarkerFaceColor','b');
plot(nodes_pos_x(E<0),nodes_pos_y(E<0),'o','MarkerFaceColor','r');
plot(2*nodes_pos_x(E>0),nodes_pos_y(E>0),'o','MarkerFaceColor','b');
plot(2*nodes_pos_x(E<0),nodes_pos_y(E<0),'o','MarkerFaceColor','r');
plot([X Xn]');
% add interactive functionality perhaps!

