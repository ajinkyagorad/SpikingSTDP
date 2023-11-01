%% Display network using matlab graph
function display_network_dir_graph(X,Xn,R,W);
if(nargin<4 || isempty(W)); W = ones(length(X),1); end
if(issparse(W)); W = full(W); end
G = digraph(X,Xn,W);
if(nargin<4 || isemtpy(W))
    plot(G,'XData',R(:,1),'YData',R(:,2),'ZData',R(:,3));
else
    plot(G,'XData',R(:,1),'YData',R(:,2),'ZData',R(:,3),'EdgeLabel',G.Edges.Weight);
end
end