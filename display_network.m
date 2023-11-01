% function display network 
function display_network(X,Xn,T,W,R,E,Nout,sim_fig)
if(nargin<7)
    Nout = 0
end
N = length(E);
E(end-Nout+1:end) = 2;
 Color = zeros(N,3); % color coding of neurons
   [~,ver] = version;
    if(posixtime(datetime(ver))<posixtime(datetime('March 27, 2017')))
        Color(find(E==1),:) =repmat([0 0 1],length(find(E==1)),1);% ones(length(find(E>0)),1).*[0 0 1]; % excitatory neuronsosn coded as blue
        Color(find(E<0),:) = repmat([1 0 0],length(find(E<0)),1); % inhibitory neurons coded red
        Color(find(E==2),:) =repmat([0 0 1],length(find(E==2)),1);
    else
        Color(find(E==1),:) = ones(length(find(E==1)),1).*[0 0 1]; % excitatory neuronsosn coded as blue
        Color(find(E<0),:) = ones(length(find(E<0)),1).*[1 0 0]; % inhibitory neurons coded red
        Color(find(E==2),:) = ones(length(find(E==2)),1).*[0 1 0]; % Output neurons coded green
    end
    % Display Network
    if(nargin<8)
        figure('name','Neuron Positions');
    end
    scatter3(R(:,1),R(:,2),R(:,3),20,Color,'filled'); hold on;% R is for color
    if(length(X)<10000)
    h=plot3([R(Xn,1)'; R(X,1)'],[R(Xn,2)' ;R(X,2)'],[R(Xn,3)'; R(X,3)'],'-');
    else
        fprintf('Too many synapses to display');
    end
    
    drawnow;
    hold off;