% Display a dynamic network with randomly moving nodes(neurons),
% where connectivity is gaussian probabilistic
%                   createNetwork(resSize,w,r0,k0,f_inhibit,tau,show)
% motivation  :
% https://www.jove.com/video/2056/how-to-culture-record-stimulate-neuronal-networks-on-micro-electrode
RECORD_VIDEO = 1;
%%
if(RECORD_VIDEO)
    writerObj = VideoWriter([int2str(now*1E6)],'MPEG-4'); % Name it.
    fprintf('Video writing to file %s\r\n',[int2str(now) '.avi']);
    open(writerObj);
end
lmax =4;
[X,Xn,Tau,W,R,E]=createNetwork([4,4,4],[3 6; -2 -2],2,[0.45 0.3;0.6 0.3],0.2,1E-3);

N = length(E);
sim_fig= figure('name','Simulation');
G = digraph(sparse(X,Xn,1,N,N));
h =plot(G,'XData',R(:,1),'YData',R(:,2),'ZData',R(:,3),'EdgeAlpha',0.2); 
%h =plot(G,'Layout','layered','Sinks',40,'EdgeAlpha',0.2); hold on;
%https://in.mathworks.com/help/matlab/ref/graph.plot.html
bgColor = [0.01 0.01 0.01];
set(gca,'Color',bgColor);set(gcf,'color',bgColor);
set(0,'DefaultFigureWindowStyle','normal');
set(gcf, 'Position', get(0, 'Screensize'));
%set(gcf,'units','normalized','outerposition',[0 0 1 1]);
hold on;
%light('Position',[-1 0 0],'Style','local')
axis equal
camproj('perspective')
%
for k = 1:1000
    tic
    if(mod(k,1)==0)
       %display_network(X,Xn,Tau,W,R,E,N,sim_plot);
       G = digraph(sparse(X,Xn,1,N,N)); hold off;
        figure(sim_fig);
        h =plot(G,'XData',R(:,1),'YData',R(:,2),'ZData',R(:,3),'EdgeAlpha',0.2);
        
        xlim([0 lmax]); ylim([0 lmax]); zlim([0 lmax]);
        set(gca,'color',bgColor);
        %h =plot(G,'Layout','force3','Iterations',10,'EdgeAlpha',0.2);

        if(RECORD_VIDEO)
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
        end
        drawnow;
    end
    R = R+randn(N,3)*0.05;  
    R = mod(R,5);
    
    [X,Xn,Tau] = modifyNetwork2(X,Xn,R,[1 2]);
    toc
    % plot histogram of input connections
end

%% MOVIE END
if(RECORD_VIDEO)
    close(writerObj);
end
