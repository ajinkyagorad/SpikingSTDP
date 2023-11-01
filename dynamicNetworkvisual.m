% Display a dynamic network with randomly moving nodes(neurons),
% where connectivity is gaussian probabilistic
%                   createNetwork(resSize,w,r0,k0,f_inhibit,tau,show)
% motivation  :
% https://www.jove.com/video/2056/how-to-culture-record-stimulate-neuronal-networks-on-micro-electrode 
RECORD_VIDEO = 0;
%%
if(RECORD_VIDEO)
    writerObj = VideoWriter([int2str(now) '.avi'],'MPEG-4'); % Name it.
    fprintf('Video writing to file %s\r\n',[int2str(now) '.avi']);
    open(writerObj);
end

[X,Xn,Tau,W,R,E]=createNetwork([4,4,4],[3 6; -2 -2],2,[0.45 0.3;0.6 0.3],0.2,1E-3);

N = length(E);
figure('name','Simulation');
sim_plot = subplot(121);    title('Network  Graph');
hist_plot = subplot(122);   title('Distance vs Synapses');

for k = 1:1000    
    tic
    if(mod(k,1)==0)
        subplot(sim_plot);
        display_network(X,Xn,Tau,W,R,E,N,sim_plot);
        subplot(hist_plot);
        h = histogram(Tau,40);
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
