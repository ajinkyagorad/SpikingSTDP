% Raster plot of given 2D spike binary array (with time)
%tt,tSpikeLog
function raster_plot(tSpikeLog,tt,E)
[N,L] = size(tSpikeLog);
if(nargin<3 || isempty(E))
   E = ones(1,N);
end

if(nargin<2 || isempty(tt))
    tt = 1:L;
    xlabel('samples');
    xlim([1 L]);
else 
    xlabel('t(s)');
    xlim([min(tt) max(tt)]);
end
cla;hold on;
for i = 1:N
    spiked_id = find(tSpikeLog(i,:));
    if(E(i)>2) color = 'm'; elseif(E(i)>1) color = 'g'; elseif(E(i)>0) color = 'b'; else color = 'r'; end
    plot(tt(spiked_id),i*ones(1,length(spiked_id)),'.','Color',color,'MarkerSize',5);
end
ylabel('Neuron #');
hold off;
ylim([1 N]);
end