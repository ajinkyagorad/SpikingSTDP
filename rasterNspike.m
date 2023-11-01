% raster plot with nth spikes of each neuron joined by a line
function rasterNspike(tSpikeLog,tt,E,Nspikes)
[Nres,L] = size(tSpikeLog);
if(nargin<3 || isempty(E))
   E = ones(1,Nres);
end

if(nargin<2 || isempty(tt))
    tt = 1:L;
end
if(nargin<4 || isempty(Nspikes))
    Nspikes = min(sum(tSpikeLog,2));
end
hold on;
for i = 1:Nres
    spiked_id = find(tSpikeLog(i,:));
    if(E(i)>0) color = 'b'; else color = 'r'; end
    plot(tt(spiked_id),i*ones(1,length(spiked_id)),'.','Color',color,'MarkerSize',0.01);
end
tNSpike = zeros(Nres,Nspikes);
for i = 1:Nres
    tNSpike(i,:) = tt(find(tSpikeLog(i,:),Nspikes));
end
for i = 1:Nspikes
    plot(tNSpike(:,i),1:Nres,'-');
end
    
xlabel('t(s)'); ylabel('Neuron #');
hold off;
end