%% Plot weight vs deltaT of spike Time arrival

% get spike timings and weight change logged data (in LOG and tSpikeLog
% rasterplot)
tspikepre_arrive(tspikepre_arrive>length(LOG.W)) = [];
tspikepost (tspikepre_arrive>length(LOG.W)) = [];
prepostid = [ones(1,length(tspikepre_arrive)) 2*ones(1,length(tspikepost))];
prear_postspikes = [reshape(tspikepre_arrive,1,numel(tspikepre_arrive)) reshape(tspikepost,1,numel(tspikepost))];
[prear_postspikes,sort_id] = sort(prear_postspikes);
prepostid = prepostid(sort_id);
deltaT = diff(prepostid).*diff(prear_postspikes);
deltaW = diff(LOG.W(prear_postspikes));
figure('name','STDP');hold on;
%plot(deltaT,deltaW,'*');

plot(deltaT(logical((deltaW>0) .* (deltaT>0))),deltaW(logical((deltaW>0) .* (deltaT>0))),'*','Color','g');
plot(deltaT(logical((deltaW<0) .* (deltaT>0))),deltaW(logical((deltaW<0) .* (deltaT>0))),'*','Color','y');
plot(deltaT(logical((deltaW>0) .* (deltaT<0))),deltaW(logical((deltaW>0) .* (deltaT<0))),'*','Color','c');
plot(deltaT(logical((deltaW<0) .* (deltaT<0))),deltaW(logical((deltaW<0) .* (deltaT<0))),'*','Color','r');
line([0 0],ylim,'Color','k'); line(xlim,[0 0],'Color','k');
xlabel('{\DeltaT_{post-pre}} (ms)'); ylabel('{\DeltaW}');
%%
