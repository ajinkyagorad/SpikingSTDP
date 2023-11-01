% Show strongly connected components in a graph

function graphSCCdisplay

[X,Xn,Tau,W,R,E]=createNetworkM([10 10 10],[3 6; -2 -2],2,0.1*[0.5 0.5; 0.5 0.5],0.2);

N = length(E);
Nsyn = length(Xn);
G = digraph(sparse(X,Xn,1,N,N));
SCC = tarjan(sparse(X,Xn,1,N,N)');
Loop3 = get_loops(X,Xn,N);
h =plot(G,'XData',R(:,1),'YData',R(:,2),'ZData',R(:,3),'EdgeAlpha',0.2,'ButtonDownFcn',@cbHighlightNode); hold on;
%https://in.mathworks.com/help/matlab/ref/graph.plot.html
%h =plot(G,'Layout','subspace3','Dimension',length(SCC),'EdgeAlpha',0.2,'ButtonDownFcn',@cbHighlightNode); % autoposition plot
bgColor = 1.01-[0.01 0.01 0.01]; 
set(gca,'Color',bgColor);set(gcf,'color',bgColor);
set(0,'DefaultFigureWindowStyle','normal');
%set(gcf, 'Position', get(0, 'Screensize'));
%set(gcf,'units','normalized','outerposition',[0 0 1 1]);

hold on;
%light('Position',[-1 0 0],'Style','local')
axis equal
camproj('perspective')
Colors = distinguishable_colors(length(SCC),'b');

slhan=uicontrol('style','slider','String','Alpha','position',[100 0 200 20],'min',0,'max',0.4,'BackgroundColor',bgColor,'ForegroundColor',1-bgColor);
s2han=uicontrol('style','slider','String','SCC select','position',[100 50 200 20],'min',1,'max',length(SCC)+1,'Value' ,1,'BackgroundColor',bgColor,'ForegroundColor',1-bgColor);
s3han=uicontrol('style','slider','String','Loop3 select','position',[100 90 200 20],'min',1,'max',length(Loop3(1,:))+1,'Value' ,1,'BackgroundColor',bgColor,'ForegroundColor',1-bgColor);

persistChkBox = uicontrol('style','checkbox','String','Persist','position',[10 90 80 20],'Value',0,'BackgroundColor',bgColor,'ForegroundColor',1-bgColor);

chkBoxSCCseparate = uicontrol('style','checkbox','String','Separate SCC','position',[10 70 80 20],'Value',0,'callback',@cbSeparateSCC,'BackgroundColor',bgColor,'ForegroundColor',1-bgColor);
textIn = uicontrol('style','edit','position',[100 25 100 20],'string','NeuronID','callback',@cbHighlightNeuron,'BackgroundColor',bgColor,'ForegroundColor',1-bgColor);
hSliderListener = addlistener(slhan, 'Value','PostSet',@callbackfn);
hSliderListener2 = addlistener(s2han,'Value','PostSet',@cbSelectSCC);
hSliderListener3 = addlistener(s3han,'Value','PostSet',@cbSelectLoop3);
  
    function cbSeparateSCC(source,eventdata)
        isChecked = source.Value;
        R2 = R;
        if(isChecked)
        [step,dimension] = min(max(R));
        for i = 1:length(SCC) % separate functional units into different cubes
            if(length(SCC{i})>4) % could add a button for spiltting it accordingly
                R2(SCC{i},dimension) = R(SCC{i},dimension)+step*(i-1);
            end
        end
        end
        set(h,'XData',R2(:,1),'YData',R2(:,2),'ZData',R2(:,3));
    end
    function cbHighlightNode(source,eventdata)
        dist = (source.XData-eventdata.IntersectionPoint(1)).^2+...
            (source.YData-eventdata.IntersectionPoint(2)).^2+...
            (source.ZData-eventdata.IntersectionPoint(3)).^2;
        if(dist>2)
            return
        end
        [~,neuronID] = min(dist);
        % fprintf("Number : %i\r\n",neuronID);
        set(textIn,'string',num2str(neuronID));
        % highlight all the incoming and outgoing connections and this
        % neuron
        if(~persistChkBox.Value) ;set(h,'NodeColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'LineWidth',0.5,'MarkerSize', 2);end
        highlight(h,neuronID,'MarkerSize',4);
        highlight(h,X(Xn==neuronID),Xn(Xn==neuronID),'EdgeColor','g','LineWidth',3); % green are incoming synapses
        highlight(h,X(X==neuronID),Xn(X==neuronID),'EdgeColor','r','LineWidth',3); % red are outgoing synapses
    end
    function callbackfn(source,eventdata)
        d = get(eventdata.AffectedObject, 'Value');
        set(h,'EdgeAlpha',d);
    end
    function cbSelectLoop3(source,eventdata)
        d = get(eventdata.AffectedObject,'Value');
        d = floor(d);
        nodes = Loop3(:,d);
        if(~persistChkBox.Value);set(h,'NodeColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'LineWidth',0.5,'MarkerSize', 2);end
        highlight(h,nodes,'NodeColor',Colors(2,:),'MarkerSize',4); % could be done in one step also by creating NodeCData
        highlight(h,nodes,[nodes(2:end);nodes(1)],'EdgeColor',Colors(mod(d,length(Colors(:,1)))+1,:),'LineWidth',3); % also highlight input edges to show SCC
    end
    function cbSelectSCC(source,eventdata)
        d = get(eventdata.AffectedObject,'Value');
        d = floor(d);
        nodes = SCC{d};
        if(~persistChkBox.Value);set(h,'NodeColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'LineWidth',0.5,'MarkerSize', 2);end
        highlight(h,nodes,'NodeColor',Colors(2,:),'MarkerSize',4); % could be done in one step also by creating NodeCData
        highlight(h,X(ismember(Xn,nodes)),Xn(ismember(Xn,nodes)),'EdgeColor',Colors(d,:),'LineWidth',3); % also highlight input edges to show SCC
    end
end

%% Temp code
%{
   % for id = 1:length(SCC)
    %nodes = SCC{id};
    %https://in.mathworks.com/help/matlab/ref/matlab.graphics.chart.primitive.graphplot.highlight.html
    highlight(h,nodes,'NodeColor',Colors(id,:)); % could be done in one step also by creating NodeCData
    highlight(h,X(ismember(Xn,nodes)),Xn(ismember(Xn,nodes)),'EdgeColor',Colors(id,:));
end
%}