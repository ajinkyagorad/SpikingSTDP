%% Simulate izhikevich model neuron 
% for parameters a,b,c,d
% parameters a,b are selected from the graph GUI
% implement equations
%   v' = 0.04v^2+5v+140-u+I
%   u' = a(bv-u)
% if v = 30mv
%   v<-c, u<- u+d

% simulate  single  neuron with constant current input
figure('name','Izhikevich');
set(gcf,'Color','w');
ab=subplot(2,2,1); title('a,b'); xlabel('a'); ylabel('b'); xlim([0 0.12]); ylim([0 0.35]);
hold on;
text(0.02,0.25,'LTS,TC'); 
text(0.02,0.2,'RS,IB,CH');
text(0.1,0.25,'RZ');
text(0.1,0.2,'o');
plot(0.02,0.25,'o'); 
plot(0.02,0.2,'o');
plot(0.1,0.25,'o');
plot(0.1,0.2,'o');

cd=subplot(2,2,2); title('c,d'); xlabel('c'); ylabel('d'); xlim([-70 -40]); ylim([0 9]);
hold on;
text(-65,8,'RS');
text(-55,4,'IB');
text(-50,2,'CH');
text(-65,2,'FS,LTS,RZ');
text(-65,0.05,'TC');
plot(-65,8,'o');
plot(-55,4,'o');
plot(-50,2,'o');
plot(-65,2,'o');
plot(-65,0.05,'o');
response = subplot(2,2,[3 4]); title('Neuron response'); xlabel('time(ms)'); ylabel('V(mV)');
while(1)
    subplot(ab); fprintf('Select (a,b)...,');
    [a,b]= ginput(1);
    subplot(cd);fprintf('Select (c,d)\r\n');
    [c,d] = ginput(1);
    v = -65;
    u = b.*v;
    I = 20;
    V = [];
    U = [];
    K = 500;
    for k = 1:K
    if(v>30)
        V(k) = 100;
    else 
        V(k) = v;
    end
    U(k) = u;
    u(v>30) = u+d;
    v(v>30) = c;
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
    u=u+a.*(b.*v-u);                 % stability
    end
    subplot(response);
    plot(V); hold on; plot(U); hold off;
    title(['(a,b,c,d)=(' num2str(a) ',' num2str(b) ',' num2str(c) ',' num2str(d) ')']);
end