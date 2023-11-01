% simulate N particles brownian motion
close all;
N = 10;
R = rand(N,3); % initial position of particles
dt = 1E-3; % time step
T = 0.2 ; % Total simulation time
sim_fig = figure('name','simulation');
xlim([0 1]); ylim([0 1]); zlim([0 1]);
for t = 0:dt:T
    
    % display particle positions
    figure(sim_fig);
    scatter3(R(:,1)',R(:,2)',R(:,3)');
    drawnow;
    % compute new
    v = 10*randn(N,3); % velocity vector
    R = R+v*dt;   
end
%% Simulate with trajectory
% http://lewysjones.com/how-to-create-view-3d-graphs-in-matlab-with-3d-glasses/
% for 3d vision
clear all;
N = 100;
R = rand(N,3); % initial position of particles
dt = 1E-3; % time step
T = 1 ; % Total simulation time
sim_fig = figure('name','simulation');
xlim([0 1]); ylim([0 1]); zlim([0 1]);
for t = 0:dt:T
    % display particle positions
    figure(sim_fig);
    plot3(squeeze(R(:,1,:))',squeeze(R(:,2,:))',squeeze(R(:,3,:))');
    drawnow;
    % compute new
    v = 10*randn(N,3); % velocity vector
    R(:,:,end+1) = R(:,:,end)+v*dt;
end