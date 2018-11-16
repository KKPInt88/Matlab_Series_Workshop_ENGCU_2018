clc;
close all;
clear all;

fs = 1000;
Ts = 1/fs;


%% Open Simulink Model
L1 = 20;  % [cm]
L2 = 12;  % [cm]

open('Part2_6_Four_Bar.slx');


%% Update Model
set_param(gcs,'SimulationCommand','Update'); % Update Model


%% Simulate the Simulink Model
sim('Part2_6_Four_Bar.slx');


%% Extract Data

Time = simout.Time;
q    = simout.Data(:,1);
wy   = simout.Data(:,2);
x    = simout.Data(:,3);
y    = simout.Data(:,4);
z    = simout.Data(:,5);

figure;
set(gcf, 'Position', [0 0 2560 1280]/2);

subplot(5,1,1); plot(Time, q, 'LineWidth', 2, 'Color', 'b');
xlabel('Time [s]'); 
ylabel('q [rad]');
grid on;
set(gca, 'FontSize', 14);

subplot(5,1,2); plot(Time, wy, 'LineWidth', 2, 'Color', 'r');
xlabel('Time [s]'); 
ylabel('wy [rad/s]');
grid on;
set(gca, 'FontSize', 14);

subplot(5,1,3); plot(Time, x, 'LineWidth', 2, 'Color', [0 0.5 0]);
xlabel('Time [s]'); 
ylabel('x [m]');
grid on;
set(gca, 'FontSize', 14);

subplot(5,1,4); plot(Time, y, 'LineWidth', 2, 'Color', [1 0 1]);
xlabel('Time [s]'); 
ylabel('y [m]');
grid on;
set(gca, 'FontSize', 14);

subplot(5,1,5); plot(Time, z, 'LineWidth', 2, 'Color', [0 1 1]);
xlabel('Time [s]'); 
ylabel('z [m]');
grid on;
set(gca, 'FontSize', 14);


%% 
figure;
plot(x, z, 'LineWidth', 2, 'Color', 'k');
xlabel('X [cm]'); ylabel('Z [cm]'); grid on;



%% Optimization

Data_opt = load('workspace_opt.mat');
x_ref = Data_opt.x; % L1 = 18, L2 = 12 [cm]
y_ref = Data_opt.y;
z_ref = Data_opt.z;


%% Optimization Part
clc;
close all;

X0 = [17.5   11.2]; 
% lower_bound = [13, 11];
% upper_bound = [19, 13];

lower_bound = [17, 11];
upper_bound = [19, 13];

%close_system('Part2_5_Four_Bar2.slx');

objective = @(X)cost(X, [x_ref z_ref]);
options = optimset('Display','iter','TolX',1e-8, 'MaxIter', 50, 'PlotFcns', {'optimplotx', 'optimplotfval' }); % optimplotfval
% options = optimoptions('patternsearch','PollMethod','GSSPositiveBasis2N', ...
%     'Display','iter','PlotFcn', @psplotbestf,'MaxIterations',40, ...
%     'MeshTolerance',0.001,'UseCompletePoll',true,'InitialMeshSize',mean(x0), ...
%     'MeshContractionFactor',0.25);

mdl = 'Part2_6_Four_Bar.slx';
open(mdl);


figure(101);
set(gcf, 'Position', [1000 200 2560 1280]/2);
plot(x_ref, z_ref, 'LineWidth', 4, 'Color', 'k');
hold on;
xlabel('X [m]');
ylabel('Z [m]');
title('Parametric Design Optimization for Trajectory Following Application');
%legend('Reference Trajectory');
grid on;
set(gca, 'FontSize', 16);

% Solve Optimization Problem
[X_opt, fval] = fmincon(objective, X0, [], [], [], [], lower_bound, upper_bound,[], options)
L1_opt = X_opt(1);
L2_opt = X_opt(2);

% Display Result
disp(sprintf('Optimal Design Parameters\nL1 = %.4f [cm]\nL2 = %.4f [cm]', L1_opt, L2_opt));



