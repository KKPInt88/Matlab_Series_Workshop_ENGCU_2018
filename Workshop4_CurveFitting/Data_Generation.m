%% Data Generation for Curve Fitting Examples
% Kan Kanjanapas (Ph.D.)

clc;
close all;
clear all;


%% Data set for Linear Regression 

x1  = [0:0.1:10]';
a1  = [1 2];
y1  = a1(1) + a1(2)*x1 + 0.5*randn(size(x1));

figure,plot(x1,y1,'o','Color','b'); xlabel('x1'); ylabel('y1'); title('Dataset #1'); grid on;



%% Data set for Polynomial Fit 

x2 = [0:0.1:10]';
polyfun2 = @(x) 1*x.^2 + 2*x + 3;
y2 = polyfun2(x2) + 3*randn(size(x1));

figure,plot(x2,y2,'o','Color','b'); xlabel('x2'); ylabel('y2'); title('Dataset #2'); grid on;



%% Data set for Arbitrary Function Fit 

% x3 = [0.01:0.01:3]';
% polyfun3 = @(x) sin(2*pi*1*x).*exp(-2*x) + log10(x3);
% y3 = polyfun3(x3) + 0*randn(size(x3));
% 
% figure,plot(x3,y3,'o','Color','b'); xlabel('x3'); ylabel('y3'); title('Dataset #3'); grid on;

Data3 = load('workspace_dataset3.mat');
x3 = Data3.x3;
y3 = Data3.y3;
figure,plot(x3,y3,'o','Color','b'); xlabel('x3'); ylabel('y3'); title('Dataset #3'); grid on;



%% Data set for Surface Fit 

x4_vec = [0 : 0.2 : 10]';
y4_vec = [0 : 0.2 : 10]';

[x4, y4] = meshgrid(x4_vec, y4_vec);
z4 = 1 + 2*x4 + 3*y4 + 4*x4.*y4 + 5*x4.^2 - 6*y4.^2 - 7*x4.*y4;
z4 = z4 + 20*randn(size(z4));

figure;
set(gcf, 'Position', [0 0 2560 1280]/2);
surf(x4, y4, z4, 'FaceAlpha', 0.4, 'Marker', 'o'); 
xlabel('x4'); ylabel('y4'); zlabel('z4'); title('Dataset #4'); grid on; view(-122,24);

