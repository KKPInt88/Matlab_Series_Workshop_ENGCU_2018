%% Solve Linear/Nonlinear System of Equation and Differential Equation
% Kan Kanjanapas (Ph.D.)
% Fri Nov 2, 2018

clc;
close all;
clear all;



%% Part 0: Rank and Solving System of Linear Equations

% Reference: The Rank of a Matrix
% Francis J. Narcowich, Department of Mathermatics, Texas A&M University,
% Jan 2015

% Theorem: Consider the system Ax = b, with coefficient matrix A and
% augmented matrix [A|b]. As above, the sizes of 
%
% dim(b) = m x 1,
% dim(A) = m x n,
% dim([A|b]) = m x (n+1), respectively; 
%
% in addition, the number of unknowns is n. Therefore:
%
% 1. Ax = b is inconsistent (i.e. no solution exists) if and only if
% rank(A) < rank([A|b])
%
% 2. Ax = b has a unique solution if and only if rank(A) = rank([A|b]) = n
% and
% 3. Ax = b has infinite many solutions if and only if rank(A) = rank([A|b]) < n


% Let's take a look at these systems of linear quations

%  System 1,         System 2,            System 3
%  x1 + x2 = 2      3*x1 + 2*x2 = 3      3*x1 + 2*x2 =  3 
%  x1 - x2 = 0     -6*x1 - 4*x2 = 0     -6*x1 - 4*x2 = -6

% Solution
% x1 = 1, x2 = 1,     no solution,        infintie solution, i.e. x1=0, x2=3/2


% [1  1][x1] = [ 2 ]     [ 3  2 ][x1] = [3]     [ 3  2 ][x1] = [ 3]
% [1 -1][x2]   [ 0 ]     [-6 -4 ][x2]   [0]     [-6 -4 ][x2] = [-6]
%   A01 * X01 = b01        A02 * X02 = b02       A03 * X03 = b03

A01 = [1 1; 1 -1];
b01 = [2; 0];
rank(A01)

A02 = [3 2; -6 -4];
b02 = [3; 0];
rank(A02)

A03 = [3 2; -6 -4];
b03 = [3; -6];
rank(A03)

% Augmented matrix [A|b]
% [1  1 2]        [ 3  2 3]        [ 3  2  3]
% [1 -1 0]        [-6 -4 0]        [-6 -4 -6]
Ab01 = [A01 b01];
Ab02 = [A02 b02];
Ab03 = [A03 b03];

rank(Ab01)
rank(Ab02)
rank(Ab03)


% Summary     rank(A)     rank([A|b])    n     # of solutions
% sys1        2           2              2       1
% sys2        1           2              2       0 (no solution)
% sys3        1           1              2       infinite




%% Part 1: Solve System of Linear Equations

% Ex 1.1
% x - y = 1
% x + y = 2
%
% [1 -1][ x ] = [ 1 ]
% [1  1][ y ]   [ 2 ]
%   A1     X       B1  --> A1*X = B1

A1 = [1 -1;
      1  1];
B1 = [1;
      2];
 
% Solve for X 
X1_1 = A1\B1;  % --> Solve A*X = B

% X =
% 
%     1.5000
%     0.5000

% Verify
A1*X1_1;

% or if A is nonsingular, det(A) is non-zero or full rank
% check full rank A
rank(A1); % or 
det(A1);

X1_2 = inv(A1)*B1;

% Verify
A1*X1_2;


%% 1.2) Non square matrix A
A2 = [1 2  3  4;   % 3x4
      5 6  7  8;
      9 10 11 12];
  
B2 = [13;  %3x1
      14;
      15];
  
% Suppose we want to solve:  A2   * X2   = B2 
%                            (3x4) (?)    (3x1) 
rank(A2)
rank([A2 B2])

% Method 1
X2 = A2\B2;  %[-4 0 0 4.25]'
% A2*X2


%% Method 2


% The dimension of A2 matrix is 3x4 --> No direct form of inverse, can't use inv(A2)
% A2*X2 = B2
% A2'*(A2*X2) = A2'*(B2)
% (A2'*A2)*X2 = (A2'*B2)
% Hopefully  X2 = inv(A2'*A2) * (A2'*B2)

% check rank or det of (A2'*A2)
rank(A2'*A2)  % =2  % <3, not full rank --> inv(A2'*A2) = [Inf ....]

[U,S,V] = svd(A2,'econ'); % A2 = U*S*V'
% U*S*V'*X = B2 --> V*S+*U'*B2

U*S*V' % --> A2

format long

% %X2_LS = VR * inv(S) * UR' * B2
% %U*S*V'*V*inv(S)*U'*B2
% %A2        X2_LS
% X2_LS     = V*inv(S)*U'*B2;
% A2*X2_LS
% 
% %X = linsolve(A,B)
% (U*(S*(V'*V)*inv(S))*U')*B2

X2_Pinv = pinv(A2)*B2;
A2*X2_Pinv % = B2?

V*inv(S)*U'*B2 % = X2_Pinv ?


% Reference: Matlab Help: pinv
%U(:,1:2)*(S(1:2,1:2))*V(:,1:2)'
U1 = U(:,1:2);
S1 = (S(1:2,1:2));
V1 = V(:,1:2);

X2_Pinv2 = V1*inv(S1)*U1'*B2;  % **************************g
A2*X2_Pinv2;


%%
A3 = reshape( [1:1000], 100, 10 );
size(A3);
B3 = [201:300]';

rank(A3) % == 2
rank([A3 B3]) % == 2

X3 = A3\B3;
% X3 =
% 
%    0.777777777777778
%                    0
%                    0
%                    0
%                    0
%                    0
%                    0
%                    0
%                    0
%    0.222222222222222

% check 
A3*X3 - B3

% Or use pinv(.)
X3_2 = pinv(A3)*B3
A3*X3_2 - B3


%% 1.3 linsolve(A,B)

% Reference Matlab Help linsolve(A,B)
X3_3 = linsolve(A3,B3);  % --> equivalent to A3\B3





%% Part 2: Solve System of Nonlinear Equation 

% The system of nonlinear equation is given by:
% exp(-exp(-x1 + x2))     = x2*(1 + x1^2)
% x1*cos(x2) + x2*sin(x1) = 1/2

% Convert the equations to the form F(X) = 0
%
% exp(-exp(-x1 + x2))     - x2*(1 + x1^2) = 0
% x1*cos(x2) + x2*sin(x1) - 1/2           = 0

% Write a function that describes F(X)

fun = @sys_NL_1;
X0  = [0; 0];  % Inital Condition
X_sol = fsolve(fun, X0);
X_sol

% Check answer, plug X_sol into fun
sys_NL_1(X_sol)



%% Add on optimization setting in fsolve
% [x, fval, exitflag, output, jacobian] = fsolve(fun, x0, options)

% Optimization Setting
options = optimoptions('fsolve', 'MaxIterations', 100, ...
                                 'Display', 'iter', ...
                                 'PlotFcn', {@optimplotfval, @optimplotx, @optimplotfirstorderopt} );
     
X0 = [0;0];
[X_sol_opt, fval, exitflag, output, jacobian] = fsolve(fun, X0, options)
% Help First-Order Optimality Measure, exitflag                             
% Jacobian at Optimal Point

% % % X_sol_opt =
% % % 
% % %    0.353246561920553
% % %    0.606082026502285
% % % 
% % % 
% % % fval =
% % % 
% % %    1.0e-06 *
% % % 
% % %   -0.240688008967815
% % %   -0.038253165413060
% % % 
% % % 
% % % exitflag =
% % % 
% % %      1
% % % 
% % % 
% % % output = 
% % % 
% % %   struct with fields:
% % % 
% % %        iterations: 4
% % %         funcCount: 15
% % %         algorithm: 'trust-region-dogleg'
% % %     firstorderopt: 2.023187897629516e-07
% % %           message: 'Equation solved.??fsolve completed because the vector of function values is near zero?as measured by the default value of the function tolerance, and?the problem appears regular as measured by the gradient.??Stopping criteria details:??Equation solved. The sum of squared function values, r = 5.939402e-14, is less than?sqrt(options.FunctionTolerance) = 1.000000e-03. The relative norm of the gradient of r,?2.023188e-07, is less than options.OptimalityTolerance = 1.000000e-06.??Optimization Metric                                         Options?relative norm(grad r) =   2.02e-07              OptimalityTolerance =   1e-06 (default)?r =   5.94e-14                              sqrt(FunctionTolerance) =  1.0e-03 (default)'
% % % 
% % % 
% % % jacobian =
% % % 
% % %   -0.166995346546173  -0.863585688173771
% % %    1.390545405447483   0.144718222320080



%% Solve a Problem using Structure Form

problem = [];


% Recall Options
% options = optimoptions('fsolve', 'MaxIterations', 100, ...
%                                  'Display', 'iter', ...
%                                  'PlotFcn', {@optimplotfval, @optimplotx, @optimplotfirstorderopt} );
problem.options     = options;
problem.objective   = @sys_NL_1;  % = fun
problem.x0          = X0;
problem.solver      = 'fsolve';

% Then solve for solution
x_sol_struc = fsolve(problem);



%% Ex2: solve
% 2*x1 - x2  = exp(-x1);
% -x1 + 2*x2 = exp(-x2);

% rewrite F(X) = 0
% 2*x1 - x2   - exp(-x1) = 0
%  -x1 + 2*x2 - exp(-x2) = 0;
 
% options; 
clc;
X0 = [0;0];
[X_sol_opt, fval, exitflag, output, jacobian] = fsolve( @(X)sys_NL_2(X), X0, options)

%check solution
sys_NL_2(X_sol_opt);

fval;
exitflag;
output;
jacobian;




%% Ex3: Matrix Equation

% Find X*X*X = [1 2;
%               3 4];
      
clc;
objective_fun = @(X) X*X*X - [1 2; 3 4];
options;

%X0 = zeros(2,2);
X0 = ones(2,2);

[X_sol_opt, fval, exitflag, output, jacobian] = fsolve(objective_fun, X0, options);
sum(sum(fval.*fval))  





%% Part 3: Solve Differential Equation --> ode45 (Linear/Nonlinear), lsim (Linear)

% Sol 1:
% Ex second order system
% Mass   M = 10 Kg
% Spring K = 100 N/m
% Damper B = 10  N.m/s

% M*x_ddot + B*x_dot + K*x = F(t)

Ts = 10^-3;         % [s]
t  = [0:Ts:10]';    % Time vector 0 to 10 second

F = 1*sin(2*pi*1*t) + 0;

% Plot to check F(t)
%figure,plot(t,F,'b', 'LineWidth', 2); xlabel('Time [s]'); ylabel('Force [N]'); grid on;

% Calculate force response of the given 2nd order system

M = 10;
B = 10;
K = 100;

params = [];
params.M = M;
params.B = B;
params.K = K;

t_F = t;
fun = @(t,X)diff_eqn_1(t, X, t_F, F, params);
X0 = [0;0];
[t_out, Xout] = ode45(fun, t, X0)

x1 = Xout(:,1);
x2 = Xout(:,2);

figure;
subplot(3,1,1);
plot(t, F); ylabel('Force [N]'); grid on;
subplot(3,1,2); plot(t, x1); ylabel('Position [m]'); grid on;
subplot(3,1,3); plot(t, x2); xlabel('Time [s]'); ylabel('Velocity [m/s]'); grid on;




%% Sol 2 --> lsim (Only for linear system) --> Further Read Control System
clc;

% Transform Differential equation to transfer function or state space form

% M*x_ddot + B*x_dot + K*x = F
% (M*s^2 + B*s + K)*X(s) = F(s)
% G(s) = X(s)/F(s) = 1/(M*s^2 + B*s + K);

s = tf('s');  % s variable in s domain (applying Laplace Transform)
G1 = 1/(M*s^2 + B*s + K);

X0 = [0;0];

% Matlab Help lsim
% [y,t,x] = lsim(sys,u,t,x0,method)
[Yout2, tout2, Xout2] = lsim( G1, F, t, X0);

figure;
plot(t, x1, 'b', 'LineWidth', 2);
hold on;
plot(tout2, Yout2, 'r--', 'LineWidth', 2); xlabel('Time [s]'); ylabel('Posiiton [m]'); grid on;



% If define G2 in state space form
A4 = [ 0    1;
      -K/M -B/M];
  
B4 = [ 0; 1/M ];
C4 = [1 0];
D4 = [0];

G2 = ss(A4, B4, C4, D4);


[Yout3, tout3, Xout3] = lsim( G2, F, t, X0);
figure;
plot(t, x1, 'b', 'LineWidth', 2);
hold on;
plot(tout2, Yout2, 'r--', 'LineWidth', 2); 
plot(tout3, Yout3, 'g-.', 'LineWidth', 2); xlabel('Time [s]'); ylabel('Posiiton [m]'); grid on;





%% Sol 3: Simulink


simin_F.time           = t;
simin_F.signals.values = F;

open('Ex_Solve_DiffEqn_Simulink.slx');
sim('Ex_Solve_DiffEqn_Simulink');

tout4   = simout.Time;
x4      = simout.Data(:,1);
x_dot4  = simout.Data(:,2);
x_ddot4 = simout.Data(:,3);


figure;
subplot(4,1,1); plot(tout4, F, 'LineWidth', 2, 'Color', 'k');             ylabel('Force [N]'); grid on;
subplot(4,1,2); plot(tout4, x4, 'LineWidth', 2, 'Color', 'b');            ylabel('Position [m]'); grid on;
subplot(4,1,3); plot(tout4, x_dot4, 'LineWidth', 2, 'Color', [0 0.5 0]);  ylabel('Velocity [m/s]'); grid on;
subplot(4,1,4); plot(tout4, x_ddot4, 'LineWidth', 2, 'Color', 'r');       ylabel('Acceleration [m/s^2]'); xlabel('Time [s]'); grid on;
















