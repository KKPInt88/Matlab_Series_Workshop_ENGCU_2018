function X_dot = diff_eqn_1(t, X, t_F, F, params)

% M*x_ddot + B*x_dot + K*x = F

% Define State of the System X = [x1, x2]
% x1 = position
% x2 = velocity

% d/dt (x1) = velocity     = x2
% d/dt (x2) = acceleration = x_ddot = 1/M*( F - B*x_dot - K*x );

% Extract State Variables
x1 = X(1,1);
x2 = X(2,1);

% Extract Parameters of the System
M = params.M;
B = params.B;
K = params.K;

% Interpolation F(t) at time t
F_interp = interp1(t_F, F, t);

X_dot(1,1) = x2;
X_dot(2,1) = 1/M*( F_interp - K*x1 - B*x2 );