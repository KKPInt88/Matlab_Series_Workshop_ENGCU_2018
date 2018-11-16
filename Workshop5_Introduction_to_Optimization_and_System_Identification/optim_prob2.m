function [f,g] = optim_prob2(x)

% Extract State Variables
x1 = x(1);
x2 = x(2);

% objective function
f = 20*(-x1^2 + 2*x2)^2 + (x1 - 1)^2 + (x2 - 2)^2;

% Gradient
if nargout > 1
    df_dx1 = 40*(-x1^2 + 2*x2)*(-2*x1) + 2*(x1 - 1)*(1) + 0;
    df_dx2 = 40*(-x1^2 + 2*x2)*(2)     + 0              + 2*(x2 - 2)*(1);
    g = [df_dx1;
         df_dx2];
end

