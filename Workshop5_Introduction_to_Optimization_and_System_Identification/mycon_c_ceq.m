function [c,ceq] = mycon_c_ceq(x)

% x1^2 + x^2 <= 3^2
c   = x(1)^2 + x(2)^2 - 9;    % Compute nonlinear inequalities at x.
ceq = [];                     % Compute nonlinear equalities at x.