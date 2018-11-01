function F = sys_NL_1(x)

% exp(-exp(-x1 + x2))     - x2*(1 + x1^2) = 0
% x1*cos(x2) + x2*sin(x1) - 1/2           = 0

F(1,1) = exp(-exp(-(x(1)+x(2)))) - x(2)*(1+x(1)^2);
F(2,1) = x(1)*cos(x(2)) + x(2)*sin(x(1)) - 0.5;




