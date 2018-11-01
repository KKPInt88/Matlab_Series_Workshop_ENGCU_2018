function F = sys_NL_2(x)

% exp(-exp(-x1 + x2))     - x2*(1 + x1^2) = 0
% x1*cos(x2) + x2*sin(x1) - 1/2           = 0
x1 = x(1,1);
x2 = x(2,1);


F(1,1) = 2*x1 - x2   - exp(-x1); 
F(2,1) =  -x1 + 2*x2 - exp(-x2);




