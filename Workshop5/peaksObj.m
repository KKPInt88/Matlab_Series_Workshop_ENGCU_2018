function f = peaksObj(x)
%PEAKSOBJ casts PEAKS function to a form accepted by optimization solvers.
%   PEAKSOBJ(X) calls PEAKS for use as an objective function for an
%   optimization solver.  X must conform to a M x 2 or N x 2 array to be
%   valid input.
%
%   Syntax
%      f = peaksObj(x)
%
%   Example
%      x = [ -3:1:3; -3:1:3]
%      f = peaksObj(x)
%
%   See also peaks

% Check x size to pass data correctly to PEAKS
[m,n] = size(x);

if (m*n) < 2
    error('peaksObj:inputMissing','Not enough inputs');
elseif (m*n) > 2 && (min(m,n) == 1) || (min(m,n) > 2)
    error('peaksObj:inputError','Input must have dimension m x 2');
elseif n ~= 2
    x = x';
end 

% Objective function
f = peaks(x(:,1),x(:,2));