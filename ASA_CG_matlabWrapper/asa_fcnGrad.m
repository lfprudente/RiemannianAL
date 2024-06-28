function [f,g] = asa_fcnGrad(x,varargin)
% objective and gradient of objective function,
%   for the problem defined in driver1.c for ASA v2.2

f = 0;
g = zeros(size(x));
for i = 1:length(x)
    t = sqrt(i);
    g(i) = exp( x(i) ) - t;
    f = f + exp(x(i)) - t*x(i);
end

if nargin > 1
    opts = varargin{1};
    % Use this to pass in extra information
%     disp(opts)
end