function g = asa_grad(x,varargin)
% gradient function for the problem defined in driver1.c for ASA v2.2

for i = 1:length(x)
    t = sqrt(i);
    g(i) = exp( x(i) ) - t;
end

if nargin > 1
    opts = varargin{1};
    % Use this to pass in extra information
%     disp(opts)
end