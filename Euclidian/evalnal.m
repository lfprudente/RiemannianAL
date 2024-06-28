function [nal] = evalnal(x,varargin)

global evalfun evalgrad nevalnal

nevalnal = nevalnal + 1;

fgopt = varargin{1};

n = fgopt.n;
m = fgopt.m;
equatn = fgopt.e;
rho    = fgopt.r;
lambda = fgopt.l;
scaling = fgopt.scaling;
sc      = fgopt.sc;

[~,gs,flag] = sevalg(n,x,sc,scaling);
evalgrad = evalgrad + 1;

cs = [];
for i = 1:m
    [~,cs(i),flag] = sevalc(n,x,i,sc,scaling);
    evalfun = evalfun + 1;
end
cs = cs';

nal = gs;
for i = 1:m
    if ( equatn(i) == true || lambda(i) + rho * cs(i) > 0 )
        [~,ncs,flag] = sevalnc(n,x,i,sc,scaling);
        evalgrad = evalgrad + 1;

        nal = nal + ( lambda(i) + rho * cs(i) ) * ncs;
    end
end