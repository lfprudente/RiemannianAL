function [al] = evalal(x,varargin)

global evalfun nevalal

nevalal = nevalal + 1;

fgopt = varargin{1};

n = fgopt.n;
m = fgopt.m;
equatn = fgopt.e;
rho    = fgopt.r;
lambda = fgopt.l;
scaling = fgopt.scaling;
sc      = fgopt.sc;

[~,fs,flag] = sevalf(n,x,sc,scaling);
evalfun = evalfun + 1;

cs = [];
for i = 1:m
    [~,cs(i),flag] = sevalc(n,x,i,sc,scaling);
    evalfun = evalfun + 1;
end
cs = cs';

al = fs;
for i = 1:m
    if ( equatn(i) == true || lambda(i) + rho * cs(i) > 0 )
        al = al + cs(i) * ( lambda(i) + 0.5 * rho * cs(i) );
    else
        al = al - 0.5 * lambda(i)^2 / rho;
    end
end