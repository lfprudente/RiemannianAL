function [al,nal] = evalalnal(x,varargin)

global evalfun evalgrad nevalal nevalnal

nevalal  = nevalal + 1;
nevalnal = nevalnal + 1;

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

[~,gs,flag] = sevalg(n,x,sc,scaling);
evalgrad = evalgrad + 1;

cs = [];
for i = 1:m
    [~,cs(i),flag] = sevalc(n,x,i,sc,scaling);
    evalfun = evalfun + 1;
end
cs = cs';

al  = fs;
nal = gs;
for i = 1:m
    if ( equatn(i) == true || lambda(i) + rho * cs(i) > 0 )
        al = al + cs(i) * ( lambda(i) + 0.5 * rho * cs(i) );

        [~,ncs,flag] = sevalnc(n,x,i,sc,scaling);
        evalgrad = evalgrad + 1;

        nal = nal + ( lambda(i) + rho * cs(i) ) * ncs;
    else
        al = al - 0.5 * lambda(i)^2 / rho;
    end
end