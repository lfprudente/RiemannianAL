function [al,store] = evalal(x,store,n,m,equatn,rho,lambar,sc,scaling)

global evalfun nevalal 

nevalal = nevalal + 1;

[~,fs,flag] = sevalf(n,x,sc,scaling);
evalfun = evalfun + 1;

if ~isfield(store, 'c')
    cs = [];
    for i = 1:m
        [~,cs(i),flag] = sevalc(n,x,i,sc,scaling);
        evalfun = evalfun + 1;
    end
    cs = cs';

    store.c = cs;       
end

cs = store.c;

al = fs;
for i = 1:m
    if ( equatn(i) == true || lambar(i) + rho * cs(i) > 0 )
        al = al + cs(i) * ( lambar(i) + 0.5 * rho * cs(i) );
    else
        al = al - 0.5 * ( lambar(i) )^2 / rho;
    end
end