function [nal,store] = evalnal(x,store,n,m,equatn,rho,lambar,sc,scaling)

global evalfun evalgrad nevalnal 

nevalnal = nevalnal + 1;

[~,gs,flag] = sevalg(n,x,sc,scaling);
evalgrad = evalgrad + 1;

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

nal = [];

nal = gs;
for i = 1:m
    if ( equatn(i) == true || lambar(i) + rho * cs(i) > 0 )
        [~,ncs,flag] = sevalnc(n,x,i,sc,scaling);
        evalgrad = evalgrad + 1;

        nal = nal + ( lambar(i) + rho * cs(i) ) * ncs;
    end
end

[nal] = reshapevector(nal);