function [nl,ni] = evalnlni(n,m,x,g,c,lambda,I,sc,scaling)

nl = [];
nl = g;

ni = zeros(n,1);

c(I) = max(c(I),0);

for i = 1:m
    [~,ncs,flag] = sevalnc(n,x,i,sc,scaling);

    nl = nl + lambda(i) * ncs;

    ni = ni + c(i) * ncs;
end