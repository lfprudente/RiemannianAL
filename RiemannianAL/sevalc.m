function [c,cs,flag] = sevalc(n,x,ind,sc,scaling)

[c,flag] = evalcc(n,x,ind);

if ( scaling == true )
    cs = c * sc.c(ind);
else
    cs = c;
end