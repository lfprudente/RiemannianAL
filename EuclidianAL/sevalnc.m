function [nc,ncs,flag] = sevalnc(n,x,ind,sc,scaling)

[nc,flag] = evalnc(n,x,ind);

if ( scaling == true )
    ncs = nc * sc.c(ind);
else
    ncs = nc;
end