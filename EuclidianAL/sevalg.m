function [g,gs,flag] = sevalg(n,x,sc,scaling)

[g,flag] = evalg(n,x);

if ( scaling == true )
    gs = g * sc.f;
else
    gs = g;
end