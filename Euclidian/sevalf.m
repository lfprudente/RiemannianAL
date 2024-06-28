function [f,fs,flag] = sevalf(n,x,sc,scaling)

[f,flag] = evalf(n,x);

if ( scaling == true )
    fs = f * sc.f;
else
    fs = f;
end