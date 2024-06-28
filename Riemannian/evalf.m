function [f,flag] = evalf(n,x)

global problemID rankY D

% ==================================================================

if ( problemID == 1 )
    
    flag = 0;
    
    f = - sum( x.^8 + x );

    return
end

% ==================================================================

if ( problemID == 2 )
    
    flag = 0;

    if ( ~isstruct(x) )
        [x] = reshapevector(x);
    end

    f = - x.r;

    return
end

% ==================================================================

if ( problemID == 3 )
    
    flag = 0;
    
    if ( size(x,2) == 1 )
        [x] = reshapevector(x);
    end

    f = 0;
    for j = 1:rankY
        f = f + dot(x(:,j),D*x(:,j));
    end
    f =  - f / 2;

    return
end