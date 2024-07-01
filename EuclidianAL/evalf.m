function [f,flag] = evalf(n,x)

global problemID N rankY D

% ==================================================================

if ( problemID == 1 )
    
    flag = 0;
    
    f = - sum( x.^8 + x );

    return
end

% ==================================================================

if ( problemID == 2 )
    
    flag = 0;

    f = - x(end);

    return
end

% ==================================================================

if ( problemID == 3 )
    
    flag = 0;

    Y = reshape(x, [N, rankY]);

    f = 0;
    for j = 1:rankY
        f = f + dot(Y(:,j),D*Y(:,j));
    end
    f = - f / 2;

    return
end