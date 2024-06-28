function [g,flag] = evalg(n,x)

global problemID N rankY D

% ==================================================================

if ( problemID == 1 )

    flag = 0;

    g = - 8 * x.^7 - 1;

    return
end

% ==================================================================

if ( problemID == 2 )
    
    flag = 0;
    
    g = zeros(n,1);

    g(end) = - 1;

    return
end

% ==================================================================

if ( problemID == 3 )
    
    flag = 0;

    Y = reshape(x, [N, rankY]);

    g = D * Y;
    g = - g(:);
    
    return
end