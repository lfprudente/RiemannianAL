function [g,flag] = evalg(n,x)

global problemID D


% ==================================================================

if ( problemID == 1 )

    flag = 0;

    g = - 8 * x.^7 - 1;

    return
end

% ==================================================================

if ( problemID == 2 )
    
    flag = 0;

    if ( ~isstruct(x) )
        [x] = reshapevector(x);
    end
    
    g = zeros(n, 1);

    g(end) = - 1;

    return
end

% ==================================================================

if ( problemID == 3 )
    
    flag = 0;
    
    if ( size(x,2) == 1 )
        [x] = reshapevector(x);
    end

    g = D * x;
    g = - g(:);

    return
end