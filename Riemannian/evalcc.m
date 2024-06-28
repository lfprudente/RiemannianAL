function [c,flag] = evalcc(n,x,ind)

global problemID pairs nballs a b N

% ==================================================================

if ( problemID == 1 )

    flag = 0;

    if ( ind == 1 ) 
        c =  sum(x) + x(2);
    end

end

% ==================================================================

if ( problemID == 2 )

    flag = 0;

    if ( ~isstruct(x) )
        [x] = reshapevector(x);
    end

    nck = nchoosek(nballs,2);
    cte = ( b / a )^2;

    if ( ind >= 1 && ind <= nballs )

        i = ind; 

        c = x.r^2 - ( x.s(i) - 1 )^2 * b^2 * ( cte * x.uv{i}(1)^2 + x.uv{i}(2)^2 );

    elseif ( ind >= nballs + 1 && ind <= nballs + nck )

        i = pairs(ind-nballs, 1);
        j = pairs(ind-nballs, 2);

        c = 4 * x.r^2 - a^2 * ( ( 1 + ( x.s(i) - 1 ) * cte ) * x.uv{i}(1) - ( 1 + ( x.s(j) - 1 ) * cte ) * x.uv{j}(1) )^2 ...
            - b^2 * ( x.s(i) * x.uv{i}(2) - x.s(j) * x.uv{j}(2) )^2;

    elseif ( ind >= nballs + nck + 1 && ind <= 2 * nballs + nck )

        i = ind - ( nballs + nck );

        c = - x.s(i);

    elseif ( ind >= 2 * nballs + nck + 1 && ind <= 3 * nballs + nck )

        i = ind - ( 2 * nballs + nck );

        c = x.s(i) - 1;

    else

        c = - x.r;

    end

    return
end

% ==================================================================

if ( problemID == 3 )

    flag = 0;

    % Reshape x into the matrix Y with dimensions [N, rankY]

    if ( size(x,2) == 1 )
        [x] = reshapevector(x);
    end

    if ( ind <= N )
        % Loop through rows for constraints YY^T e = e

        c = sum(x(ind, :) * x') - 1;
    else
        i = ind - N;
        x = x(:);
        c = - x(i);
    end

end