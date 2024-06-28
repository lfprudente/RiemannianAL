function [nc,flag] = evalnc(n,x,ind)

global problemID pairs nballs a b N

% ==================================================================

if ( problemID == 1 )

    flag = 0;

    if  ( ind == 1 )

        nc(1:n) =  1;
        nc(2) =  2;

        nc = nc';
        return
    end
end

% ==================================================================

if ( problemID == 2 )

    flag = 0;

    if ( ~isstruct(x) )
        [x] = reshapevector(x);
    end

    nc = zeros(n, 1);

    nck = nchoosek(nballs,2);
    cte = ( b / a )^2;

    if ( ind >= 1 && ind <= nballs )

        i = ind;

        Pui = i;
        Pvi = i + nballs;
        Psi = i + 2 * nballs;
        Pr  = 3 * nballs + 1;

        nc(Pui) = - ( x.s(i) - 1 )^2 * b^2 * cte * 2 * x.uv{i}(1);
        nc(Pvi) = - ( x.s(i) - 1 )^2 * b^2 * 2 * x.uv{i}(2);
        nc(Psi) = - 2 * ( x.s(i) - 1 )* b^2 * ( cte * x.uv{i}(1)^2 + x.uv{i}(2)^2 );
        nc(Pr)  = 2 * x.r;

    elseif ( ind >= nballs + 1 && ind <= nballs + nck )

        i = pairs(ind-nballs, 1);
        j = pairs(ind-nballs, 2);

        Pui = i;
        Pvi = i + nballs;
        Psi = i + 2 * nballs;
        Puj = j;
        Pvj = j + nballs;
        Psj = j + 2 * nballs;
        Pr  = 3 * nballs + 1;

        xxi = ( 1 + ( x.s(i) - 1 ) * cte ) * x.uv{i}(1);
        xxj = ( 1 + ( x.s(j) - 1 ) * cte ) * x.uv{j}(1);
        yyi = x.s(i) * x.uv{i}(2);
        yyj = x.s(j) * x.uv{j}(2);

        nc(Pui) = - 2 * a^2 * ( xxi - xxj ) * ( 1 + ( x.s(i) - 1 ) * cte );
        nc(Pvi) = - 2 * b^2 * ( yyi - yyj ) * x.s(i);
        nc(Psi) = - 2 * a^2 * ( xxi - xxj ) * cte * x.uv{i}(1) - 2 * b^2 * ( yyi - yyj ) * x.uv{i}(2);
        nc(Puj) =   2 * a^2 * ( xxi - xxj ) * ( 1 + ( x.s(j) - 1 ) * cte );
        nc(Pvj) =   2 * b^2 * ( yyi - yyj ) * x.s(j);
        nc(Psj) =   2 * a^2 * ( xxi - xxj ) * cte * x.uv{j}(1) + 2 * b^2 * ( yyi - yyj ) * x.uv{j}(2);
        nc(Pr)  =   8 * x.r;

    elseif ( ind >= nballs + nck + 1 && ind <= 2 * nballs + nck )

        i   = ind - ( nballs + nck );
        Psi = i + 2 * nballs;

        nc(Psi) = - 1;

    elseif ( ind >= 2 * nballs + nck + 1 && ind <= 3 * nballs + nck )

        i   = ind - ( 2 * nballs + nck );
        Psi = i + 2 * nballs;

        nc(Psi) = 1;
    else
        Pr  = 3 * nballs + 1;

        nc(Pr) = - 1;
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

    % Initialize the gradient vector with zeros
    nc = zeros(n, 1);

    if ( ind <= N )
        % Gradient of the constraints YY^T e = e

        val = repmat(x(ind, :), N, 1);
        val(ind, :) = val(ind, :) + ones(N,1)'*x;
        nc = val(:);

    else
        i = ind - N;
        nc(i) = - 1;
    end
end