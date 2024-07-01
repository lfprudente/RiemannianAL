function [nc,flag] = evalnc(n,x,ind)

global problemID centers dim nSph pairs nballs a b N rankY

% ==================================================================

if ( problemID == 1 )

    flag = 0;
    
    if  ( ind == 1 )

        nc = 2 * x;
        return
    end

    if  ( ind == 2 )

        nc(1:n) =  1;
        nc(2) =  2;

        nc = nc';
        return
    end
   
end

% ==================================================================

if ( problemID == 2 )

    flag = 0;

    u(1:nballs,1) = x(1:nballs,1);
    v(1:nballs,1) = x(nballs+1:2*nballs,1);
    s = x(2*nballs+1:3*nballs,1);
    r = x(3*nballs+1,1);

    nc = zeros(n, 1);

    nck = nchoosek(nballs,2);
    cte = ( b / a )^2;

    if ( ind >= 1 && ind <= nballs )

        i = ind;

        Pui = i;
        Pvi = i + nballs;
        Psi = i + 2 * nballs;
        Pr  = 3 * nballs + 1;

        nc(Pui) = - ( s(i) - 1 )^2 * b^2 * cte * 2 * u(i);
        nc(Pvi) = - ( s(i) - 1 )^2 * b^2 * 2 * v(i);
        nc(Psi) = - 2 * ( s(i) - 1 )* b^2 * ( cte * u(i)^2 + v(i)^2 );
        nc(Pr)  = 2 * x(end);

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

        xxi = ( 1 + ( s(i) - 1 ) * cte ) * u(i);
        xxj = ( 1 + ( s(j) - 1 ) * cte ) * u(j);
        yyi = s(i) * v(i);
        yyj = s(j) * v(j);

        nc(Pui) = - 2 * a^2 * ( xxi - xxj ) * ( 1 + ( s(i) - 1 ) * cte );
        nc(Pvi) = - 2 * b^2 * ( yyi - yyj ) * s(i);
        nc(Psi) = - 2 * a^2 * ( xxi - xxj ) * cte * u(i) - 2 * b^2 * ( yyi - yyj ) * v(i);
        nc(Puj) =   2 * a^2 * ( xxi - xxj ) * ( 1 + ( s(j) - 1 ) * cte );
        nc(Pvj) =   2 * b^2 * ( yyi - yyj ) * s(j);
        nc(Psj) =   2 * a^2 * ( xxi - xxj ) * cte * u(j) + 2 * b^2 * ( yyi - yyj ) * v(j);
        nc(Pr)  =   8 * r;

    
    else
        i = ind - ( nballs + nck );

        Pui = i;
        Pvi = i + nballs;

        nc(Pui) = 2 * u(i);
        nc(Pvi) = 2 * v(i);
    end

    return
end

% ==================================================================


if ( problemID == 3 )

    flag = 0;
    
    % Reshape x into the matrix Y with dimensions [N, rankY]
    Y = reshape(x, [N, rankY]);
    
    % Calculate the number of constraints from Y^T Y = I
    Stconst = rankY * (rankY + 1) / 2;
    
    % Initialize the gradient vector with zeros
    nc = zeros(n, 1);
    
    if ind <= Stconst
        % Initialize the counter for the constraints
        current_constraint = 0;
        
        % Loop through columns for constraints Y^T Y = I
        for j = 1:rankY
            % Loop through rows (only considering upper triangular part including diagonal)
            for i = 1:j
                current_constraint = current_constraint + 1;
                if ind == current_constraint
                    if i == j
                        % Normalization constraint
                        for k = 1:N
                            nc((i-1)*N + k) = 2 * Y(k, i);
                        end
                    else
                        % Orthogonality constraint
                        for k = 1:N
                            nc((i-1)*N + k) = Y(k, j);
                            nc((j-1)*N + k) = Y(k, i);
                        end
                    end
                    return
                end
            end
        end
    else
        
        % Loop through rows for constraints YY^T e = e

        i = ind - Stconst;

        val = repmat(Y(i, :), N, 1);
        val(i, :) = val(i, :) + ones(N,1)'*Y;
        nc = val(:);
    end
   
end