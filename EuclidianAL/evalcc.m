function [c,flag] = evalcc(n,x,ind)

global problemID pairs nballs a b N rankY

% ==================================================================

if ( problemID == 1 )

    flag = 0;
    
    if  ( ind == 1 )
        c = dot(x,x) - 1;
        return
    end

    if  ( ind == 2 )
        c =  sum(x) + x(2);
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

    nck = nchoosek(nballs,2);
    cte = ( b / a )^2;

    if ( ind >= 1 && ind <= nballs )

        i = ind; 

        c = r^2 - ( s(i) - 1 )^2 * b^2 * ( cte * u(i)^2 + v(i)^2 );

    elseif ( ind >= nballs + 1 && ind <= nballs + nck )

        i = pairs(ind-nballs, 1);
        j = pairs(ind-nballs, 2);

        c = 4 * r^2 - a^2 * ( ( 1 + ( s(i) - 1 ) * cte ) * u(i) - ( 1 + ( s(j) - 1 ) * cte ) * u(j) )^2 ...
            - b^2 * ( s(i) * v(i) - s(j) * v(j) )^2;

    else

        i = ind - ( nballs + nck );
        c = u(i)^2 + v(i)^2 - 1;
    end

    return
end

% ==================================================================

if ( problemID == 3 )

    flag = 0;

    % Reshape x into the matrix Y with dimensions [N, rankY]

    Y = reshape(x, [N, rankY]);

    % Calculate the number of constraints from Y^T Y = I

    Stconst = rankY * ( rankY + 1 )/2;

    if ( ind <= Stconst )

        % Initialize the counters for the constraints
        current_constraint = 0;
    
        % Loop through columns for constraints Y^T Y = I
        for j = 1:rankY
            % Loop through rows (only considering upper triangular part including diagonal)
            for i = 1:j
                current_constraint = current_constraint + 1;
                if ( ind == current_constraint )
                    if ( i == j )
                        % Normalization constraint
                        c = norm(Y(:, i))^2 - 1;
                    else
                        % Orthogonality constraint
                        c = dot(Y(:, i), Y(:, j));
                    end
                    return
                end
            end
        end

    else

        % Loop through rows for constraints YY^T e = e

        c = sum(Y(ind - Stconst, :) * Y') - 1;
    end

end