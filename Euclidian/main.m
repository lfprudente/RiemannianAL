clear all
% clc

global problemID pairs a b nballs N rankY D

problemID = 3;

% ==================================================================

if ( problemID == 1 )

    rng(2024);

    % Number of variables
    
    n = 50;
    
    % Constraints

    m = 2;
    
    equatn(1) = true;
    equatn(2) = false;

    lambda = zeros(m,1);
    
    E = (equatn==true);
    I = (equatn==false);
    
    % Set lower bounds and upper bounds

    l = - Inf(n,1);
    u =   Inf(n,1);

    % Checking derivatives?
    
    checkder = true;
    
    % Scale the problem?
    
    scaling = true;

    % Set initial guess

    x = randn(n,1);
    x = x / norm(x);

    if ( checkder )
        checkd(n,m,x,l,u)
    end

    % Call Euclidian-AL

    [x,lambda,f,iter,time,nevalal,nevalnal,alinfo] = auglag(n,m,l,u,x,equatn,lambda,scaling);

    return
end

% ==================================================================

if ( problemID == 2 )

    rng(2024);

    % Set the number of balls to be packed

    nballs = 20;

    fprintf('nballs = %i \n',nballs)

    % Print solution?
    
    print = true;

    % Problem data
        
    a = 2;
    b = 1;

    % Number of variables

    n = 3 * nballs + 1;

    % Constraints

    nck = nchoosek(nballs,2);

    m = nck + 2 * nballs;
            
    equatn(1:nck+nballs)   = false;
    equatn(nck+nballs+1:m) = true;

    lambda = zeros(m,1);
    
    E = (equatn==true);
    I = (equatn==false);

    % Set lower bounds and upper bounds

    % ui, vi variables

    l(1:2*nballs) = - 1;
    u(1:2*nballs) =   1;

    % si variables

    l(2*nballs+1:3*nballs) =  0;
    u(2*nballs+1:3*nballs) =  1;

    % r variable

    l(3*nballs+1) = 0;
    u(3*nballs+1) = Inf;

    l = l';
    u = u';

    pairs = nchoosek(1:nballs, 2);

    % Checking derivatives?
    
    checkder = false;
    
    % Scale the problem?
    
    scaling = true;

    % Set initial guess

    for i = 1:nballs
        xa = 2 * rand(2,1) - 1;
        xa = xa / norm(xa);
        x0.uv{i}(:,1) = xa;
    end
    x0.uv = x0.uv';
    x0.s = rand(nballs,1);
    x0.r = rand;

    xa = [];
    xa = cell2mat(x0.uv);
    xb = [];
    xb(1:nballs,1) = xa(1:2:end,1);
    xb(nballs+1:2*nballs,1) = xa(2:2:end,1);
    x = [xb; x0.s; x0.r];
        
    if ( checkder )
        checkd(n,m,x,l,u)
    end

    % Call Euclidian-AL
    
    [x,lambda,f,iter,time,nevalal,nevalnal,alinfo] = auglag(n,m,l,u,x,equatn,lambda,scaling);
        
    if ( print )

        u = [];
        v = [];
        s = [];
        r = [];
        
        u(1:nballs,1) = x(1:nballs,1);
        v(1:nballs,1) = x(nballs+1:2*nballs,1);
        s = x(2*nballs+1:3*nballs,1);
        r = x(3*nballs+1,1);
        
        xc = [];
        yc = [];
        for i = 1:nballs
            xc(i) = ( 1 + ( s(i) - 1 ) * b^2/a^2 ) * a * u(i);
            yc(i) = s(i) * b * v(i);
        end
    
        centers = [xc' yc'];
    
        plot_ellipse_and_circles(a, b, nballs, centers, r,[0.94,0.831,0.827])
    end
end

% ==================================================================

if ( problemID == 3 )

    % Number of points

    N = 300;

    if N == 200
        rng(2025);
    elseif N == 300
        rng(2026);
    elseif N == 400
        rng(2027);
    elseif N == 500
        rng(2028);
    elseif N == 800
        rng(2029);
    elseif N == 1000
        rng(2030);
    else
        rng(2024);
    end
        
    % Generate data points
    
    a = 1;
    b = 4;
    
    data = [];
    
    if ( N == 200 || N == 300 || N == 400 || N == 500 )

        % Number of clusters
        k = 3;

        % References points
        centroids = [3, 3; -3, -3; 6, -6];
    else
        % Number of clusters
        k = 2; 

        % References points
        centroids = [3, 3; -3, -3];
    end

    points_per_cluster = floor(N / k);
    
    for i = 1:k
        dispersion_factor = ( b - a ) * rand(points_per_cluster,1) + a;
        data = [data; centroids(i, :) + dispersion_factor .* randn(points_per_cluster, 2)];
    end
    
    % Add extra points if N is not exactly divisible by k

    remaining_points = N - points_per_cluster * k;
    if remaining_points > 0
        dispersion_factor = ( b - a ) * rand(remaining_points,1) + a;
        data = [data; centroids(1:remaining_points, :) + dispersion_factor .* randn(remaining_points, 2)];
    end
    
    rng(2024);

    % Set matrix D
   
    D = data * data.';

    % Print solution?

    print = true;
    
    % Number of variables

    rankY = k;
    
    n = N * rankY;
    
    % Checking derivatives?
    
    checkder = false;
    
    % Scale the problem?
    
    scaling = true;
    
    % Constraints
    
    m = rankY * ( rankY + 1 )/2 + N;
    
    equatn(1:m) = true;
    lambda = zeros(m,1);
    
    E = (equatn==true);
    I = (equatn==false);

    % Set lower bounds and upper bounds

    l = zeros(n,1);
    u = Inf(n,1);

    % Set initial guess

    M = stiefelfactory(N, rankY);
    problem.M = M;
    x = M.rand();

    x = x(:);

    if ( checkder )
        checkd(n,m,x,- Inf(n,1),Inf(n,1))
    end

    x = max( l, min( x, u ) );

    % Call Euclidian-AL
            
    [x,lambda,f,iter,time,nevalal,nevalnal,alinfo] = auglag(n,m,l,u,x,equatn,lambda,scaling);

    if ( print )
        Y = reshape(x, [N, rankY]);
        
        for i = 1:N
          [~,j] = max(Y(i,:));
          Y(i,:) = 0;
          Y(i,j) = 1;
        end

        centroids = (Y' * data) ./ sum(Y', 2);
  
        figure;
        hold on;
    
        colors = {'k','r','b'};
        marker = {'o','s','^'};
    
        axis off

        for j = 1:rankY
            cluster_points = data(Y(:, j) == 1, :);
            scatter(cluster_points(:, 1), cluster_points(:, 2), 36, colors{j},marker{j},'filled');
            %scatter(centroids(j, 1), centroids(j, 2), 300, 'y', 'pentagram', 'filled');
        end

    end
end