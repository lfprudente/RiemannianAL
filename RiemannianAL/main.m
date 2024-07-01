 clear all
% clc

global problemID pairs a b nballs N rankY D

format long

problemID = 2;

% ==================================================================

if ( problemID == 1 )

    rng(2024);

    % Number of variables
    
    n = 50;
    
    % Constraints

    m = 1;
    
    equatn(1:m) = false;
    lambda = zeros(m,1);
    
    E = (equatn==true);
    I = (equatn==false);

    manifold = spherefactory(n);
    problem.M = manifold;
    
    % Checking derivatives?
    
    checkder = false;
    
    % Scale the problem?
    
    scaling = false;

    % Set initial guess

    x = randn(n,1);
    x = x / norm(x);

    if ( checkder )
        l(1:n) = - Inf;
        u(1:n) =   Inf;
        checkd(n,m,x,l,u)
    end

    % Call Riemannian-AL

    [x,lambda,f,iter,time,nevalal,nevalnal,alinfo] = auglag(n,m,x,equatn,lambda,scaling,problem);

end

% ==================================================================

if ( problemID == 2 )

    rng(2024);

    % Set the number of balls to be packed

    nballs = 10;

    fprintf('nballs = %i \n',nballs)

    % Print solution?
    
    print = true;

    % Problem data
        
    a = 2;
    b = 1;

    % Number of variables

    n = 3 * nballs + 1;

    % Constraints

    m = nchoosek(nballs,2) + 3 * nballs + 1;
    
    equatn(1:m) = false;
    lambda = zeros(m,1);
    
    E = (equatn==true);
    I = (equatn==false);

    sphere = spherefactory(2);
    manifold = productmanifold(struct('uv',powermanifold(sphere, nballs),'s', euclideanfactory(nballs),'r', euclideanfactory(1)));
    problem.M = manifold;
    
    pairs = nchoosek(1:nballs, 2);

    % Checking derivatives?
    
    checkder = false;
    
    % Scale the problem?
    
    scaling = true;

    % Set initial guess

    for i = 1:nballs
        x0 = 2 * rand(2,1) - 1;
        x0 = x0 / norm(x0);
        x.uv{i}(:,1) = x0;
    end
    x.uv = x.uv';
    x.s = rand(nballs,1);
    x.r = rand;
    
    if ( checkder )
        xa = cell2mat(x.uv);
        xb = [];
        xb(1:nballs,1) = xa(1:2:end,1);
        xb(nballs+1:2*nballs,1) = xa(2:2:end,1);
        xz = [xb; x.s; x.r];
        checkd(n,m,xz,-Inf(n,1),Inf(n,1))
    end

    % Call Riemannian-AL

    [x,lambda,f,iter,time,nevalal,nevalnal,alinfo] = auglag(n,m,x,equatn,lambda,scaling,problem);

    if ( print )
        
        xa = [];
        xa = cell2mat(x.uv);
    
        u(1:nballs,1) = xa(1:2:end,1);
        v(1:nballs,1) = xa(2:2:end,1);
        s = x.s;
        r = x.r;
        
        xc = [];
        yc = [];
        for i = 1:nballs
            xc(i) = ( 1 + ( s(i) - 1 ) * b^2/a^2 ) * a * u(i);
            yc(i) = s(i) * b * v(i);
        end
    
        centers = [xc' yc'];
    
        plot_ellipse_and_circles(a, b, nballs, centers, r, [0.94,0.831,0.827])
    end
    
end

% ==================================================================

if ( problemID == 3 )

    % Number of points

    N = 200;

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
    
    m = N + rankY * N;
    
    equatn(1:N)   = true;
    equatn(N+1:m) = false;
    lambda = zeros(m,1);
    
    E = (equatn==true);
    I = (equatn==false);
    
    M = stiefelfactory(N, rankY);
    problem.M = M;

    % Set initial guess

    x = M.rand();
    
    if ( checkder )
        checkd(n,m,x(:),- Inf(n,1),Inf(n,1))
    end

    % Call Riemannian-AL
            
    [x,lambda,f,iter,time,nevalal,nevalnal,alinfo] = auglag(n,m,x,equatn,lambda,scaling,problem);
    
    if ( print )
    
        if ( size(x,2) == 1 )
            Y = reshape(x, [N, rankY]);
        else
            Y = x;
        end
    
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