clear all
% clc

global problemID pairs a b nballs N rankY D

format long

problemID = 3;

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

    [x,lambda,f,iter,csupn,time,nevalal,nevalnal,alinfo] = auglag(n,m,x,equatn,lambda,scaling,problem);

end

% ==================================================================

if ( problemID == 2 )

    rng(2024);

    % Set the number of balls to be packed

    nballs = 3;

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

    [x,lambda,f,iter,csupn,time,nevalal,nevalnal,alinfo] = auglag(n,m,x,equatn,lambda,scaling,problem);

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

    problem_name = 'synthetic';

    if ( strcmp(problem_name, 'wdbc') == 1 )
        N = 569;
        k = 2;
    elseif ( strcmp(problem_name, 'cloud') == 1 )
        N = 2048;
        k = 2;
    elseif ( strcmp(problem_name, 'ecoli') == 1 )
        N = 336;
        k = 8;
    elseif ( strcmp(problem_name, 'ionosphere') == 1 )
        N = 351;
        k = 2;
    elseif ( strcmp(problem_name, 'iris') == 1 )
        N = 150;
        k = 3;
    elseif ( strcmp(problem_name, 'parkinsons') == 1 )
        N = 195;
        k = 2;
    elseif ( strcmp(problem_name, 'pima') == 1 )
        N = 768;
        k = 2;
    elseif ( strcmp(problem_name, 'raisin') == 1 )
        N = 900;
        k = 2;
    elseif ( strcmp(problem_name, 'seeds') == 1 )
        N = 210;
        k = 3;
   elseif ( strcmp(problem_name, 'SPECTF') == 1 )
        N = 267;
        k = 2;  
    elseif ( strcmp(problem_name, 'transfusion') == 1 )
        N = 748;
        k = 2;
    elseif ( strcmp(problem_name, 'thyroid') == 1 )
        N = 215;
        k = 3; 
    elseif ( strcmp(problem_name, 'wine') == 1 )
        N = 178;
        k = 3;
    elseif ( strcmp(problem_name, 'synthetic') == 1 )
        N = 500;
        k = 3;

        % Print solution?
        print = true;
    end

    addpath(fullfile('..','kmeansdata'));

    [data,true_labels] = datacleaner(problem_name,N,k);

    [nrow,ncol] = size(data);

    fprintf('Problem    : %s\n', problem_name);
    fprintf('Data(N)    : %d\n', nrow);
    fprintf('Features(l): %d\n', ncol);
    fprintf('Clusters(k): %d\n\n', k);

    rng(2024);

    % Set matrix D
   
    D = data * data.';

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

    pos = randi(k,N,1);

    x = zeros(N,k);

    for i = 1:k
        x(pos==i,i) = 1;
    end

    for i = 1:k
        x(:,i) = x(:,i)/norm(x(:,i));
    end        
    
    if ( checkder )
        checkd(n,m,x(:),- Inf(n,1),Inf(n,1))
    end

    % Call Riemannian-AL
            
    [x,lambda,f,iter,csupn,time,nevalal,nevalnal,alinfo] = auglag(n,m,x,equatn,lambda,scaling,problem);

    if ( alinfo == 0 && strcmp(problem_name, 'synthetic') ~= 1 )
        [accuracy] = analyses(x,N,k,true_labels);
        fprintf('Feas    : %7.1e \n', csupn);
        fprintf('Accuracy: %.1f \n', accuracy);
    end

    if ( strcmp(problem_name, 'synthetic') == 1 && print )
        plot_kmeans(x,N,k,data)
    end

    rmpath(fullfile('..','kmeansdata'));
end