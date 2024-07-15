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
    
    checkder = false;
    
    % Scale the problem?
    
    scaling = true;

    % Set initial guess

    x = randn(n,1);
    x = x / norm(x);

    if ( checkder )
        checkd(n,m,x,l,u)
    end

    % Call Euclidian-AL

    [x,lambda,f,iter,csupn,time,nevalal,nevalnal,alinfo] = auglag(n,m,l,u,x,equatn,lambda,scaling);

    return
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
    
    [x,lambda,f,iter,csupn,time,nevalal,nevalnal,alinfo] = auglag(n,m,l,u,x,equatn,lambda,scaling);
            
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
    
    m = rankY * ( rankY + 1 )/2 + N;
    
    equatn(1:m) = true;
    lambda = zeros(m,1);
    
    E = (equatn==true);
    I = (equatn==false);

    % Set lower bounds and upper bounds

    l = zeros(n,1);
    u = Inf(n,1);

    % Set initial guess

    pos = randi(k,N,1);

    x = zeros(N,k);

    for i = 1:k
        x(pos==i,i) = 1;
    end

    for i = 1:k
        x(:,i) = x(:,i)/norm(x(:,i));
    end   

    x = x(:);

    if ( checkder )
        checkd(n,m,x,- Inf(n,1),Inf(n,1))
    end

    x = max( l, min( x, u ) );

    % Call Euclidian-AL
            
    [x,lambda,f,iter,csupn,time,nevalal,nevalnal,alinfo] = auglag(n,m,l,u,x,equatn,lambda,scaling);

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