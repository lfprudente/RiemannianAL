function [x,lambda,f,outiter,time,nevalal,nevalnal,alinfo] = auglag(n,m,l,u,x,equatn,lambda,scaling)

global evalfun evalgrad nevalal nevalnal

%  Parameters

epsopt    = 10^(-5);
epsfeas   = 10^(-5);
efstain   = sqrt(epsfeas);
eostain   = epsopt^(3/2);
epsfeas12 = sqrt( epsfeas );
epsopt12  = sqrt( epsopt );

lammin  = - 10^20; 
lammax  =   10^20;

rhomax  = 10^20;
rhofrac = 0.5;
rhomult = 10;

maxoutit = 50;
maxinnit = 50000;

E = ( equatn == true  );
I = ( equatn == false );

if ( scaling == true )
    [g,flag] = evalg(n,x);
    sc.f = 1 / max( 1, norm(g,"inf") );
    
    for i = 1:m
        [nc,flag] = evalnc(n,x,i);
        sc.c(i) = 1 / max( 1, norm(nc,"inf") );
    end
else
    sc.f = 1;
    sc.c(1:m) = 1;
end

% Option structures from ASA interface

fcn  = 'evalal';
grad = 'evalnal'; 
fcnGrad = 'evalalnal';

opts = [];
opts.PrintParms = false;
opts.PrintLevel = 0;
opts.PrintFinal = false;

CGopts.PrintParms = false;

fgopt.n = n;
fgopt.m = m;
fgopt.e = equatn;
fgopt.scaling = scaling;
fgopt.sc = sc;

% ==================================================================
% Print initial information
% ==================================================================

fprintf('----------------------------------------------------------------------\n')
fprintf('                  Augmented Lagrangian method                         \n')
fprintf('----------------------------------------------------------------------\n')
fprintf('Number of variables                : %i \n',n)
fprintf('Number of constraints              : %i \n',m)
fprintf('Optimality tolerance               : %.0e \n',epsopt)
fprintf('Feasibility tolerance              : %.0e \n',epsfeas)
if ( scaling == true )
    fprintf('Objective function scale factor    : %.0e \n',sc.f)
    fprintf('Smallest constraints scale factor  : %.0e \n',min(sc.c))
end

% Start timing

tic;

%  ==================================================================
%  Initialization
%  ==================================================================

%  Counters

outiter  = 0;
evalfun  = 0;
evalgrad = 0;

nevalal  = 0;
nevalnal = 0;

cifacnt = 0;

% Compute objective function value and gradient

[f,fs,flag] = sevalf(n,x,sc,scaling);
evalfun = evalfun + 1;

[~,gs,flag] = sevalg(n,x,sc,scaling);
evalgrad = evalgrad + 1;

% Compute constraints 

c  = [];
cs = [];
for i = 1:m
    [c(i),cs(i),flag] = sevalc(n,x,i,sc,scaling);
    evalfun = evalfun + 1;
end
c  = c';
cs = cs';

% Compute safeguarded Lagrange multipliers

lambar = max( lammin, min( lambda, lammax ) );

% Compute complementarity and feasibility violations

snormprev = Inf;

csupn  = max( norm(abs(c(E)),Inf), norm(max(c(I),0),Inf) );
cssupn = max( norm(abs(cs(E)),Inf), norm(max(cs(I),0),Inf) );
snorm  = max( norm(abs(cs(E)),Inf), norm(min(-cs(I),lambda(I)),Inf) );

% Compute continuous projected Lagrangian gradient norm and squared-infeasibility gradient norm

[nl,ni] = evalnlni(n,m,x,gs,cs,lambda,I,sc,scaling);
evalgrad = evalgrad + m;

nlpsupn = norm( max( l, min( x  - nl, u ) ) - x, Inf);

% Compute squared-infeasibility gradient norm

nisupn = norm( max( l, min( x  - ni, u ) ) - x, Inf);

% Save best solution

fb       = f;
csupnb   = csupn;
snormb   = snorm;
nlpsupnb = nlpsupn;
xb       = x;
lambdab  = lambda;

% ==================================================================
% Main loop
% ==================================================================

while (1)

    % ==================================================================
    % Print information of this iteration
    % ==================================================================
    
    % Print information
    
    if ( scaling == true )
        if ( mod(outiter,10) == 0 )
            fprintf('\n')
            fprintf('%3s %-6s  %-9s %-6s  %-6s  %+9s %-6s %-4s   %-6s %-5s \n','out','penalt','objective','infeas','scaled','scaled','infeas','norm','|Grad|','inner')
            fprintf('%3s %-6s  %-9s %-6s  %-6s  %-6s %-6s %-4s %-6s %-5s \n','ite','param ','function ','ibilty','obj-funct','infeas','+compl','graLag','infeas','stop ')
        end
        if ( outiter == 0 )
            fprintf('%3i   -    %10.3e %6.0e %10.3e %6.0e %6.0e %6.0e %6.0e  - \n',outiter,f,csupn,fs,cssupn,snorm,nlpsupn,nisupn)
        else
            fprintf('%3i %6.0e %10.3e %6.0e %10.3e %6.0e %6.0e %6.0e %6.0e  %i \n',outiter,rho,f,csupn,fs,cssupn,snorm,nlpsupn,nisupn,flagIS)
        end

    else
        if ( mod(outiter,10) == 0 )
            fprintf('\n')
            fprintf('%3s %-6s  %-9s %-6s %-6s %-4s   %-6s %-5s \n','out','penalt','objective','infeas','infeas','norm','|Grad|','inner')
            fprintf('%3s %-6s  %-9s %-6s %-6s %-4s %-6s %-5s \n','ite','param ','function ','ibilty','+compl','graLag','infeas','stop ')
        end
        if ( outiter == 0 )
            fprintf('%3i   -    %10.3e %6.0e %6.0e %6.0e %6.0e - \n',outiter,f,csupn,snorm,nlpsupn,nisupn)
        else
            fprintf('%3i %6.0e %10.3e %6.0e %6.0e %6.0e %6.0e %i \n',outiter,rho,f,csupn,snorm,nlpsupn,nisupn,flagIS)
        end
    end
    
    % ==================================================================
    % Test stopping criteria
    % ==================================================================
    
    % Test feasibility, optimality and complementarity
    
    if ( max( snorm, csupn ) <= epsfeas && nlpsupn <= epsopt ) 
        alinfo = 0;
        
        % Stop timing 
    
        time = toc;
    
        % Print information
    
        fprintf('\n')
        fprintf('Solution was found.\n')
        fprintf('Number of function evaluations   : %i\n',evalfun)
        fprintf('Number of gradient evaluations   : %i\n',evalgrad)
        fprintf('Number of AL evaluations         : %i\n',nevalal)
        fprintf('Number of AL gradient evaluations: %i\n',nevalnal)
        fprintf('CPU time(s)                      : %.1f \n',time)

        return
    end
    
    % Test whether we are at an infeasible point that is stationary for
    % the sum of the squared infeasibilities
    
%     if ( csupn > efstain && nisupn <= eostain ) 
%     
%         alinfo = 1;
%         
%         % Stop timing 
%         
%         time = toc;
%         
%         % Print information
%         
%         fprintf('\n')
%         fprintf('It seems that a stationary-of-the-infeasibility probably infeasible point was found.\n')
%         fprintf('Whether the final iterate is a solution or not requires further analysis.\n')
%         fprintf('Number of function evaluations: %i\n',evalfun)
%         fprintf('Number of gradient evaluations: %i\n',evalgrad)
%         fprintf('CPU time(s)                   : %.1f \n',time)
%         
%         return
%     end
    
    % Test whether the penalty parameter is too large
    
    if ( outiter > 0 && rho > rhomax )
    
        alinfo = 2;
        
        % Stop timing 
        
        time = toc;
        
        % Print information
        
        fprintf('\n')
        fprintf('The penalty parameter is too large. The problem may be \n')
        fprintf('infeasible or badly scaled. Further analysis is required.\n')
        fprintf('Number of function evaluations   : %i\n',evalfun)
        fprintf('Number of gradient evaluations   : %i\n',evalgrad)
        fprintf('Number of AL evaluations         : %i\n',nevalal)
        fprintf('Number of AL gradient evaluations: %i\n',nevalnal)
        fprintf('CPU time(s)                      : %.1f \n',time)
        
        return

    end
    
    % Test whether the number of iterations is exhausted
    
    if ( outiter >= maxoutit )
        alinfo = 3;
        
        % Stop timing 
        
        time = toc;
        
        % Print information
        
        fprintf('\n')
        fprintf('Maximum of iterations reached. The feasibility-complementarity and \n')
        fprintf('optimality tolerances could not be achieved.\n')
        fprintf('Number of function evaluations   : %i\n',evalfun)
        fprintf('Number of gradient evaluations   : %i\n',evalgrad)
        fprintf('Number of AL evaluations         : %i\n',nevalal)
        fprintf('Number of AL gradient evaluations: %i\n',nevalnal)
        fprintf('CPU time(s)                      : %.1f \n',time)
        
        return
    end
    
    % ==================================================================
    % Iteration
    % ==================================================================
    
    outiter = outiter + 1;
    
    % ==================================================================
    % Set penalty parameter
    % ==================================================================
    
    %if ( outiter == 1 || outiter == 2 ) 
    if ( outiter == 1 || outiter == 2 )
        [rho] = comprhoini(c,f,E,I);
    elseif ( max( snorm, csupn ) > epsfeas && snorm > rhofrac * snormprev )
        rho = rhomult * rho;
    elseif ( cifacnt >= 2 )
        rhoa = min( rhoa * rhomult, 1 );
        rhob = max( rhob / rhomult, 1 );
        rhoc = rhoc * rhomult;

        [rhoini] = comprhoini(cs,fs,E,I);

        rhoini = max( rhoa, min( rhoini, rhob ) );

        rho = rhoini;
    end
    
    % ==================================================================
    % Solve the augmented Lagrangian subproblem
    % ==================================================================
    
    % Set optimality requeriment for the subproblem
    
    if ( outiter == 1 ) 
        epsopk = sqrt( epsopt );
    elseif ( snorm <= epsfeas12 && nlpsupn <= epsopt12 )
        epsopk = min( rhofrac * nlpsupn, 0.1 * epsopk );
        epsopk = max( epsopk, epsopt );
    end
    
    if ( outiter == 1 )
        CGopts.maxit = 100;
    else
        CGopts.maxit = maxinnit;
    end
    
    % Call the inner-solver

    fgopt.r = rho;
    fgopt.l = lambar;

    [x,flagIS,~] = asa_wrapper(x,l,u,epsopk,fcn,grad,fcnGrad,opts,CGopts,fgopt);

    if ( outiter == 1 )
        flagIS = 0;
    end
    
    % ==================================================================
    % Prepare for the next iteration
    % ==================================================================
    
    % Compute objective function value and gradient
    
    [f,fs,flag] = sevalf(n,x,sc,scaling);
    evalfun = evalfun + 1;

    [~,gs,flag] = sevalg(n,x,sc,scaling);
    evalgrad = evalgrad + 1;
    
    % Compute constraints
    
    c  = [];
    cs = [];
    for i = 1:m
        [c(i),cs(i),flag] = sevalc(n,x,i,sc,scaling);
        evalfun = evalfun + 1;
    end
    c  = c';
    cs = cs';
    
    % Compute feasibility violation
    
    csupn  = max( norm(abs(c(E)),Inf), norm(max(c(I),0),Inf) );
    cssupn = max( norm(abs(cs(E)),Inf), norm(max(cs(I),0),Inf) );
    
    % Update Lagrange multipliers approximation

    lambda = lambar + rho * cs;
    lambda(I) = max( 0, lambda(I) );

    lambar = max( lammin, min( lambda, lammax ) );
    
    % Compute complementarity violation
    
    snormprev = snorm;
    
    snorm = max( norm(abs(cs(E)),Inf), norm(min(-cs(I),lambda(I)),Inf) );

    % Compute continuous projected Lagrangian gradient and squared-infeasibility gradient

    [nl,ni] = evalnlni(n,m,x,gs,cs,lambda,I,sc,scaling);
    evalgrad = evalgrad + m;

     % Compute continuous projected Lagrangian gradient norm
    
    nlpsupn = norm( max( l, min( x  - nl, u ) ) - x, Inf);
    
    % Compute squared-infeasibility gradient norm
    
    nisupn = norm( max( l, min( x  - ni, u ) ) - x, Inf);

    % Track consecutve failures of the inner solver at feasible points

    if ( max( snorm, csupn ) <= epsfeas && flagIS ~= 0 )
        cifacnt = cifacnt + 1;
    else
        cifacnt = 0;
    end
    
    % Save best solution
    
    if ( ( csupnb > epsfeas && csupn < csupnb ) || ...
     ( csupnb <= epsfeas && csupn <= epsfeas && f < fb ) ) 
    
        fb       = f;
        csupnb   = csupn;
        snormb   = snorm;
        nlpsupnb = nlpsupn;
        xb       = x;
        lambdab  = lambda;
    
    end
    
    % ==================================================================
    % Iterate
    % ==================================================================

end

% ==================================================================
% End of main loop
% ==================================================================