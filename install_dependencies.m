function install_dependencies()

    % Installing dependencies

    % Define paths
    asaWrapperDir = fullfile(pwd, 'ASA_CG_matlabWrapper');
    asaCGDir = fullfile(pwd, 'ASA_CG-3.0');
    manoptDir = fullfile(pwd, 'manopt');
    
    % Compile asa_wrapper.c
    fprintf('Compiling ASA interface with Matlab...\n');
    mexCmd = sprintf('mex -DMEXPRINTF -DVER30 %s/asa_wrapper.c %s/asa_cg.c -I%s -largeArrayDims', ...
        asaWrapperDir, asaCGDir, asaCGDir);
    eval(mexCmd);
    
    % Add ASA_CG_matlabWrapper to MATLAB path
    fprintf('Adding ASA_CG_matlabWrapper to MATLAB path...\n\n');
    addpath(asaWrapperDir);

     % Test ASA installation
    test_asa = input('\n Do you want to test ASA installation with a sample problem? [Y/N] ', 's');
    if strcmpi(test_asa, 'Y') || strcmpi(test_asa, 'y')
        fprintf('Running sample problem for ASA...\n');
        cd('ASA_CG_matlabWrapper');
        run_sample_problem_ASA();
        cd('..');
    end
    
    % Run importmanopt.m to set up Manopt
    fprintf('\n\nInstalling Manopt...\n');
    currentDir = pwd;
    cd(manoptDir);
    importmanopt();
    cd(currentDir);

    % Test Manopt installation
    test_manopt = input('\n Do you want to test Manopt installation with a sample problem? [Y/N] ', 's');
    if strcmpi(test_manopt, 'Y') || strcmpi(test_manopt, 'y')
        fprintf('Running sample problem for Manopt...\n');
        run_sample_problem_manopt();
    end
    
    % Save the MATLAB path
    fprintf('Saving the MATLAB path...\n');
    if savepath() == 0
        fprintf('Path saved successfully.\n');
    else
        fprintf('Failed to save the path. You might need administrator privileges.\n');
    end
    
    % Confirm installation
    fprintf('\n Dependencies installed successfully.\n');
end

function run_sample_problem_ASA()
    fprintf('Running a sample problem for ASA...\n');
    n = 100;
    lo = zeros(n,1);    % lower bound
    hi = ones(n,1);     % upper bound
    x  = ones(n,1);     % initial guess; necessary so that it knows the size of the problem
    tol = 10^(-8);      % tolerance to optimality

    fcn  = 'asa_fcn';
    grad = 'asa_grad'; 
    fcnGrad = 'asa_fcnGrad'; % optional

    opts = [];
    opts.PrintParms = false;
    opts.PrintLevel = 0;
    opts.PrintFinal = false;
    CGopts = struct('PrintParms',false);

    fgopt.n = n;

    [x,status,statistics] = asa_wrapper( x, lo, hi, tol,fcn, grad, fcnGrad, opts, CGopts, fgopt);
    
    fprintf('Sample problem for ASA completed.\n');
    fprintf('Status: %d\n', status);
    fprintf('Function value: %f\n', statistics.f);
    fprintf('Gradient norm: %e\n', statistics.gnorm);
end

function run_sample_problem_Manopt()
    fprintf('Running a sample problem for Manopt...\n');
    % Generate random problem data.
    n = 1000;
    A = randn(n);
    A = .5*(A+A');

    % Create the problem structure.
    manifold = spherefactory(n);
    problem.M = manifold;

    % Define the problem cost function and its Euclidean gradient.
    problem.cost  = @(x) -x'*(A*x);
    problem.egrad = @(x) -2*A*x;      % notice the 'e' in 'egrad' for Euclidean

    % Numerically check gradient consistency (just once, optional).
    checkgradient(problem); pause;

    % Solve.
    [x, xcost, info, options] = rlbfgs(problem);

    % Display some statistics.
    figure;
    semilogy([info.iter], [info.gradnorm], '.-');
    xlabel('Iteration number');
    ylabel('Norm of the Riemannian gradient of f');
    
    fprintf('Sample problem for Manopt completed.\n');
end
