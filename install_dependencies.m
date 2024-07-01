function install_dependencies()

    % Installing dependencies

    % Define paths
    asaWrapperDir = fullfile(pwd, 'ASA_CG_matlabWrapper');
    asaCGDir = fullfile(pwd, 'ASA_CG-3.0');
    manoptDir = fullfile(pwd, 'manopt');
    
    % Compile asa_wrapper.c
    fprintf('Compiling ASA interface with Matlab...\n');
    try
        mexCmd = sprintf('mex -outdir %s -DMEXPRINTF -DVER30 %s/asa_wrapper.c %s/asa_cg.c -I%s -largeArrayDims', ...
            asaWrapperDir, asaWrapperDir, asaCGDir, asaCGDir);
        eval(mexCmd);
        fprintf('\nASA interface compiled successfully.\n');

        % Test ASA installation
        test_asa = input('Do you want to test ASA installation with a sample problem? [Y/N] ', 's');
        if strcmpi(test_asa, 'Y') || strcmpi(test_asa, 'y')
            fprintf('Running sample problem for ASA...\n');
            cd('ASA_CG_matlabWrapper');
            [flagASA] = run_sample_problem_ASA();

            if ( flagASA == 0 )
                fprintf('\nASA converged successfully\n');
            else
                fprintf('\nASA did not converge as it should.\n');
            end

            cd('..');  
            fprintf('\nPress any key to continue...\n');
            pause;
        end

        % Add ASA_CG_matlabWrapper to MATLAB path
        fprintf('Adding ASA_CG_matlabWrapper to MATLAB path...\n\n');
        addpath(asaWrapperDir);
    catch
        fprintf('\nFailed to compile ASA interface.\n');
    end

     
    
    % Run importmanopt.m to set up Manopt
    fprintf('\n\n=================================\n');
    fprintf('Installing Manopt...\n');

    try
        currentDir = pwd;
        cd(manoptDir);
        importmanopt();
        cd(currentDir);

        % Test Manopt installation
        test_manopt = input('\nDo you want to test Manopt installation with a sample problem? [Y/N] ', 's');
        if strcmpi(test_manopt, 'Y') || strcmpi(test_manopt, 'y')
            fprintf('Running sample problem for Manopt...\n');
            [flagManopt] = run_sample_problem_manopt();

            if ( flagManopt == 0 )
                fprintf('\nManopt converged successfully\n');
            else
                fprintf('\nManopt did not converge as it should.\n');
            end
        end

    catch
        fprintf('\nFailed to install Manopt.\n');
    end
    
    % Save the MATLAB path
    response = input('Save paths for future Matlab sessions? [Y/N] ', 's');
    if strcmpi(response, 'Y')
        failed = savepath();
        if ~failed
            fprintf('Path saved successfully: no need to call install_dependencies next time.\n');
        else
            fprintf(['Failed to save the path. You might need administrator privileges \nto write on pathdef.m?\nPath not saved: ' ...
                     'please re-call install_dependencies next time.\n']);
        end
    else
        fprintf('Path not saved: please re-call install_dependencies next time.\n');
    end
    
    % Confirm installation
    if ( flagASA == 0 && flagManopt == 0 )
        fprintf('\n\nDependencies installed successfully.\n');
    elseif ( flagASA ~= 0 && flagManopt == 0 )
        fprintf('\nSomething went wrong with ASA. Manopt is OK.\n');
    elseif ( flagASA == 0 && flagManopt ~= 0 )
        fprintf('\nSomething went wrong with Manopt. ASA is OK.\n');
    else
        fprintf('\nSomething went wrong both with ASA and Manopt.\n');
    end
end

function [status] = run_sample_problem_ASA()
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
    opts.PrintLevel = 1;
    opts.PrintFinal = true;
    CGopts = struct('PrintParms',false);

    fgopt.n = n;
    
    [x,status,statistics] = asa_wrapper( x, lo, hi, tol,fcn, grad, fcnGrad, opts, CGopts, fgopt);
    
    fprintf('\n');
    fprintf('Sample problem for ASA completed.\n');
end

function [status] = run_sample_problem_manopt()
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
    % checkgradient(problem); pause;

    % Solve.
    [x, xcost, info, options] = rlbfgs(problem);
    
    status = - 1;
    if ( info(end).gradnorm <= options.tolgradnorm )
        status = 0;
    end

    % Display some statistics.
%     figure;
%     semilogy([info.iter], [info.gradnorm], '.-');
%     xlabel('Iteration number');
%     ylabel('Norm of the Riemannian gradient of f');
    
    fprintf('\n');
    fprintf('Sample problem for Manopt completed.\n');
end
