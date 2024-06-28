function realtest_AD1()
% Test AD for a real optimization problem on a product manifold (struct)

    % Verify that Manopt was indeed added to the Matlab path.
    if isempty(which('spherefactory'))
        error(['You should first add Manopt to the Matlab path.\n' ...
		       'Please run importmanopt.']);
    end
    
    % Verify that the deep learning tool box was installed
    assert(exist('dlarray', 'file') == 2, ['Deep learning tool box is '... 
    'needed for automatic differentiation.\n Please install the'...
    'latest version of the deep learning tool box and \nupgrade to Matlab'...
    ' R2021b if possible.'])
    
    % Generate the problem data.
    n = 100;
    A = randn(n);
    A = .5*(A+A');
    
    % Create the product manifold
    S = spherefactory(n);
    manifold.x = S;
    manifold.y = S;
    problem.M = productmanifold(manifold);
    
    % Define the problem cost function
    problem.cost  = @(X) -X.x'*(A*X.y);
    
    % Define the gradient and the hessian via automatic differentiation
    problem = manoptAD(problem);

    % Numerically check gradient and Hessian consistency.
    figure;
    checkgradient(problem);
    figure;
    checkhessian(problem);
    
    % Solve.
    [x, xcost, info] = trustregions(problem);          %#ok<ASGLU>
    
    % Test
    ground_truth = svd(A);
    distance = abs(ground_truth(1) - (-problem.cost(x)));
    fprintf('The distance between the ground truth and the solution is %e \n',distance);

    
end