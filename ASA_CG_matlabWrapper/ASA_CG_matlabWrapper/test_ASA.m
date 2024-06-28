%{
Stephen Becker, March 27 2012
Updated Jan 22 2018, downloads ASA_CG ver 3.0 instead of 2.2 now
  and fixes old dead links
  (not sure if it is 100% compatible with ver 3.0, but seems that
   at least the basic usage works)

Installation and sample usage of the asa_wrapper mex file "asa_wrapper.c"
We show how to use the helper files asa_fcn, asa_grad and asa_fcnGrad

Info about ASA is here:
    http://users.clas.ufl.edu/hager/papers/CG/
    W. W. Hager and H. Zhang,
    "A new active set algorithm for box constrained optimization", 
    SIAM Journal on Optimization, 17 (2006), pp. 526-557
    http://users.clas.ufl.edu/hager/papers/CG/asa.pdf

source code for ASA version 2.2 (April 14, 2011) is here:
http://users.clas.ufl.edu/hager/papers/CG/Archive/ASA_CG-2.2.tar.gz

This file will download ASA, install it, and compile
the mex files. We also show a simple example
(the same example as driver1.c in the ASA package, except
now in Matlab), as well as a non-negative least-squares
problem example.

Stephen.Becker@Colorado.edu

%}

% download the code
USE_NEW_VERSION = true;
if USE_NEW_VERSION
    % Code was developed for ver 2.2, but ver 3.0 was released in 2013 (main change is using an updated CG subroutine)
    % so you can try your luck with that too:
    url ='http://users.clas.ufl.edu/hager/papers/CG/Archive/ASA_CG-3.0.tar.gz';
else
    url ='http://users.clas.ufl.edu/hager/papers/CG/Archive/ASA_CG-2.2.tar.gz';
end

srcDir = untar(url,'./');

if isunix && exist('modify_printf.sh','file')
    % This will modify the ASA .c file to allow it
    % to print to the Matlab screen (it replaces
    % "printf" with "mexPrintf")
    unix('bash ./modify_printf.sh');
    % (Update: I don't think Mathworks let's us upload .sh files)
    % The file should be:
%{
    
#!/usr/bin/env bash

# the idea: replace all "printf" statements with "Printf"
#   Then, if we compile with the -DMEXPRINTF option,
#   "Printf" is linked to "mexPrintf" so that output
#   appears on the Matlab prompt.

# Stephen Becker, March 22 2012; Jan 22 2018

#file='ASA_CG-2.2/asa_cg.c' # change appropriately
file='ASA_CG-3.0/asa_cg.c'


cat > tempFile << EOF
#ifdef MEXPRINTF
#define Printf mexPrintf
#else
#define Printf printf
#endif
EOF

sed 's/printf/Printf/g' $file >> tempFile

mv tempFile $file    
    
    
%}
else
    % If you don't have linux/unix, just edit asa_cg.c
    %   by hand and do a search/replace with a text
    %   editor, and add these lines to the top of the file:
%{
#ifdef MEXPRINTF
#define Printf mexPrintf
#else
#define Printf printf
#endif
%}
    
    % This isn't necessary, but it's nice to have output
    % on the screen.
end
%%
% compile with
if USE_NEW_VERSION
    % Define the variable VER30 to make it work:
    mex -DMEXPRINTF -DVER30  asa_wrapper.c ASA_CG-3.0/asa_cg.c -IASA_CG-3.0 -largeArrayDims
else
    mex -DMEXPRINTF asa_wrapper.c ASA_CG-2.2/asa_cg.c -IASA_CG-2.2 -largeArrayDims
end

%% Run a sample problem
n = 100;
lo = zeros(n,1);    % lower bound
hi = ones(n,1);     % upper bound
x  = ones(n,1);     % initial guess; necessary so that it knows the size of the problem

% To use this on your own, make these strings point
%   to your function. The "fcn" is a function that computes
%   your objective function; "grad" computes its gradient.
%   The "fcnGrad" is optional; if provided,
%   it computes both the objective and the gradient
%   (since sometimes this saves computation over computing
%    them separately). Its outputs should be [f,g] (f=objective,
%    g=gradient) in that order).

fcn  = 'asa_fcn';
grad = 'asa_grad'; 
fcnGrad = 'asa_fcnGrad'; % optional

% add some ASA options (these are optional). See driver1.c for examples,
%   and see asa_user.h for all possible values
opts = [];
opts.PrintParms = false;
CGopts = struct('PrintParms',false);


% run the function
[out,status] = asa_wrapper( x, lo, hi,fcn,grad, fcnGrad, opts, CGopts);

%{
Your display should look like this

Using combined fcn_grad function 'asa_fcnGrad'

For code ver 2.2:

Final convergence status = 0
Convergence tolerance for gradient satisfied
projected gradient max norm:  3.560360e-09
function value:              -4.011035e+02

Total cg  iterations:                   11
Total cg  function evaluations:         17
Total cg  gradient evaluations:         16
Total cbb iterations:                    3
Total cbb function evaluations:          4
Total cbb gradient evaluations:          4
------------------------------------------
Total function evaluations:             21
Total gradient evaluations:             20
==========================================

For code ver 3.0:

Final convergence status = 0
Convergence tolerance for gradient satisfied
projected gradient max norm:  2.618104e-10
function value:              -4.011035e+02

Total cg  iterations:                    8
Total cg  function evaluations:          9
Total cg  gradient evaluations:         15
Total cbb iterations:                    3
Total cbb function evaluations:          4
Total cbb gradient evaluations:          4
------------------------------------------
Total function evaluations:             13
Total gradient evaluations:             19
==========================================

%}



%% Test on a non-negative least-squares (NNLS) problem
%{
The NNLS problem:

min_x  .5*||Ax-b||^2  = .5x'A'Ax - x'A'b + .5b'b
subject to
 x >= 0

%}
n = 100;
lo = zeros(n,1);    % lower bound
hi = inf(n,1);     % upper bound
x  = ones(n,1);     % initial guess; necessary so that it knows the size of the problem
A = randn(n,n); %A = A*A'; % make it pos def.
b = randn(n,1);
%% ... first, verify via lsqnonneg
% lsqnonneg is very slow for large problems, but we can
% use it to help verify that we don't have any huge bugs:
x_lsq = lsqnonneg( A, b );
%% ... second, run with ASA
% I have created these 3 special functions just for quadratic problems.
% You don't need to modify every time you change your quadratic
% parameters; rather, pass in the parameters via the "param" structure.
fcn  = 'asa_quadratic_fcn';
grad = 'asa_quadratic_grad'; 
fcnGrad = 'asa_quadratic_fcnGrad'; % optional
param = struct('A',A,'b',b);
% an alternative way:
% param = struct('Q',A'*A,'c',-A'*b,'offset',norm(b)^2/2);


% add some options (these are optional). See driver1.c for examples,
%   and see asa_user.h for all possible values
[opts,CGopts] = deal(struct('PrintParms',false));

% zero-out the counters
asa_quadratic_fcnGrad();


% run the function
[out,status,statistics] = asa_wrapper( x, lo, hi,fcn,grad, fcnGrad, opts, CGopts, param);
fprintf('Error compared to LSQNONNEG version is %.2e\n', ...
    norm( out - x_lsq) );

% View the function values
fcnHistory = asa_quadratic_fcnGrad();
semilogy( fcnHistory - min(fcnHistory), 'o-' );
xlabel('call to the objective function');
ylabel('value of objective function (minus true answer)');
