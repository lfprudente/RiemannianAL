function checkd(n,m,xini,l,u)

macheps12 = 10^(-8);

rng(123456);

for i = 1:n
    if ( l(i) < xini(i) && xini(i) < u(i) )
        x(i) = xini(i) + macheps12 * ( 2.0d0 * rand - 1.0d0 ) * max( 1.0d0, abs( xini(i) ) );
    elseif ( xini(i) == l(i) )
        x(i) = xini(i) + macheps12 * rand * max( 1.0d0, abs( xini(i) ) );
    else
        x(i) = xini(i) - macheps12 * rand * max( 1.0d0, abs( xini(i) ) );
    end
    x(i) = max( l(i), min( x(i), u(i) ) );
end
x = x';

fprintf(' Derivatives will be tested at the perturbed initial guess:\n')
for i = 1:n
    fprintf(' x( %6i ) = %15.8f \n',i,x(i))
end

% CHECK OBJECTIVE FUNCTION GRADIENT

fprintf('\n Would you like to check subroutine evalf?\n')
prompt = " Type Y(es), N(o), A(bort checking) or S(kip checking this subroutine): ";
answer = input(prompt,"s");
if isempty(answer)
    answer = 'S';
end

if ( answer == 'A' || answer == 'a' )
    return
elseif ( answer == 'Y' || answer == 'y' )
    checkg(n,x)
end

% CHECK CONSTRAINT GRADIENTS

if ( m == 0 ) 
    fprintf('\n')
    fprintf('Press any key to continue.\n')
    fprintf('\n\n')
    pause
    return
end

fprintf('\n Would you like to check subroutine evalnc?\n')
prompt = " Type Y(es), N(o), A(bort checking) or S(kip checking this subroutine): ";
answer = input(prompt,"s");
if isempty(answer)
    answer = 'S';
end

if ( answer == 'A' || answer == 'a' )
    return
elseif ( answer == 'Y' || answer == 'y' )
    for i = 1:m
        fprintf('\n Check gradient of constraint %5i ?\n',i)
        prompt = " Type Y(es), N(o), A(bort checking) or S(kip checking this subroutine): ";
        answerc = input(prompt,"s");
        if isempty(answer)
            answerc = 'S';
        end
        
        if ( answerc == 'A' || answerc == 'a' )
            return
        elseif ( answerc == 'Y' || answerc == 'y' )
            checknc(n,x,i)
        end
    end
end

fprintf('\n')
fprintf('Press any key to continue.\n')
fprintf('\n\n')
pause

end

% ==================================================================

function checkg(n,x)

[g,flag] = evalg(n,x);

fprintf('\n')
fprintf(' Gradient vector of the objective function.\n')
fprintf(' Index             evalg  Central diff (two different steps)    Absolute error\n')


maxerr = 0.0d0;

eps = 10^(-16/3);

for i = 1:n
	 tmp  = x(i);

	 step1 = eps * max( abs( tmp ), 1.0d0 );

	 x(i) = tmp + step1;
     [fplus,flag] = evalf(n,x);
	 
	 x(i) = tmp - step1;
     [fminus,flag] = evalf(n,x);

	 gdiff1 = ( fplus - fminus ) / ( 2.0d0 * step1 );

	 step2 = eps * max( abs( tmp ), 1.0d-03 );

	 x(i) = tmp + step2;
     [fplus,flag] = evalf(n,x);

	 x(i) = tmp - step2;
	 [fminus,flag] = evalf(n,x);
	 
	 x(i) = tmp;

	 gdiff2 = ( fplus - fminus ) / ( 2.0d0 * step2 );

	 tmp = min( abs( g(i) - gdiff1 ), abs( g(i) - gdiff2 ) );
     fprintf(' %5i   %15.8f   %15.8f   %15.8f   %15.8f\n',i,g(i),gdiff1,gdiff2,tmp)

	 maxerr = max( maxerr, tmp );

end

fprintf('\n Maximum absolute error = %15.8f \n',maxerr)

end

% ==================================================================

function checknc(n,x,ind)

[g,flag] = evalnc(n,x,ind);

fprintf('\n')
fprintf(' Gradient vector of the objective function.\n')
fprintf(' Index             evalg  Central diff (two different steps)    Absolute error\n')


maxerr = 0.0d0;

eps = 10^(-16/3);

for i = 1:n
	 tmp  = x(i);

	 step1 = eps * max( abs( tmp ), 1.0d0 );

	 x(i) = tmp + step1;
     [fplus,flag] = evalcc(n,x,ind);
	 
	 x(i) = tmp - step1;
     [fminus,flag] = evalcc(n,x,ind);

	 gdiff1 = ( fplus - fminus ) / ( 2.0d0 * step1 );

	 step2 = eps * max( abs( tmp ), 1.0d-03 );

	 x(i) = tmp + step2;
     [fplus,flag] = evalcc(n,x,ind);

	 x(i) = tmp - step2;
     [fminus,flag] = evalcc(n,x,ind);
	 
	 x(i) = tmp;

	 gdiff2 = ( fplus - fminus ) / ( 2.0d0 * step2 );

	 tmp = min( abs( g(i) - gdiff1 ), abs( g(i) - gdiff2 ) );
     fprintf(' %5i   %15.8f   %15.8f   %15.8f   %15.8f\n',i,g(i),gdiff1,gdiff2,tmp)

	 maxerr = max( maxerr, tmp );

end

fprintf('\n Maximum absolute error = %15.8f \n',maxerr)

end