function [rhoini] = comprhoini(c,f,E,I)

rhoinimin = 10^(-8);
rhoinimax = 10^8;

sumc = 0.5 * ( sum(c(E).^2) + sum(max(c(I),0).^2) );

rhoini = 10 * max( 1, abs( f ) ) / max( 1, sumc );

rhoini = max( rhoinimin, min( rhoini, rhoinimax ) );