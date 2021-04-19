% Calculates the total van der Waals interaction energy of a single nanoparticle using	LJ pair potential
% p	:	index of a particle
% n	:	total number of species
% eps :	maximum attractive interaction energy between two nanoparticles
% sig :	second LJ parameter
% rx, ry, rz	:	vectors containing positions of species
% L	:	length of side of cubic central box
% rcut2	:	square of cut-off distance


function [pot_single] = vanderwaals_single( p, rx, ry, rz, eps, sig, rcut2, L, n )

	pot_single = 0.0;

	rxi = rx( p ) ; 
	ryi = ry( p ) ; 
	rzi = rz( p ) ;

	for j = 1 : p - 1
		rxij = rxi - rx( j ) ; 
		ryij = ryi - ry( j ) ; 
		rzij = rzi - rz( j ) ;

		rxij = rxij - round( rxij / L ) * L ; 
		ryij = ryij - round( ryij / L ) * L ; 
		rzij = rzij - round( rzij / L ) * L ;

		rij2 = rxij^2 + ryij^2 + rzij^2 ; 

		if rij2 <= rcut2
			sr = ( sig^2 / rij2 )^3 ;
			pot_single = pot_single + sr * ( sr - 1 ) ;

		end
	end

	for j = p+1:n

		rxij= rxi-rx(j); 
		ryij= ryi-ry(j); 
		rzij= rzi-rz(j);

		rxij = rxij - round( rxij / L ) * L ; 
		ryij = ryij - round( ryij / L ) * L ;
		rzij = rzij - round( rzij / L ) * L ; 

		rij2 = rxij^2 + ryij^2 + rzij^2 ;

		if rij2 <= rcut2
			sr = ( sig^2 / rij2 )^3 ;
			pot_single = pot_single +	sr * ( sr - 1 ) ;

		end
	end

	pot_single = 4 * pot_single * eps ;

end	
