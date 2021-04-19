% Calculates the total van der Waals interaction energy using LJ pair potential
% n : 		total number of species
% eps : 	maximum attractive interaction energy between two nanoparticles
% sig : 	second LJ parameter
% rx, ry, rz :	vectors containing positions of species
% L :		length of side of cubic central box
% rcut2 :	square of cut-off distance

function [pot_total] = vanderwaals_total( rx, ry, rz, eps, sig, rcut2, L, n)
	
	pot_total = 0.0;

	for i = 1 : n
		rxi = rx( i ) ;
		ryi = ry( i ) ;
		rzi = rz( i ) ;

		for j = i + 1 : n
			rxij = rxi - rx( j ) ;
			ryij = ryi - ry( j ) ;
			rzij = rzi - rz( j ) ;

			rxij = rxij - round( rxij / L ) * L ;
			ryij = ryij - round( ryij / L ) * L ;
			rzij = rzij - round( rzij / L ) * L ;

			rij2 = rxij^2 + ryij^2 + rzij^2 ;

			if rij2 <= rcut2
				sr = ( sig^2 / rij2 )^3 ;
				pot_total = pot_total + sr * ( sr - 1 ) ;
			end
		end
	end

	pot_total = 4 * pot_total * eps ;
end
