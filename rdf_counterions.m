%************ Radial Distribution Function ***************%


% n_neg		number of negative ions
% n_pos		number of positive ions
% r		probing distance vector
% dr		step size of r
% rho		number density of system
% n		total no. of species

function [ gr ] = rdffunction( rx, ry, rz, L, r, dr, n_neg, n_pos, n, rho )

	gr = zeros( numel(r), 1) ;

	for k = 1 : numel( r ) 
		rk = r( k ) ;
		NGRk = 0 ;

		for i = n_neg + n_pos + 1 : n 
			rxi = rx( i ) ;
			ryi = ry( i ) ; 
			rzi = rz( i ) ;

			ngri = 0 ;

			for j = n_neg + 1 : n_neg + n_pos 
				rxij = rxi - rx( j ) ;
				ryij = ryi - ry( j ) ; 
				rzij = rzi - rz( j ) ;

				rxij = rxij - round( rxij / L ) * L ; 
				ryij = ryij - round( ryij / L ) * L ; 
				rzij = rzij - round( rzij / L ) * L ;

				rij2 = rxij^2 +	ryij^2 + rzij^2 ; 
				roundij2 = 0.001 * round( 1000 * rij2 ) ; 
				rij = sqrt(roundij2) ;

				if rij < rk + dr && rij > rk – dr 
					ngri =	ngri + 1 ;
				end
			end

			NGRk = NGRk + ngri ;
		end

		NGRk = 0.25 * NGRk / n / pi / dr / rho ; 
		gr(k) = NGRk / (rk^2) ;
	end
end
