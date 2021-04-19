%%******************* Real space sum for whole system *******************%%

% Calculates the real part of the Ewald sum
% q	:	charge vector containing charge on each species
% eps :	dielectric constant
% rx, ry, rz	:	vectors containing positions of species
% L	:	length of side of cubic central box
% n	:	number of species
% sigma_gauss	:	width of Gaussian
% sigma :	sigma for macro - macro van der Waal interactions
% epsilon :	epsilon for macro - macro van der Waal interactions


function [ pot_real ] = ewald_real_new( q, eps, k, rx, ry, rz, L, n_neg, n_pos, n_macro, kappa_gauss, rcut2_ewald, R_macro, R_micro)

	u = 0.0 ;

	for i = 1 : n_neg + n_pos + n_macro - 1 
		rxi = rx(i) ;
		ryi = ry(i) ; 
		rzi = rz(i) ;

		for j = i + 1 : n_neg + n_pos + n_macro
			rxij = rxi - rx(j)
			ryij = ryi - ry(j)
			rzij = rzi - rz(j)
	
			rxij = rxij - round( rxij / L )*L
			ryij = ryij - round( ryij / L )*L
			rzij = rzij - round( rzij / L )*L

	rij2 = rxij^2 +	ryij^2 + rzij^2 ; 

		if rij2 <= rcut2_ewald
			rij = sqrt( rij2 ) ;
			krij = kappa_gauss * rij ;

			if i <= n_neg
				if j <= n_neg
					if rij < 2 * R_micro
						pot_real = inf ; 
						return
					end

					u = u + q(1) * q(1) * erfc( krij ) / rij ; %% interaction between two negative ions ( electrostatic )

				elseif j <= n_neg + n_pos
					if rij < 2 * R_micro
						pot_real = inf ;
						return
					end

					u = u + q(2) * q(1) * erfc( krij ) / rij ; %% interaction between negative and positive ions ( electrostatic )

				else
					if rij < ( R_macro + R_micro )
						pot_real = inf ;
						return
					end

					u = u + q(3) * q(1) * erfc( krij ) / rij ; %% interaction between colloidal particle and negative ion ( electrostatic )

				end

			elseif i <= n_neg + n_pos
				if j <= n_neg + n_pos
					if rij < 2 * R_micro
						pot_real = inf ;
						return
					end

					u = u + q(2) * q(2) * erfc( krij ) / rij ; %% interaction between two positive ions ( electrostatic )

				else
					if rij < ( R_macro + R_micro )
						pot_real = inf ;
						return
					end

					u = u + q(3) * q(2) * erfc( krij ) / rij ; %% interaction between colloidal particle and positive ion ( electrostatic )
				end
			else
				if rij < ( R_macro + R_micro )
					pot_real = inf ;
					return
				end

				u = u + q(3) * q(3) * erfc( krij ) / rij ; %% interaction between two colloidal particles ( electrostatic + van der Waals )
			end
		end
	end

	pot_real = k * u / eps ;

end
