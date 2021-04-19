%%**************** Reciprocal space sum *******************%%
% Calculates the reciprocal part of the Ewald sum
% q:		charge vector containing charge on each species
% eps:		dielectric constant of background medium
% rx, ry, rz :	vectors containing positions of species
% L :		length of side of cubic central box
% kappa_gauss : width of Gaussian
% k :		1 / 4 / pi / epsilon zero
% kvector :	parameter storing the term 2*pi*exp(-k_sqr/4kappa_sqr)/ksqr
% eik : 	vector containing exp(ik.(rj)) for each species j
% K :		all components of each wave vector kx ky kz
% kvector :	vector containing exp(-b*k^2)/k^2 for each wave vector

function [ pot_resip ] = ewald_reciprocal_new( q, eps, k, rx, ry, rz, eik, L, n_neg, n_pos, n, kappa_gauss, K, kvector )

	pot_resip = 0.0 ;
	
	%****** sum over all k-vectors ******%

	rows = size( K, 1 ) ;	% rows of matrix K
	
	for i = 1 : rows	%% loop for each k
		sum = 0 ;
		for j = 1 : n
			rx( j ) = rx( j ) - round( rx( j ) / L ) * L ;
			ry( j ) = ry( j ) - round( ry( j ) / L ) * L ;
			rz( j ) = rz( j ) - round( rz( j ) / L ) * L ;


			eik( j ) = exp( 1i*(K(i, 1)*rx(j) + K(i, 2) * ry(j) + K(i, 3)*rz(j) )) ;
			
			if j <= n_neg
				add = q(1) * eik ( j ) ; 
			elseif j <= n_neg + n_pos
				add = q(2) * eik ( j ) ; 
			else
				add = q(3) * eik ( j ) ;
			end

			sum = sum + add ;
		end
		
		pot_resip = pot_resip + sum * conj( sum ) * kvector( i );
	end

	pot_resip = 2 * pi * k * pot_resip / L / L / L / eps ;


	%********** Self part of the reciprocal space sum ********** %% pot_self = 0 ;
	
	for temp = 1 : n

		if temp <= n_neg
			add = q(1) * q(1) ; 
		elseif temp <= n_neg + n_pos
			add = q(2) * q(2) ; 
		else
			add = q(3) * q(3) ;

		end

		pot_self = pot_self + add ;

	end

	pot_self = (pi^-0.5) * k * kappa_gauss / eps * pot_self	;


	%********** Total reciprocal space energy ************ % 
	pot_resip =	pot_resip -	pot_self ;

end
