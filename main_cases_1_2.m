% Monte Carlo simulation of charged colloidal dispersion

%%%%%%**************** MAIN EXECUTION FILE *****************%%%%%%

% The following code calculates the total interaction energy of the system
% consisting of colloidal particles and ions of 1:1 ionic salt
% System description
% Cubic box of length 170 nm

% Species		Dimension (nm)		Number		Charge

% Colloidal particles		20		20		-2
% Counterions			0.1		40		+1
% Positive ions			0.1		1000		+1
% Negative ions			0.1		1000		-1
% Solvent (water)		-		continum	 0


% All length units in nanometers
% All energy units in kJ/mol

% R_macro	radius of each colloidal particle
% R_micro	radius of each ion ( same for cation and anion )
% n_macro	number of colloidal particles
% n_neg		number of negative ions
% n_pos		number of positive ions
% eps		dielectric constant of background medium ( water )
% delr_micro	displacement for ion
% delr_macro	displacement for colloidal particle
% epsilon,sigma	van der Waals interaction parameters
% rcut_lj	cutoff distance for van der Waals interaction
% L		length of side of cubic simulation box
% q		charge on each ion and colloidal particle
% n_loop	total Monte Carlo moves
% kappa_gauss	width of Gaussian distributions
% K		all components of each wave vector : kx, ky, kz
% eik		vector containing exp(ik.(rj)) for each species j
% kvector	vector containing exp(-b*k^2)/k^2 for each wave vector
% T		temperature
% R		vector recording the changes in the position after each move
% sigma_ewald	width of Gaussian
% k		1 / 4 / pi / epsilon zero


clc
clear all
load r_init_1000ions.txt	% contains initial configuration of system n_neg = 1000 ;
n_pos = 1040 ;
n_macro = 20 ;
n = n_neg + n_pos + n_macro ;

load K.mat
load kvector.txt
eik = zeros ( n , 1 ) ;

R_macro = 10.0 ; 
R_micro = 0.1 ; 
L = 170.0 ;


kappa_gauss = 5.714 / L ;
sigma_ewald = 1 / sqrt(2) / kappa_gauss ; 
eps = 80 ;
k = 1382.5 ;

q = [ -1 +1 -2 ] ;
rcut_ewald = 0.46 * L ; 
rcut2_ewald = ( rcut_ewald )^2 ;

rx = r_init_1000ions( : , 1 ) ; 
ry = r_init_1000ions( : , 2 ) ;
rz = r_init_1000ions( : , 3 ) ;

epsilon = 1 ;	% kJ/mol
sigma = 0.5 * ( 0.3 + 2 * R_macro ) ; % sigma is the van der Waals radius and equals half the internuclear distance between non-bonded particles

rcut_lj	=	3 * sigma ; rcut2_lj = ( rcut_lj )^2 ;

T	=	108 ;
beta	=	1.0 / ( 8.314472 * 0.001 * T ) ;

delr_micro = 1.5 ; 
delr_macro = 0.5 ;

n_loop = 5000 ;
R = zeros( n_loop , 4 ) ;


realpart = ewald_real_new( q, eps, k, rx, ry, rz, L, n_neg, n_pos, n_macro, kappa_gauss, rcut2_ewald, R_macro, R_micro ) ;

u_sys =	realpart + ewald_reciprocal_new( q, eps, k, rx, ry, rz, eik, L, n_neg, n_pos, n, kappa_gauss, K, kvector ) ;

upot = zeros( n_loop , 1 ) ;

naccept_micro = 0 ; 
nint_micro = 50; 
naccept_macro = 0 ; 
nint_macro = 50;

for nmc = 1 : n_loop

	% *********** preferential movement ************* % 
	g = rand() ;

	if g >= 0.2
		p = min( int32( rand() * (n_neg + n_pos) ) + 1	, n_neg + n_pos );	% an ion is chosen

		% ********** storing the old configuration ********** % 
		xold = rx( p ) ;
		yold = ry( p ) ; 
		zold = rz( p ) ;

		% ************* Giving a displacement **************** %
		rx( p ) = rx( p ) + delr_micro * ( 2.0 * rand() - 1.0 ) ;
		ry( p ) = ry( p ) + delr_micro * ( 2.0 * rand() - 1.0 ) ;
		rz( p ) = rz( p ) + delr_micro * ( 2.0 * rand() - 1.0 ) ;

		realpart = ewald_real_new( q, eps, k, rx, ry, rz, L, n_neg, n_pos, n_macro, kappa_gauss, rcut2_ewald, R_macro, R_micro ) ;

		if realpart == inf  
			rx( p )	= xold ; 
			ry( p )	= yold ; 
			rz( p )	= zold ;

			R(nmc,:) = [ p rx(p) ry(p) rz(p) ] ;
	else
		upnew_micro = realpart + ewald_reciprocal_new( q, eps, k, rx, ry, rz, eik, L, n_neg, n_pos, n, kappa_gauss, K, kvector ) ;

		delu_micro = upnew_micro - u_sys ;

		% ******************* Metropolis Algorithm ******************** % 
		if delu_micro <= 0
			u_sys =	upnew_micro ;
			naccept_micro =	naccept_micro +	1 ;
			R(nmc,:) = [ p rx(p) ry(p) rz(p) ] ; % recording the move

		else
			if exp( -beta * delu_micro ) > rand() 
				u_sys =	upnew_micro ;
				naccept_micro =	naccept_micro +	1 ;
				R(nmc,:) = [ p rx(p) ry(p) rz(p) ] ; % recording the move

			else
				rx( p ) = xold ;
				ry( p ) = yold ;
				rz( p ) = zold ;
				R(nmc,:) = [ p rx(p) ry(p) rz(p) ] ;
			end
		end
	end

	upot( nmc ) = u_sys + vanderwaals_total(rx(n - n_macro + 1 : n ), ry(n - n_macro + 1 : n ),rz(n - n_macro + 1 : n ), epsilon, sigma, rcut2_lj, L, n_macro ) ;

	% ******** Re-adjusting the maximum displacement ******** % 

	if mod( nmc , nint_micro ) == 0
		if naccept_micro / nint_micro >= 0.5
			delr_micro = min( 1.05 * delr_micro , 0.5 * L ) ; 
		else
			delr_micro = 0.95 * delr_micro ;

		end

		naccept	= 0 ;
	end

	else
		% if colloidal particle is chosen, then

		p = min( int32( rand() * n_macro ) + 1 + n - n_macro , n ) ;

		pot_colloid_old = vanderwaals_single( p - n + n_macro, rx(n - n_macro + 1 : n ), ry(n - n_macro + 1 : n ),rz(n - n_macro + 1 : n ), epsilon, sigma, rcut2_lj, L, n_macro);

		% ********** storing the old configuration ********** % 
		xold	=	rx( p ) ;
		yold	=	ry( p ) ; 
		zold	=	rz( p ) ;

		% ************* Giving a displacement**************** %
		rx( p ) = rx( p ) + delr_macro * ( 2.0 * rand() - 1.0 ) ;
		ry( p ) = ry( p ) + delr_macro * ( 2.0 * rand() - 1.0 ) ;
		rz( p ) = rz( p ) + delr_macro * ( 2.0 * rand() - 1.0 ) ;

		pot_colloid_new = vanderwaals_single( p - n + n_macro, rx(n - n_macro + 1 : n ), ry(n - n_macro + 1 : n ),rz(n - n_macro + 1 : n ), epsilon, sigma, rcut2_lj, L, n_macro);


		realpart = ewald_real_new( q, eps, k, rx, ry, rz, L, n_neg, n_pos, n_macro, kappa_gauss, rcut2_ewald, R_macro, R_micro ) ;

		if realpart == inf
			rx( p ) = xold ; 
			ry( p ) = yold ; 
			rz( p ) = zold ;

			R(nmc,:) = [ p rx(p) ry(p) rz(p) ] ;

		else
			upnew_macro =	realpart + ewald_reciprocal_new( q, eps, k, rx, ry, rz, eik, L, n_neg, n_pos, n, kappa_gauss, K, kvector ) ;

			delu_macro = upnew_macro - u_sys + pot_colloid_new - pot_colloid_old;


			% ******************* Metropolis ******************** % 

			if delu_macro <= 0
				u_sys =	upnew_macro ;
				naccept_macro =	naccept_macro +	1 ;

				R(nmc,:) = [ p rx(p) ry(p) rz(p) ] ; % recording the move

			else
				if exp( -beta * delu_macro ) > rand()
					u_sys =	upnew_macro ;
					naccept_macro =	naccept_macro +	1 ;

					R(nmc,:) = [ p rx(p) ry(p) rz(p) ] ; % recording move
				else
					rx( p ) = xold ; 
					ry( p ) = yold ; 
					rz( p ) = zold ;

					R(nmc,:) = [ p rx(p) ry(p) rz(p) ] ;
				end
			end
		end

		upot( nmc ) = u_sys + vanderwaals_total(rx(n - n_macro + 1 : n ), ry(n - n_macro + 1 : n ), rz(n - n_macro + 1 : n ), epsilon, sigma, rcut2_lj, L,n_macro) ;

		% ******* Re-adjusting the maximum displacement ********* %

		if mod( nmc , nint_macro ) == 0
			if naccept_macro / nint_macro >= 0.5
				delr_macro = min( 1.05 * delr_macro , 0.5 * L ) ;
			else
				delr_macro = 0.95 * delr_macro ;
			end
			naccept_macro =	0 ;

		end
	end
end
