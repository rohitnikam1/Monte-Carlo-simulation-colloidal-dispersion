clear all 
clc

load r_init_500ions.txt 	%% initial configuration of system 
load R_500ions_5000moves.txt

r_init = r_init_500ions ; 
R = R_500ions_5000moves ;

rx = r_init( : , 1 ) ;
ry = r_init( : , 2 ) ;
rz = r_init( : , 3 ) ; 
n = numel( rx ) ;
L	= 170 ; 
dr = 0.1 ;
r	= 0 : dr : 0.5 * L ; 
n_neg = 500 ;
n_pos = 520 ;

n_loop = 5000 ;
n_sampling = 4900 ; % number from which the sampling starts 
nint	= 100 ;

sum = zeros( numel(r), 1 ); 
count = 0 ;

for i = 1 : n_sampling 
	p = R( i , 1 ) ;
	rx( p ) = R( i, 2 ) ;
	ry( p ) = R( i, 3 ) ;
	rz( p ) = R( i, 4 ) ;

end


for j = n_sampling + 1 : n_loop
	p = R( j, 1 ) ;
	rx( p ) = R( j, 2 ) ;
	ry( p ) = R( j, 3 ) ;
	rz( p ) = R( j, 4 ) ;

	if mod( j , nint ) == 0
		gr = rdffunction( rx, ry, rz, L, r, dr, n_neg, n_pos, n, rho ) ; 
		sum = sum + gr ;
		count = count + 1 ;
	end
end

avg_gr = sum / count ;
