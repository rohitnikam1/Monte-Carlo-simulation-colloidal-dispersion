clear all 
close all 
clc

L = 170 ;

load r_colloid.txt 
r = r_colloid ;
n = 2040 ; % total number of ions 
h = zeros(n,3);

for j = 1 : n
	x = 0.5*L*( 2 * rand() -1 ) ;
	y = 0.5*L*( 2 * rand() -1 ) ;
	z = 0.5*L*( 2 * rand() -1 ) ;
	h(j,:) = [x y z] ;

end

a = size(r,1) ; 

for i = 1 : a
	for k = 1 : n

		rxij = r(i,1) - h(k,1) ;
		ryij = r(i,2) - h(k,2) ;
		rzij = r(i,3) - h(k,3) ;

		rxij = rxij - round( rxij / L ) * L ; 	
		ryij = ryij - round( ryij / L ) * L ; 
		rzij = rzij - round( rzij / L ) * L ;

		rij2 = rxij^2 +	ryij^2 + rzij^2 ; 
		rij = sqrt(rij2);

		if rij <= 10 + 0.1
			h(k,:) = [0 0 0] ;
		end
	end
end


k = 1 ;
for i = 1 : n
	if h(i,2) ~= 0
		r_ions(k,:) = h(i,:);
		k = k + 1 ;
	end
end
