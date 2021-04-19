% K: the matrix containing all three components of each wave vector
% kvector : the vector containing the term exp(-ksqr)/ksqr for each wavevector
% size_kvector = (2*nmax + 1) * (2*nmax + 1) * (2*nmax + 1) - 1
% L: length of box decided by number of secies and density
% kappa_gauss controls the width of Gaussian distribution

function [ kvector K ] = setup_kvector( kappa_gauss, L )

	nmax = 4 ;	%% max integer component of k-vector
	b = 1.0 / 4.0 / kappa_gauss / kappa_gauss ;	
	count = 0 ;	%% total number of k-vectors stored

	for nx = -nmax : nmax
		kx = 2 * pi * nx / L ;
		for ny = -nmax : nmax
			ky = 2 * pi * ny / L ;
			for nz = -nmax : nmax
				kz = 2 * pi * nz / L ;
				n2 = nx * nx + ny * ny + nz * nz ;
				if n2 ~= 0
					count	= count + 1 ;
					K( count, : ) = [ kx	ky	kz ] ;
					k2 =	kx * kx	+	ky * ky	+	kz * kz ;
					kvector( count ) = exp( -b * k2 ) / k2 ;
				end
			end
		end
	end
	fprintf(' Number of wave vectors is %d.\n',count)
end
