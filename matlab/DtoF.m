function [ hetaF, hphiF ] = DtoF( hetaD, kx, ky, g )
%   Transforms directional spectrum into Fourier spectrum.
%
%   [ hetaF, hphiF ] = DtoF( hetaD, kx, ky, g );

    % Set g=1 unless given as input
    if (nargin==3)
	g=1;
    end
    
    Nx = length(kx);
    Ny = length(ky);
    
    % A loop is highly inefficient, but in this case we prefer to 
    % make the code human readable.
    for ikx = 1:length(kx)

	for iky = 1:length(ky)

	    % Mirrored indexes
	    ikxm = Nx - ikx + 2;
	    ikym = Ny - iky + 2;
	    
	    % Special threatment for the Nyquist
	    % frequency. We throw away its
	    % contribution, which should be
	    % nothing but small noise.
	    if (ikx==1 || iky==1)
		
	        hetaD(ikx,iky) = 0;
		ikxm = ikx;
	        ikym = iky;
	    
	    end
	   
	    k = abs( kx(ikx) + 1i*ky(iky) );
	    alpha = sqrt(g/k);
	    
	    hetaF(ikx,iky) = 0.5*hetaD(ikx,iky) + 0.5*conj(hetaD(ikxm, ikym));
	    hphiF(ikx,iky) = -1i*0.5*alpha*hetaD(ikx,iky) + 1i*0.5*alpha*conj(hetaD(ikxm, ikym));

	    % STILL NEED TO NORMALIZE !!!
	    norm = sqrt( g^2/2/(g*k)^(3/2) );
	    hetaF(ikx,iky) = norm*hetaF(ikx,iky); 
	    hphiF(ikx,iky) = norm*hphiF(ikx,iky);
	    
	    
	    if (kx(ikx)==0 & ky(iky)==0)
	    
		hetaF(ikx,iky) = 0;
		hphiF(ikx,iky) = 0;
		
	    end
	    
	end
    
    end

    
    
end

