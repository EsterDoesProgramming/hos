import numpy as np

def sech(x):
    return np.cosh(x)**(-1)

def FtoD(hetaF, hphiF, kx, ky, **kwargs):
#   Transforms Fourier coeffs into directional spectrum coeffs.
#
#   hetaD = FtoD( hetaF, hphiF, kx, ky, g=.. );
#
#   NOTE: the Fourier coeffs are assumed to be
#   ordered with the zero-frequency component in
#   the center of the spectrum. If hu is the output of 
#   fft2, i.e., hu = fft2(u), it is necessary to
#   perform fftshift on hu before using it as input. 

    # Set g=1 unless given as input
    g = kwargs.get('g',1.0)

    Nx = len(kx)
    Ny = len(ky)
    
    # A loop is highly inefficient, but in this case we prefer to 
    # make the code human readable.

    hetaD = np.zeros((Nx,Ny))

    for ikx in np.arange(1,Nx+1):

	for iky in np.arange(1,Ny+1):
	   
	    k     = abs( kx[ikx-1] + 1j*ky[iky-1])
	    alpha = np.sqrt(g/k)
            h     = 1000.0
            kh    = k*h;
            jac   = np.sqrt(0.5*np.sqrt(g)*(np.sqrt(np.tanh(kh)/k**3) + h*(1/np.sqrt(k*np.tanh(kh)))*sech(kh)**2));
            #norm = 1
            alpha = np.sqrt(g/k/np.tanh(kh));


	    hetaD[ikx-1,iky-1] = 1j*hphiF[ikx-1,iky-1]/alpha[0] + hetaF[ikx-1,iky-1]/jac[0]

	#    norm               = np.sqrt( g**2/2/(g*k)**(3/2) )
	#    hetaF[ikx-1,iky-1] = hetaF[ikx-1,iky-1]/float(norm)
	#    hphiF[ikx-1,iky-1] = hphiF[ikx-1,iky-1]/float(norm)
	    
        
    return hetaD

