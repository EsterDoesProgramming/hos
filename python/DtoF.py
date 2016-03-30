import numpy as np

def sech(x):
    return 1/np.cosh(x)

def DtoF( hetaD, kx, ky, **kwargs ):
#   Transforms directional spectrum into Fourier spectrum.
#
#   hetaF, hphiF = DtoF( hetaD, kx, ky, g=.. );

    # Set g=1 unless given as input
    depth_flag = 0
    if 'g' in kwargs: depth_flag = 1

    g  = kwargs.get('g', 9.81)
    h  = kwargs.get('h',1000.0)
        
    Nx = np.size(kx)
    Ny = np.size(ky)
    
    hetaF = np.zeros((Nx,Ny))
    hphiF = np.zeros((Nx,Ny))
    
    # A loop is highly inefficient, but in this case we prefer to 
    # make the code human readable.
    for ikx in np.arange(0,Nx):

        for iky in np.arange(0,Ny):

            # Mirrored indexes
            ikxm = Nx - ikx
            ikym = Ny - iky
    
            # Special threatment for the Nyquist
            # frequency. We throw away its
            # contribution, which should be
            # nothing but small noise.
            if ((ikx == 0 ) | (iky == 0)): #0=1

                hetaD[ikx,iky] = 0
                ikxm = ikx
                ikym = iky
   
            k = np.abs( kx[ikx] + 1j*ky[iky] )
    
            if (depth_flag == 0):
                norm  = np.sqrt(g**(1/2)/2/k**(3/2))
                alpha = np.sqrt(g/k)
                
            elif (depth_flag == 1):
                kh    = k*h
                norm  = np.sqrt(0.5*np.sqrt(g)*(np.sqrt(np.tanh(kh)/k**3) + h*(1/np.sqrt(k*np.tanh(kh)))*sech(kh)**2))
                alpha = np.sqrt(g/k/np.tanh(kh))
                
            hetaF[ikx,iky] = 0.5*norm*hetaD[ikx,iky] + 0.5*norm*np.conj(hetaD[ikxm, ikym])
            hphiF[ikx,iky] = -1j*0.5*norm*alpha*hetaD[ikx,iky] + 1j*0.5*alpha*norm*np.conj(hetaD[ikxm, ikym])

            if ((kx[ikx] == 0 ) & (ky[iky] == 0)):
                hetaF[ikx,iky] = 0.
                hphiF[ikx,iky] = 0.


    return hetaF, hphiF
