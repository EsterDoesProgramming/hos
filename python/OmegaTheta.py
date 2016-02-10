import numpy as np

def  OmegaTheta( kx, ky, **kwargs ):
# Returns polar coordinates on each node of the
# rectangular k-grid in two matrices.
#
# omega, theta = OmegaTheta( kx, ky, g=.. );
#
# NOTE: since we use meshgrid() we obtain
# size(omega) = [ length(ky) length(kx) ]
# and likewise for theta.

    # Set g=1 unless given as input
    g = kwargs.get('g', 1)    

    [KX,KY] = np.meshgrid(kx,ky)
    
    omega = np.sqrt(g*abs(KX + 1j*KY))
    theta = np.angle(KX + 1j*KY)

    # Enforce 0 < theta <= 2*pi
    for k in np.arange(len(kx)):
        for l in np.arange(len(ky)):
            theta[k,l] = theta[k,l] + 2*np.pi*theta[k,l] if (theta[k,l] < 0) else theta[k,l]

    return omega, theta

