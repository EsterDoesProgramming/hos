import numpy as np

def theta1(theta):
    if (theta < 0.0):
         return theta + 2*np.pi
    else: return theta
    
fetchtheta = np.vectorize(theta1) 

def  OmegaTheta( kx, ky, **kwargs ):
# Returns polar coordinates on each node of the
# rectangular k-grid in two matrices.
#
# omega, theta = OmegaTheta( kx, ky, g=.., h=.. );
#
# NOTE: since we use meshgrid() we obtain
# size(omega) = [ length(ky) length(kx) ]
# and likewise for theta.

    depth_flag = 0  #1 for finite depth, 0 for deep water 
    if 'g' in kwargs: depth_flag = 1
        
    print 'depth_flag', depth_flag

    # Set g=1 unless given as input
    g = kwargs.get('g', 1) 
    h = kwargs.get('h', 1000.0)
   
    KX,KY = np.meshgrid(kx,ky)
    theta   = np.zeros((len(kx),len(ky)))
    omega   = np.zeros((len(kx),len(ky)))
    
    K = KX + 1j*KY
    
    gK = np.multiply(g,np.abs(K))

    if (depth_flag == 0): omega = np.sqrt(gK)
    if (depth_flag == 1): omega = np.sqrt(np.multiply(gK,np.tanh(np.multiply(np.abs(K),h))))
    
    theta = np.angle(K)
    theta = fetchtheta(theta)

    return omega, theta

