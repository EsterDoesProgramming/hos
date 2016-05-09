import numpy as np

# Function implementing the JONSWAP directional wave spectrum
#
# P = JONSWAP(omega, theta, g=.., alpha_p=.., omega_p=.., gamma=.., ThetaPar=.., theta_shift=..);

def sigma1(omega, omega_p):
    if (omega <= omega_p):
        return 0.07 
    else: return 0.09 

sigma = np.vectorize(sigma1)

def theta1(theta):
    if (theta > np.pi):
         return theta - 2*np.pi
    else: return theta
    
fetchtheta = np.vectorize(theta1) 

def dirspread1(theta,ThetaPar):
    if (np.abs(theta) < 0.5*ThetaPar):
        return 2/ThetaPar*np.cos(np.pi*theta/ThetaPar)**2
    else: return 0.

dirspread = np.vectorize(dirspread1)
   
def JONSWAP(omega, theta, **kwargs):

    # default parameters
    g           = kwargs.get('g', 1.0)
    alpha_p     = kwargs.get('alpha_p', 2*np.pi*0.014/10) # 0.0081
    omega_p     = kwargs.get('omega_p', 2*np.pi)          # 1.0
    gamma       = kwargs.get('gamma', 6)
    ThetaPar    = kwargs.get('ThetaPar', 12*np.pi/180)
    theta_shift = kwargs.get('theta_shift', 0.0)
    
    # angular shift if needed
    theta = np.subtract(theta, theta_shift)
    theta = fetchtheta(theta)

    omegadiff = np.subtract(omega, omega_p)
    sigmavec  = sigma(omega,omega_p)

    delta = np.exp( np.divide(omegadiff**2, np.multiply(np.power(sigmavec,2),-2.*omega_p**2)) )
    const = np.multiply(alpha_p*g**2,np.power(omega,-5))
    
    # Filter Nyquist frequency - it's only small noise
    N = np.size(omega,1)
    const[N/2,N/2] = 0.0

    exponent = np.multiply(-1.25,np.power(np.divide(omega,omega_p),-4))

    # JONSWAP spectrum (Xiao et al. 2013: eq.(2.4))
    Sj = np.multiply(np.multiply(const,np.exp(exponent)),np.power(gamma, delta))
    
    # The directional spreading function (Xiao et al. 2013: eq.(2.5))
    D = dirspread(theta,ThetaPar)

    # Directional spectrum
    Pk = np.multiply(Sj,D)

    return Pk, Sj, D

