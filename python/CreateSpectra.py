import numpy as np
import DtoF as DF
import JONSWAP as J
import OmegaTheta as OT


def create_gauss():
    return 0

def create_jonswap(kx,ky):
    #-- JONSWAP Spectrum
    g     = 9.8               #acceleration due to gravity

    #-- Discretization: Number of cells in x and y direction
    Nx = len(kx)
    Ny = len(ky)

    #-- JONSWAP PARAMETERS
    gamma   = 3.3
    alpha_p = 0.081
    dir_spr = 40.*np.pi/180.

    kp      = 0.033            #for deep water
    lmbda   = 2*np.pi/kp
    omega_p = np.sqrt(g*kp)
    epsilon = 0.12
    theta_p = 0.0 
    h       = 50               #water depth... set high to induce deep water effects
    sigma   = epsilon/2/kp     #desired eta variance
    hs      = 4*sigma
    
    omega, theta = OT.OmegaTheta(kx, ky, g=g, h=h)
    
    Pk, Sj, D = J.JONSWAP(omega, theta, g=g, alpha_p=alpha_p, omega_p=omega_p, gamma=gamma, \
                          ThetaPar=dir_spr, theta_shift=theta_p);

    #Introduce random phase 
    hetaD = np.multiply(np.sqrt(Pk),np.exp(np.multiply(1j*2*np.pi,np.random.rand(Nx,Ny))))

    #Convert to wavenumber space
    hetar, hphir   = DF.DtoF(hetaD, kx, ky, g=g, h=h)
    hetar          = hetar.T
    hphir          = hphir.T
    
    # Compute wave field in the physical space
    #NOTE need to use fourier shift as hetar was not generated using the fft
    etar = np.fft.ifft2(np.fft.ifftshift(hetar))    
    etar = etar.real
    phir = np.fft.ifft2(np.fft.ifftshift(hphir))
    phir = phir.real
    
    # Rescale the wave field to the desired eta variance
    rescale = np.std(etar)/(hs/4) 
    etar    = etar/rescale
    phir    = phir/rescale

    return etar, phir


def create_stokes(kx,ky):
    return 0