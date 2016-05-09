import numpy as np
import DtoF as DF
import JONSWAP as J
import OmegaTheta as OT
import rieneckerfenton as RF


def create_gauss(kx,ky,**kwargs):
    return 0

def create_jonswap(kx,ky,**kwargs):

    #-- The Jonswap Parameters
    Param = kwargs.get('Param',None)

    Nx = len(kx)
    Ny = len(ky)

    omega, theta = OT.OmegaTheta(kx, ky, g=Param['g'], h=Param['h'])
    
    Pk, Sj, D = J.JONSWAP(omega, theta, g=Param['g'], alpha_p=Param['alpha_p'], omega_p=Param['omega_p'], \
                          gamma=Param['gamma'],ThetaPar= Param['dir_spr'], theta_shift=Param['theta_p']);

    #Introduce random phase 
    hetaD = np.multiply(np.sqrt(Pk),np.exp(np.multiply(1j*2*np.pi,np.random.rand(Nx,Ny))))

    #Convert to wavenumber space
    hetar, hphir   = DF.DtoF(hetaD, kx, ky, g=Param['g'], h=Param['h'])
    hetar          = hetar.T
    hphir          = hphir.T
    
    # Compute wave field in the physical space
    #NOTE need to use fourier shift as hetar was not generated using the fft
    etar = np.fft.ifft2(np.fft.ifftshift(hetar))    
    etar = etar.real
    phir = np.fft.ifft2(np.fft.ifftshift(hphir))
    phir = phir.real
    
    # Rescale the wave field to the desired eta variance
    rescale = np.std(etar)/(Param['hs']/4) 
    etar    = etar/rescale
    phir    = phir/rescale

    return etar, phir

def create_rieneckerfenton(kx,ky,**kwargs):
    #-- Use the parameters as in Rienecker & Fenton (1981)
    Param = kwargs.get('Param', None)
    
    eta, phi = RF.RieneckerFenton(kx,ky)

    return eta, phi

def create_stokes(kx,ky,**kwargs):
    return 0
