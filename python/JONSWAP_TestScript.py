# Test for the JONSWAP spectrum

g     = 9.8               #acceleration due to gravity
alpha = 1.0 #0.125

#-- Discretization: Number of cells in x and y direction
Nx = alpha*512           #fourier modes in x/y directions
Ny = alpha*512 

kp    = 0.033;                #for deep water
lmbda = 2*np.pi/kp;
Lx    = 16*lmbda              # Lx = 1*ceil(((epsilon^-2)*pi/kp)); Lx  = 16*lambda;4.*np.pi
Ly    = Lx                   # Ly = 4.*np.pi/sqrt(3.)

kx = 2*np.pi/Lx*np.arange(-Nx/2,Nx/2)
ky = 2*np.pi/Ly*np.arange(-Ny/2,Ny/2)
 
etar, phir = create_jonswap(kx,ky)

fig = plt.figure(1)
#plt.contourf(kx,ky,etar)
plt.imshow(etar)
fig = plt.figure(2)
plt.contourf(kx,ky,phir)
plt.colorbar()
plt.show()
