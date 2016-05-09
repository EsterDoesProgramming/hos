import numpy
import matplotlib.pyplot as plt

##
# @file rieneckerfenton.py
# @author N. Beisiegel
# @date May 4, 2016
# @brief Provides a routine rieneckerfenton that computes the initial condition computed in
#        Rienecker & Fenton (1981); 2N+6 unknowns are solved with Newton's method.

#-------------------------------------------------#
# Auxiliary functions for velocities etc.
#

def compute_potential(x,y,B,k,D):
    N     = len(B)
    sizeX = len(x)
    sizeY = len(y)

    erg = numpy.zeros((sizeX,sizeY))

    for i in xrange(0,sizeX):
        for j in xrange(0,sizeY):
            erg[i,j] = B[0]*y[j]
            for m in xrange(1,N):
                erg[i,j] += B[m] * numpy.sinh(float(m)*k*y[j])/ \
                          numpy.cosh(float(m)*k*D) * numpy.cos(float(m)*k*x[i])
    return erg


def u(eta,B,D,k,N):

    y    = B[0] * numpy.ones(N+1)
    mvec = numpy.linspace(0,N,N+1)

    for j in xrange(1,N+1):
        y += k * j * B[j] * numpy.cosh(numpy.dot(j * k,eta)) / numpy.cosh(j * k * D) *\
             numpy.cos(numpy.dot(mvec ,float(j) * numpy.pi / float(N)) )
    return y

def v(eta,B,D,k,N):

    y    = numpy.zeros(N+1)
    mvec = numpy.linspace(0,N,N+1) 

    for j in xrange(1,N+1):
        y += k * j * B[j] * numpy.sinh(numpy.dot(j*k,eta)) / numpy.cosh(j * k * D) *\
             numpy.sin(numpy.dot(j*numpy.pi/N,mvec))
    return y


def Sjm_one(k,eta,D,N):
    y    = numpy.zeros((N+1,N))
    mvec = numpy.linspace(0,N,N+1)
    jvec = numpy.linspace(1,N,N)

    y += numpy.multiply(numpy.divide(numpy.sinh(numpy.outer(eta,jvec*k)),numpy.cosh(jvec*k*D)),\
         numpy.cos(numpy.outer(mvec,jvec*numpy.pi/N)))

    return y

def Sjm_two(k,eta,D,N):
    y    = numpy.zeros((N+1,N))
    mvec = numpy.linspace(0,N,N+1)
    jvec = numpy.linspace(1,N,N)

    y += numpy.multiply(numpy.divide(numpy.sinh(numpy.outer(eta,jvec*k)),numpy.cosh(jvec*k*D)),\
         numpy.sin(numpy.outer(mvec,jvec*numpy.pi/N)))

    return y

def Cjm_one(k,eta,D,N):
    y    = numpy.zeros((N+1,N))
    mvec = numpy.linspace(0,N,N+1)
    jvec = numpy.linspace(1,N,N)

    y += numpy.multiply(numpy.divide(numpy.cosh(numpy.outer(eta,jvec*k)),numpy.cosh(jvec*k*D)), \
         numpy.sin(numpy.outer(mvec,jvec*numpy.pi/N)))

    return y

def Cjm_two(k,eta,D,N):
    y    = numpy.zeros((N+1,N))
    mvec = numpy.linspace(0,N,N+1)
    jvec = numpy.linspace(1,N,N)

    y += numpy.multiply(numpy.divide(numpy.cosh(numpy.outer(eta,jvec*k)),numpy.cosh(jvec*k*D)), \
         numpy.cos(numpy.outer(mvec,jvec*numpy.pi/N)))

    return y

#----------------------------------------------------#
# Nonlinear equation to be solved by Newton's method 
# func(z) = 0
# Rienecker & Fenton (1981); eq. (15)
#

def func(eta,B,c,k,Q,R,cE,H,tau,N,D):
    y = numpy.zeros(2*N+6)
    e = numpy.ones(N-1)

    #-- Determine u and v
    mvec = numpy.linspace(0,N,N+1)

    um = u(eta,B,D,k,N)
    vm = v(eta,B,D,k,N)

    y[:N+1]        = numpy.dot(B[0], eta) + Q*numpy.ones(N+1)
    y[N+1:2*(N+1)] = 0.5 * (um**2 + vm**2 ) + eta - R*numpy.ones(N+1)
    y[2*N+2]       = 1./(2*N)*(eta[0]+eta[N]+2.*numpy.dot(eta[1:N],e)) - 1.
    y[2*N+3]       = eta[0] - eta[N] - H
    y[2*N+4]       = k * c * tau - 2. * numpy.pi
    y[2*N+5]       = c - cE + B[0]
    
    for j in xrange(1,N+1):
        y[:N+1] += B[j]*numpy.cos(numpy.dot(mvec, j * numpy.pi / N)) * \
                   numpy.sinh(numpy.dot(j*k,eta)) / \
                   numpy.cosh(j*k*D)

    return y

#-------------------------------------------------#
# Jacobian of func.
#

def jacobian(eta,B,c,k,Q,R,cE,H,tau,N,D):
    y = numpy.zeros((2*N+6,2*N+6))

    #-- Determine u and v
    um = u(eta,B,D,k,N)
    vm = v(eta,B,D,k,N)

    e  = numpy.ones(N+1)
    ee = numpy.ones(N)
    tmp1 = numpy.zeros(N+1)
    tmp2 = numpy.zeros(N+1)

    jvec = numpy.linspace(1,N,N)
    mvec = numpy.linspace(0,N-1,N)

    jBtanh  = jvec * B[1:] * numpy.tanh(jvec*k*D)
    j2Btanh = jvec**2 * B[1:] * numpy.tanh(jvec*k*D)
    j2B     = jvec**2 * B[1:]

    #dfi/detai
    y[:N+1,:N+1]        += numpy.diag(u(eta,B,D,k,N))
    #dfi/dB0
    y[:N+1,N+1]         += numpy.dot(-1,eta)
    #dfi/dBj
    y[:N+1,N+2:2*N+2]   += Sjm_one(k,eta,D,N)
    #dfi/dc = 0
    #dfi/dk 
    y[:N+1,2*N+3]       += numpy.dot(eta * (um-B[0]*e),1./k) - \
                           D * numpy.dot(Sjm_one(k,eta,D,N),jBtanh)


    #dfi/dQ
    y[:N+1,2*N+4]       += e
    #dfi/dR = 0 

    #dfi/detai
    y[N+1:2*(N+1),:N+1]        += e + um * k**2 * numpy.dot(Sjm_one(k,eta,D,N),j2B) + \
                                  vm * k**2 * numpy.dot(Cjm_one(k,eta,D,N),j2B)

    #dfi/dB0
    y[N+1:2*(N+1),N+1]         += numpy.dot(-1,um) 
    #dfi/dBj
    y[N+1:2*N+2,N+2:2*N+2] += numpy.dot(numpy.dot(numpy.diag(um), Cjm_two(k,eta,D,N)), numpy.diag(jvec * k)) + \
                              numpy.dot(numpy.dot(numpy.diag(vm), Sjm_two(k,eta,D,N)), numpy.diag(jvec * k))
    #dfi/dc = 0
    #dfi/dk 
    
    tmp1 = (um - B[0]*e) / k
    tmp2 = vm / k

    tmp1 += numpy.dot(k,eta) * numpy.dot(Sjm_one(k,eta,D,N),j2B) - \
            k * D * numpy.dot(Cjm_two(k,eta,D,N),j2Btanh)
    tmp2 += numpy.dot(k,eta) * numpy.dot(Cjm_one(k,eta,D,N),j2B) - \
            k * D * numpy.dot(Sjm_two(k,eta,D,N),j2Btanh)
    
    y[N+1:2*(N+1),2*N+3]       += um*tmp1 + vm*tmp2
    #dfi/dQ=0
    #dfi/dR 
    y[N+1:2*(N+1),2*N+5]       += numpy.dot(-1,e)
    
    #df_{2N+3}/detaj
    y[2*N+2,0]   += 1 / (2 * float(N))
    y[2*N+2,1:N] += 1 / float(N)
    y[2*N+2,N]   += 1 / (2 * float(N))

    #df_{2N+4}/detaj
    y[2*N+3,0] +=  1.0
    y[2*N+3,N] += -1.0

    #df_{2N+5}/dc
    y[2*N+4,2*N+2] += k * tau
    #df_{2N+5}/dk
    y[2*N+4,2*N+3] += c * tau

    #df_{2N+6}/dc
    y[2*N+5,2*N+2]  +=1.
    #df_{2N+6}/dB0
    y[2*N+5,N+1]    +=1.
    #df_{2N+6}/dQ
    #y[2*N+5,2*N+4]  +=-1

    return y

#-------------------------------------------------#
# The Rienecker & Fenton spectrum
# Parameters to be modified by user 
#

def RieneckerFenton(kx,ky):

    Nx = len(kx)
    Ny = len(ky)

    N = Nx-1   #-- number of unknowns
    D = 1.0  #-- arbitrary reference level


    #-- Initialize quantities PLUS initial guess
    eta  = numpy.ones(N+1)
    B    = numpy.zeros(N+1)
    mvec = numpy.linspace(0,N,N+1)
    
    H   = .4  # wave height
    c   = 2.0
    tau = 10. # wave period
    cE  = 0.0 # current speed in Eulerian frame
    
    k   = 2.*numpy.pi / (tau*c)
    c   = numpy.sqrt(numpy.tanh(k)/k)
    
    per   = 5. 
    eta  += 0.5 * H * numpy.cos(numpy.dot(mvec,per*numpy.pi/N))
    
    B[0] = -c
    B[1] = -0.25*H/c*k
    
    D   = 1.
    R   = 1. + 0.5*c**2
    Q   = c
    
    
    #-- The Newton method will compute non-dimensionalized variables
    lmbda = 2.*numpy.pi / k
    
    x = numpy.linspace(0,N,N+1)*lmbda/(2.*float(N))
    y = x.copy()
    
    z0 = numpy.concatenate((eta,B)) 
    z0 = numpy.append(z0,numpy.array([c,k,Q,R]),axis=1)
    
    #print "The initial guess is z_0 = %s." %z0
    
    zn  = numpy.zeros(2*N+6)
    #ind = numpy.amax(abs(func(z0[:N+1],z0[N+1:2*N+2],z0[2*N+2],z0[2*N+3],z0[2*N+4],z0[2*N+5],cE,H,tau,N,D)))
    
    #while (ind > 10e-2):
    for i in xrange(0,4):
        f0  = func(z0[:N+1],z0[N+1:2*N+2],z0[2*N+2],z0[2*N+3],z0[2*N+4],z0[2*N+5],cE,H,tau,N,D)
        Df0 = jacobian(z0[:N+1],z0[N+1:2*N+2],z0[2*N+2],z0[2*N+3],z0[2*N+4],z0[2*N+5],cE,H,tau,N,D)
        
        delta_z = numpy.linalg.solve(Df0,-f0)
        
        zn = z0 + delta_z

        #ind = numpy.amax(abs(func(zn[:N+1],zn[N+1:2*N+2],zn[2*N+2],zn[2*N+3],zn[2*N+4],zn[2*N+5],cE,H,tau,N,D)))
        z0  = zn
        
        #print "The final result is z_n = %s." %zn

        #-- Use z_n to compute the potential and the initial water height.
        eta = zn[0:N+1]
        B   = zn[N+1:2*N+2]
        c   = zn[2*N+2] 
        k   = zn[2*N+3]
        Q   = zn[2*N+4]
        R   = zn[2*N+5]
        
        #-- Return dimensional variables
        psi = compute_potential(x,eta,B,k,D)
        eta_mean = 0.5 * (eta[0] + eta[N])
        
        x   = x   * eta_mean
        eta = eta * eta_mean
        psi = psi * numpy.sqrt(9.81 * eta_mean**3)
        c   = c * numpy.sqrt(9.81 * eta_mean)
        k   = k / eta_mean
        Q   = Q * numpy.sqrt(9.81 * eta_mean**3)
        R   = R * (9.81 * eta_mean)
        D   = D * eta_mean

        etar = numpy.reshape(numpy.repeat(eta,numpy.size(x),axis=0),\
                             (numpy.size(x),numpy.size(x)),order="F")

        phir = numpy.reshape(numpy.repeat(psi.diagonal(),numpy.size(x),axis=0),\
                             (numpy.size(x),numpy.size(x)),order="F")

       
        return etar, phir


def print_rieneckerfenton():
    #-- Plot the initial solution
    plt.figure()
    plt.suptitle("Rienecker & Fenton Solution", fontsize = 16)
    
    #        plt.subplot(221)
    #        plt.pcolormesh(x,y,psi, cmap=plt.get_cmap("Greens"))
    #        plt.colorbar()
    #        plt.title("Initial Velocity Potential $\psi$")
    #        plt.xlabel("x")
    #        plt.xlim(numpy.amin(x), numpy.amax(x))
    #        plt.ylabel("y=$\eta$(x)")
    #        plt.ylim(numpy.amin(y), numpy.amax(y))
    
    plt.subplot(221)
    plt.plot(x,eta)
    plt.title("Initial Sea Surface Elevation $\eta$")
    plt.xlabel("x")
    plt.xlim(numpy.amin(x), numpy.amax(x))
    
    plt.subplot(222)
    plt.pcolormesh(x,y,etar,cmap=plt.get_cmap("Blues"))
    plt.colorbar()
    plt.title("Initial Sea Surface Elevation $\eta(x,y)$")
    plt.xlabel("x")
    plt.xlim(numpy.amin(x), numpy.amax(x))
    plt.ylabel("y")
    plt.ylim(numpy.amin(y), numpy.amax(y))
    
    
    plt.subplot(223)
    plt.plot(x,psi.diagonal())
    plt.title("Initial Velocity Potential $\psi$")
    plt.xlabel("x")
    plt.xlim(numpy.amin(x), numpy.amax(x))
    plt.ylabel('$\psi$(x,$\eta$(x))')
    
    plt.subplot(224)
    plt.pcolormesh(x,y,phir,cmap=plt.get_cmap("Greens"))
    plt.colorbar()
    plt.title("Initial Velocity Potential $\psi(x,y)$")
    plt.xlabel("x")
    plt.xlim(numpy.amin(x), numpy.amax(x))
    plt.ylabel("y")
    plt.ylim(numpy.amin(y), numpy.amax(y))
    
    
    plt.tight_layout()
    #plt.show()

    return 0
