import numpy

##
# @file bathymetry.py
# @author N. Beisiegel
# @date March 31, 2016
# @brief Provides a routine create_bathymetry that computes an array of length (len(kx),len(ky))
#        which contains the bathymetry.

#-------------------------------------------------#
# Options for different bathymetries.
#

##
# Initializes a zero bathymetry
def bathy_zero(kx,ky):
    return numpy.zeros(len(kx), len(ky))

##
# Initializes a constant bathymetry
def bathy_constant(kx,ky):
    const = 34.0
    return numpy.multiply(const, numpy.ones(len(kx), len(ky)))

##
# Initializes a (in ky direction) linearly sloping bathymetry
def bathy_linear(kx,ky):
    base  = 10.
    slope = 5.
    return numpy.ones((len(kx),1))*(slope*(ky+ky[-1]))+base

##
# Creates an array of bathymetry data corresponding to the input argument Type.
# Zero is default.
#
# @param[in] kx    float array of length Nx; contains wave numbers  
# @param[in] ky    float array of length Ny; contains wave numbers
# @param[in, optional] Type
# @param[out] bathy float array of length Nx times Ny; contains bathymetry information
def create_bathymetry(kx,ky,**kwargs):

    Nx = len(kx)
    Ny = len(ky)

    bathy_flag = kwargs.get('Type',"zero")

    switch = {'zero': bathy_zero, 
              'constant': bathy_constant,
              'linear': bathy_linear}

    if bathy_flag in switch:
        bathy = switch[bathy_flag](kx,ky)
    else: print "Bathymetry type not supported. Please check create_bathymetry.py for options."
    
    return bathy

