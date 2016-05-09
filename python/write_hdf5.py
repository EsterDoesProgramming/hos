import h5py
import os
import os.path as op

def write_hdf5(SimFolder, TestName, Nx, Ny, Lx, Ly, g, h, T, dtsave, saveflg, runsubid, Tramp, etar, phir, bathy):

    filename = SimFolder + 'initpars_' + str(TestName) + '.h5'
    
    if op.exists(filename): 
        os.remove(filename)
    
    with h5py.File(filename,'w') as file:

        file.create_dataset('Nx', data = Nx)
        file.create_dataset('Ny', data = Ny)
        file.create_dataset('Lx', data = Lx)
        file.create_dataset('Ly', data = Ly)
        file.create_dataset('g', data = g)
        file.create_dataset('KP', data = 10000)
        file.create_dataset('h', data = h)
        file.create_dataset('T', data = T)
        file.create_dataset('dtsave', data = dtsave)
        file.create_dataset('saveflg', data = saveflg)
        file.create_dataset('runsubid', data = runsubid)
        file.create_dataset('rampflg', data = 1)
        file.create_dataset('Tramp', data = Tramp)
        
    filename = SimFolder + 'initdata_' + str(TestName) + '.h5'
    
    if op.exists(filename): 
        os.remove(filename)
    
    with h5py.File(filename,'w') as file:
    
        file.create_dataset('eta0', data = etar)
        file.create_dataset('phi0', data = phir)
        file.create_dataset('bat0', data = bathy)
        
    print "Success!"
    print "Output files written."
