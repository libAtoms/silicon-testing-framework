# This script defines a test case which computes one or more physical
# properties with a given model
#
# INPUTS:
#   model.calculator -- an ase.calculator.Calculator instance
#     this script can assume the calculator is checkpointed.
#
# OUTPUTS:
#   properties -- dictionary of key/value pairs corresponding
#     to physical quantities computed by this test

# standard ASE structure generation routines
from ase.lattice.cubic import Diamond
from ase.constraints import FixAtoms
import numpy as np

import ase.io, sys, os

# set of utility routines specific this this model/testing framework
from utilities import relax_atoms, relax_atoms_cell

# the current model
import model

a0 = 5.44 # initial guess at lattice constant, cell will be relaxed below
fmax = 0.001 # maximum force following relaxtion [eV/A]

if not hasattr(model, 'bulk_reference'):
    # set up the a
    bulk = Diamond(symbol='Si', latticeconstant=a0)

    # specify that we will use model.calculator to compute forces, energies and stresses
    bulk.set_calculator(model.calculator)

    # use one of the routines from utilities module to relax the initial
    # unit cell and atomic positions
    bulk = relax_atoms_cell(bulk, tol=fmax, traj_file=None)
else:
    bulk = model.bulk_reference.copy()
    bulk.set_calculator(model.calculator)

a0 = bulk.cell[0,0] # get lattice constant from relaxed bulk
print "got a0 ", a0
bulk = Diamond(symbol="Si", latticeconstant=a0, directions=[[1,-1,0],[1,0,-1],[1,1,1]])
bulk.set_calculator(model.calculator)

# compute surface formation energy as difference of periodic and curface cell
ebulk_per_atom = bulk.get_potential_energy()/bulk.get_number_of_atoms()
print 'bulk per atom energy', ebulk_per_atom
bulk *= (1,1,10)
ebulk = bulk.get_potential_energy()
print 'bulk 1x1x10 cell energy', ebulk
# now calculate the relaxed unreconstructed 111 surface energy
bulk.positions[:,2] += 2.0 # select shuffle plane to be cut
bulk.wrap()
bulk.cell[2,:] += [0.0,0.0,10.0]
np.random.seed(75)
bulkNat = bulk.get_number_of_atoms()
#bulk.positions += (np.random.rand((bulkNat*3))*0.1).reshape([bulkNat,3])
#bulk = relax_atoms(bulk, tol=fmax, traj_file="model-"+model.name+"-surface-energy-111-expanded-bulk-relaxed.opt.xyz")
eexp = bulk.get_potential_energy()
print 'expanded cell energy', eexp
e111 = 0.5*(eexp-ebulk) / np.linalg.norm(np.cross(bulk.cell[0,:],bulk.cell[1,:]))
print '111 unreconstructed surface energy', e111


def surface_energy(n, bulk, ebulk_per_atom, e111):

    tmp = ase.io.read(os.path.dirname(__file__)+'/si111_'+str(n)+'x'+str(n)+'_8layer.xyz', index=':', format='extxyz')
    surf = tmp[0]
    
    # adjust lattice constant
    L = surf.get_cell()
    bulkL = bulk.get_cell()
    scaling_factor = bulkL[0,0]*n/L[0,0]
    surf.set_cell(L*scaling_factor, scale_atoms=True)
    
    fixed_list = np.where(surf.arrays['move_mask'] == 0)[0]
    c = FixAtoms(fixed_list)
    surf.set_constraint(c)

    # relax atom positions of the reconstructed surface, holding cell fixed
    surfNat = surf.get_number_of_atoms()
    #surf.positions += (np.random.rand((surfNat*3))*0.1).reshape([surfNat,3])
    surf = relax_atoms(surf, tol=fmax, traj_file="model-"+model.name+"-surface-energy-111-"+str(n)+"x"+str(n)+"-reconstruction-relaxed.opt.xyz")
    ase.io.write(sys.stdout, surf, format='extxyz')

    esurf  = surf.get_potential_energy()
    print '111 '+str(n)+'x'+str(n)+' reconstructed surface cell energy', esurf
    surf_area = np.linalg.norm(np.cross(surf.cell[0,:],surf.cell[1,:]))
    e_form = (esurf-surfNat*ebulk_per_atom) / surf_area - e111
    print '111 '+str(n)+'x'+str(n)+' reconstructed surface formation energy', e_form
    return e_form

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = {'surface_energy_111_3x3-reconstructed_relaxed':
                surface_energy(3, bulk, ebulk_per_atom, e111) ,
              'surface_energy_111_5x5-reconstructed_relaxed':
                surface_energy(5, bulk, ebulk_per_atom, e111) ,
              'surface_energy_111_7x7-reconstructed_relaxed':
                surface_energy(7, bulk, ebulk_per_atom, e111) ,
              'surface_energy_111_9x9-reconstructed_relaxed':
                surface_energy(9, bulk, ebulk_per_atom, e111)}
