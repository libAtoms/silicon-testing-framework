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
from ase import Atom
import numpy as np

import ase.io, sys

# set of utility routines specific this this model/testing framework
from utilities import relax_atoms, relax_atoms_cell

# the current model
import model

a0 = 5.44 # initial guess at lattice constant, cell will be relaxed below
tol = 1e-3 # maximum force following relaxtion [eV/A]
N = 2 # number of unit cells in each direction

if not hasattr(model, 'bulk_reference_216'):
    # set up the a
    bulk = Diamond(symbol='Si', latticeconstant=a0)

    # specify that we will use model.calculator to compute forces, energies and stresses
    bulk.set_calculator(model.calculator)

    # use one of the routines from utilities module to relax the initial
    # unit cell and atomic positions
    print "JOE test relaxing bulk"
    bulk = relax_atoms_cell(bulk, tol=tol, traj_file=None,method='lbfgs')
    print "JOE test done relaxing bulk"
    bulk *= (N, N, N)
    print "JOE test reading bulk E"
    bulk_energy = bulk.get_potential_energy()
    print "JOE test done reading bulk E"
else:
    bulk = model.bulk_reference_216
    bulk_energy = bulk.get_potential_energy()

def fourfold_defect_energy(bulk):
    Nat = bulk.get_number_of_atoms()
    defect = bulk.copy()
    defect.set_calculator(bulk.get_calculator())

    # modify atom positions so that relaxation creates a fourfold defect
    p = defect.get_positions()
    # p[0,0] += 1.2
    # p[0,1] += 1.2
    # p[1,0] -= 1.2
    # p[1,1] -= 1.2
    # from CASTEP relaxed structure
    p [0,:] += (0.81, 0.81, -0.43)
    p [1,:] -= (0.81, 0.81, -0.43)

    defect.set_positions(p)
    
    ase.io.write(sys.stdout, defect, format='extxyz')
    # relax atom positions, holding cell fixed
    print "JOE test relaxing defect"
    defect = relax_atoms(defect, tol=tol, traj_file="model-"+model.name+"-test-fourfold-defect.opt.xyz", method="lbfgs")
    print "JOE test done relaxing defect"
    ase.io.write(sys.stdout, defect, format='extxyz')

    # compute formation energy as difference of bulk and int energies
    print 'bulk energy', bulk_energy
    print 'defect energy', defect.get_potential_energy()
    e_form = defect.get_potential_energy() - bulk_energy
    print 'defect formation energy', e_form
    return e_form


# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = {'fourfold_defect_energy':
                fourfold_defect_energy(bulk)}
