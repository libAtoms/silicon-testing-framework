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
from ase import Atoms

import ase.io, sys

import lattice_tetragonal

# the current model
import model 

a0 = (20.0*2)**(1.0/3.0) # initial guess at lattice constant, cell will be relaxed below

# set up the a
bulk = Atoms([14] * 2, positions=[(0.0, -0.25, -0.069), (0.0, 0.25, 0.069)], 
   cell=[ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.5, 0.5, 0.276]], pbc=(1,1,1))

cell = bulk.get_cell()
cell *= (20.0*2/bulk.get_volume())**(1.0/3.0)
bulk.set_cell(cell, scale_atoms=True)

print "unrelaxed bulk"
ase.io.write(sys.stdout, bulk, format='extxyz')

(E_vs_V) = lattice_tetragonal.do_lattice(bulk, use_precon=False, elastic=False)

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = {'beta-Sn_E_vs_V': E_vs_V }
