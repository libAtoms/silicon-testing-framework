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
import ase.io, os
from math import sqrt

import lattice_tetragonal

# the current model
import model 

v0 = 20.0

# set up the a
bulk = ase.io.read(os.path.join(os.path.dirname(__file__),"st12.extxyz"))
cell_factor = (len(bulk)*v0/bulk.get_volume())**(1.0/3.0)
bulk.set_cell(bulk.get_cell()*cell_factor, scale_atoms=True)

(E_vs_V) = lattice_tetragonal.do_lattice(bulk, elastic=False)

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = { 'st12_E_vs_V': E_vs_V }
