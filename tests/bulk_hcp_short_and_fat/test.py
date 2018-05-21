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
from ase.lattice.hexagonal import HexagonalClosedPacked
from math import sqrt

import lattice_tetragonal

# the current model
import model 

c_over_a=1.0
a0 = (16.0*2*2/sqrt(3.0)/c_over_a)**(1.0/3.0)# initial guess at lattice constant, cell will be relaxed below

# set up the a
bulk = HexagonalClosedPacked(symbol='Si', latticeconstant=(a0,a0*c_over_a))

(E_vs_V) = lattice_tetragonal.do_lattice(bulk, elastic=False)

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = {'hcp_short_and_fat_E_vs_V': E_vs_V }
