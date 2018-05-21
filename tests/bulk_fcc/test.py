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
from ase.lattice.cubic import FaceCenteredCubic

import lattice_cubic

# the current model
import model 

a0 = (15.0*4)**(1.0/3.0)# initial guess at lattice constant, cell will be relaxed below

# set up the a
bulk = FaceCenteredCubic(symbol='Si', latticeconstant=a0)

(E_vs_V) = lattice_cubic.do_lattice(bulk, elastic=False)

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = {'fcc_E_vs_V': E_vs_V }
