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

import lattice_cubic

# the current model
import model

a0 = (20.0*8)**(1.0/3.0) # initial guess at lattice constant, cell will be relaxed below

# set up the a
bulk = Diamond(symbol='Si', latticeconstant=a0)

(c11, c12, c44, E_vs_V) = lattice_cubic.do_lattice(bulk, elastic=True)

a0 = bulk.cell[0,0] # save lattice constant after relaxation

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = {'diamond_a0': a0, 'diamond_c11': c11, 'diamond_c12': c12, 'diamond_c44': c44, 'diamond_bulk_modulus': (c11+2.0*c12)/3.0, 'diamond_E_vs_V': E_vs_V }
