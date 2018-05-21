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

import lattice_cubic

# the current model
import model

v0 = 17.0

# set up the a
bulk = ase.io.read(os.path.join(os.path.dirname(__file__),"bc8.extxyz"))
cell_factor = (len(bulk)*v0/bulk.get_volume())**(1.0/3.0)
bulk.set_cell(bulk.get_cell()*cell_factor, scale_atoms=True)

(E_vs_V) = lattice_cubic.do_lattice(bulk, elastic=False)

a0 = bulk.cell[0,0] # save lattice constant after relaxation

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = {'bc8_a0': a0, 'bc8_E_vs_V': E_vs_V }
