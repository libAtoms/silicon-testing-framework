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

import ase.io, os

# set of utility routines specific this this model/testing framework
from utilities import relax_atoms_cell

# the current model
import model

fmax = 0.001 # maximum force following relaxtion [eV/A]

if not hasattr(model, 'bulk_reference'):
    # set up the a
    bulk = Diamond(symbol='Si', latticeconstant=5.43)

    # specify that we will use model.calculator to compute forces, energies and stresses
    bulk.set_calculator(model.calculator)

    # use one of the routines from utilities module to relax the initial
    # unit cell and atomic positions
    bulk = relax_atoms_cell(bulk, tol=fmax, traj_file=None)
else:
    bulk = model.bulk_reference.copy()
    bulk.set_calculator(model.calculator)

print "calling bulk.get_potential_energy()"
e0_per_atom = bulk.get_potential_energy()/len(bulk)
print "got e0_per_atom ", e0_per_atom

def traj_energy(e0_per_atom, traj):
    ats = ase.io.read(os.path.dirname(__file__)+"/"+traj, index=':', format='extxyz')
    orig_Es = []
    model_Es = []
    for at in ats:
        # save reference energy
        orig_Es.append(at.get_potential_energy()/len(at))

        # calc new energy
        at.set_calculator(model.calculator)
        model_Es.append(at.get_potential_energy()/len(at))

    return (orig_Es, model_Es)

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
(big_change_ref_E, big_change_model_E) = traj_energy(e0_per_atom, "gap_6_rss_like_open_ended_minima_rerelax_sd2_00099_rerun.traj.xyz")
(small_change_ref_E, small_change_model_E) = traj_energy(e0_per_atom, "gap_6_rss_like_open_ended_minima_rerelax_sd2_00084_rerun.traj.xyz")
properties = {'big_change_ref_E' : big_change_ref_E, 'big_change_model_E' : big_change_model_E,
              'small_change_ref_E' : small_change_ref_E, 'small_change_model_E' : small_change_model_E } 
