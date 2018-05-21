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
from ase.calculators.neighborlist import NeighborList
from numpy import dot, sum, amax
from math import sqrt

from quippy import Potential

# set of utility routines specific this this model/testing framework
from utilities import relax_atoms, relax_atoms_cell

# the current model
import model

a0 = 5.44 # initial guess at lattice constant, cell will be relaxed below
fmax = 1e-3 # maximum force following relaxtion [eV/A]
N = 3 # number of unit cells in each direction
remove_index = 0 # which atom to remove

if not hasattr(model, 'bulk_reference_216'):
    # set up the a
    bulk = Diamond(symbol='Si', latticeconstant=a0)

    # specify that we will use model.calculator to compute forces, energies and stresses
    bulk.set_calculator(model.calculator)

    # use one of the routines from utilities module to relax the initial
    # unit cell and atomic positions
    bulk = relax_atoms_cell(bulk, tol=fmax, traj_file=None)
    bulk *= (N, N, N)
    bulk_energy = bulk.get_potential_energy()
else:
    bulk = model.bulk_reference_216
    bulk_energy = bulk.get_potential_energy()

def vacancy_energy(bulk, remove_index=0):
    Nat = bulk.get_number_of_atoms()
    vac = bulk.copy()
    vac.set_calculator(bulk.get_calculator())

    nl = NeighborList([a0*sqrt(3.0)/4*0.6]*len(bulk), self_interaction=False, bothways=True)
    nl.update(bulk)
    indices, offsets = nl.get_neighbors(remove_index)
    offset_factor=0.13
    for i, offset in zip(indices, offsets):
       ri = vac.positions[remove_index] - (vac.positions[i] + dot(offset, vac.get_cell()))
       vac.positions[i] += offset_factor*ri
       offset_factor += 0.01

    del vac[remove_index] # remove an atom to introduce a vacancy

    # perturb positions
    vac.rattle(0.1)

    ##
    try:
        model.calculator.set(local_gap_error='local_gap_error')
        vac.get_potential_energy()
        unrelaxed_local_gap_error_max = amax(model.calculator.results['local_gap_error'])
        unrelaxed_local_gap_error_sum = sum(model.calculator.results['local_gap_error'])
    except:
        unrelaxed_local_gap_error_max = None
        unrelaxed_local_gap_error_sum = None
    if isinstance(model, Potential):
        model.calculator.set(local_gap_error='')

    # relax atom positions, holding cell fixed
    vac = relax_atoms(vac, tol=fmax, traj_file="model-"+model.name+"-test-vacancy-energy.opt.xyz")

    vac.write("model-"+model.name+"test-vacancy-energy-relaxed.xyz")

    ##
    try:
        model.calculator.set(local_gap_error='local_gap_error')
        vac.get_potential_energy()
        relaxed_local_gap_error_max = amax(model.calculator.results['local_gap_error'])
        relaxed_local_gap_error_sum = sum(model.calculator.results['local_gap_error'])
    except:
        relaxed_local_gap_error_max = None
        relaxed_local_gap_error_sum = None
    if isinstance(model, Potential):
        model.calculator.set(local_gap_error='')
    ##

    # compute vacancy formation energy as difference of bulk and vac energies
    print 'bulk cell energy', bulk_energy
    print 'vacancy cell energy', vac.get_potential_energy()
    e_form = vac.get_potential_energy() - bulk_energy*vac.get_number_of_atoms()/bulk.get_number_of_atoms()
    print 'vacancy formation energy', e_form
    return (e_form, unrelaxed_local_gap_error_max, unrelaxed_local_gap_error_sum, relaxed_local_gap_error_max, relaxed_local_gap_error_sum)

# dictionary of computed properties - this is output of this test, to
#   be compared with other models

(e_form, unrelaxed_local_gap_error_max, unrelaxed_local_gap_error_sum, relaxed_local_gap_error_max, relaxed_local_gap_error_sum) = vacancy_energy(bulk, remove_index=remove_index)
properties = {'vacancy_energy': e_form,
              'vacancy_unrelaxed_local_gap_error_max' : unrelaxed_local_gap_error_max,
              'vacancy_unrelaxed_local_gap_error_sum' : unrelaxed_local_gap_error_sum,
              'vacancy_relaxed_local_gap_error_max' : relaxed_local_gap_error_max,
              'vacancy_relaxed_local_gap_error_sum' : relaxed_local_gap_error_sum }
