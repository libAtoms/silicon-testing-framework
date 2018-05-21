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
tol = 3e-3 # maximum force following relaxtion [eV/A]
N = 3 # number of unit cells in each direction

if not hasattr(model, 'bulk_reference_216'):
    # set up the a
    bulk = Diamond(symbol='Si', latticeconstant=a0)

    # specify that we will use model.calculator to compute forces, energies and stresses
    bulk.set_calculator(model.calculator)

    # use one of the routines from utilities module to relax the initial
    # unit cell and atomic positions
    bulk = relax_atoms_cell(bulk, tol=tol, traj_file=None)
    bulk *= (N, N, N)
    bulk_energy = bulk.get_potential_energy()
else:
    bulk = model.bulk_reference_216
    bulk_energy = bulk.get_potential_energy()

def tetrahedral_interstitial_energy(bulk):
    Nat = bulk.get_number_of_atoms()
    int_struct = bulk.copy()
    int_struct.set_calculator(bulk.get_calculator())

    # add an atom to introduce an interstitial
    int_struct.append(Atom('Si', (0.001, 0.002, 5.44/2.0+0.003)))

    ase.io.write(sys.stdout, int_struct, format='extxyz')
    # relax atom positions, holding cell fixed
    int_struct = relax_atoms(int_struct, tol=tol, traj_file="model-"+model.name+"-test-interstitial-tetrahedral.opt.xyz")
    ase.io.write(sys.stdout, int_struct, format='extxyz')

    # compute formation energy as difference of bulk and int energies
    print 'bulk cell energy', bulk_energy
    print 'interstitial cell energy', int_struct.get_potential_energy()
    e_form = int_struct.get_potential_energy() - bulk_energy*((Nat+1.0)/Nat)
    print 'interstitial formation energy', e_form
    return e_form

def dumbbell_interstitial_energy(bulk):
    Nat = bulk.get_number_of_atoms()
    int_struct = bulk.copy()
    int_struct.set_calculator(bulk.get_calculator())

    # add an atom to introduce an interstitial
    int_struct.append(Atom('Si', (-0.5, 0.5, 5.44/2.0+1.0)))
    p = int_struct.get_positions()
    p[149,0] -= 1.0
    p[149,1] += 1.0
    p[149,2] -= 0.5
    int_struct.set_positions(p)

    ase.io.write(sys.stdout, int_struct, format='extxyz')
    # relax atom positions, holding cell fixed
    int_struct = relax_atoms(int_struct, tol=tol, traj_file="model-"+model.name+"-test-interstitial-dumbbell.opt.xyz")
    ase.io.write(sys.stdout, int_struct, format='extxyz')

    # compute formation energy as difference of bulk and int energies
    print 'bulk cell energy', bulk_energy
    print 'interstitial cell energy', int_struct.get_potential_energy()
    e_form = int_struct.get_potential_energy() - bulk_energy*((Nat+1.0)/Nat)
    print 'interstitial formation energy', e_form
    return e_form


def hexagonal_interstitial_energy(bulk):
    Nat = bulk.get_number_of_atoms()
    int_struct = bulk.copy()
    int_struct.set_calculator(bulk.get_calculator())
    # add an atom to introduce an interstitial
    pos = int_struct.get_positions()
    int_pos = (pos[83,:] + pos[101,:] + pos[105,:] + pos[106,:] + pos[108,:] + pos[176,:])/6.0

    print "int_pos ", int_pos
    int_struct.append(Atom('Si', int_pos))

    rot_mat = np.zeros( (3,3) )
    rot_mat[0,:] = (1.0/np.sqrt(2.0), -1.0/np.sqrt(2.0), 0.0)
    rot_mat[1,:] = (1.0/np.sqrt(6.0), 1.0/np.sqrt(6.0), -2.0/np.sqrt(6.0))
    rot_mat[2,:] = (1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0))
    cell = np.dot(int_struct.get_cell(), rot_mat.T)
    int_struct.set_cell(cell, scale_atoms=True)

    int_struct.arrays['move_mask_3'] = np.ones( (len(int_struct), 3), dtype=int )
    int_struct.arrays['move_mask_3'][len(int_struct)-1,2] = 0

    ase.io.write(sys.stdout, int_struct, format='extxyz')
    # relax atom positions, holding cell fixed
    int_struct = relax_atoms(int_struct, tol=tol, traj_file="model-"+model.name+"-test-interstitial-hexagonal.opt.xyz", method='cg_n')
    ase.io.write(sys.stdout, int_struct, format='extxyz')

    # compute formation energy as difference of bulk and int energies
    print 'bulk cell energy', bulk_energy
    print 'interstitial cell energy', int_struct.get_potential_energy()
    e_form = int_struct.get_potential_energy() - bulk_energy*((Nat+1.0)/Nat)
    print 'interstitial formation energy', e_form
    return e_form

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = {'tetrahedral_interstitial_energy':
                tetrahedral_interstitial_energy(bulk),
              'hexagonal_interstitial_energy':
                hexagonal_interstitial_energy(bulk),
              'dumbbell_interstitial_energy':
                dumbbell_interstitial_energy(bulk)}
