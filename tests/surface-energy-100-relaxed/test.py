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

import ase.io, sys

import numpy as np
# from pylab import *

# set of utility routines specific this this model/testing framework
from utilities import relax_atoms, relax_atoms_cell

# the current model
import model

a0 = 5.44 # initial guess at lattice constant, cell will be relaxed below
fmax = 0.001 # maximum force following relaxtion [eV/A]

# if not hasattr(model, 'bulk_reference'):

# set up the a
bulk = Diamond(symbol='Si', latticeconstant=a0, directions=[[1,1,0],[1,-1,0],[0,0,1]])

# specify that we will use model.calculator to compute forces, energies and stresses
bulk.set_calculator(model.calculator)

# use one of the routines from utilities module to relax the initial
# unit cell and atomic positions
bulk = relax_atoms_cell(bulk, tol=fmax, traj_file=None)

# else:
    # bulk = model.bulk_reference.copy()
    # bulk.set_calculator(model.calculator)    

# set up supercell
big_bulk = bulk*(1,1,10)
big_bulk.set_calculator(model.calculator)

def surface_energy(bulk):
    Nat = bulk.get_number_of_atoms()

    # relax atom positions, holding cell fixed
    # vac = relax_atoms(vac, fmax=fmax)

    # compute surface formation energy as half the difference of bulk and expanded cell
    print "len(bulk)",len(bulk)
    ebulk = bulk.get_potential_energy()
    print 'bulk cell energy', ebulk

    bulk.cell[2,:] -= [0.0,0.0,10.0]

    ase.io.write(sys.stdout, bulk, format='extxyz')

    # relax expanded cell
    bulk = relax_atoms(bulk, tol=fmax, method='lbfgs', traj_file="model-"+model.name+"-surface-energy-100-relaxed.opt.xyz", use_armijo=False)
    eexp  = bulk.get_potential_energy()

    ase.io.write(sys.stdout, bulk, format='extxyz')

    print 'expanded cell energy', eexp
    e_form = 0.5*(eexp - ebulk) / np.linalg.norm(np.cross(bulk.cell[0,:],bulk.cell[1,:]))

    print 'relaxed 100 surface formation energy', e_form
    return e_form

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
E = surface_energy(big_bulk)
properties = {'surface_energy_100_relaxed': E }
