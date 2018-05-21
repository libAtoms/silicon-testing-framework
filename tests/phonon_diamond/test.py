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

import model
import json

from quippy import diamond2
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

from ase import Atoms
from utilities import relax_atoms_cell, phonons
import numpy as np

# set up cell

try:
   with open("../model-{}-test-bulk_diamond-properties.json".format(model.name)) as f:
       j=json.load(f)
   a0 = j["diamond_a0"]
   bulk = Atoms(diamond2(a0,14,14))
except:
   bulk = Atoms(diamond2(5.43,14,14))
   bulk.set_calculator(model.calculator)
   bulk = relax_atoms_cell(bulk, tol=1.0e-4, traj_file=None)
   a0 = bulk.get_cell_lengths_and_angles()[0]*np.sqrt(2.0)
   bulk = Atoms(diamond2(a0,14,14))

supercell = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]

dx = 0.03

mesh = [16, 16, 16]

points = np.array([
    [0,0,0],
    [0,0.5,0.5],
    [1,1,1],
    [0.5,0.5,0.5]
    ]
)

n_points = 50

phonon_properties = phonons(model,bulk,supercell=supercell,dx=dx,mesh=mesh,points=points,n_points=n_points)

properties = {
              "a0":a0,
              "phonon_diamond_frequencies":phonon_properties["frequencies"],
              "phonon_diamond_weights":phonon_properties["weights"],
              "phonon_diamond_band_frequencies":phonon_properties["band_frequencies"],
              "phonon_diamond_band_distances":phonon_properties["band_distances"]
              }

