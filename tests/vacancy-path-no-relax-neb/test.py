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

import os
import numpy as np

# standard ASE structure generation routines
from ase.lattice.cubic import Diamond
from ase.optimize import FIRE, MDMin
from quippy.neb import NEB
from quippy.atoms import Atoms
from quippy.potential import Minim
from ase.calculators.neighborlist import NeighborList
from numpy import dot, array
from math import sqrt

import ase.io, sys

# set of utility routines specific this this model/testing framework
from utilities import relax_atoms, relax_atoms_cell

# the current model
import model

a0 = 5.44 # initial guess at lattice constant, cell will be relaxed below
tol = 1e-3 # maximum force following relaxtion [eV/A]
fmax = 0.01
N = 2 # number of unit cells in each direction
remove_index = 'center' # which atom to remove - special value 'center' looks for atom in centre of cell
n_images = 9

# set up the a
bulk = Diamond(symbol='Si', latticeconstant=a0)

# specify that we will use model.calculator to compute forces, energies and stresses
bulk.set_calculator(model.calculator)

# use one of the routines from utilities module to relax the initial
# unit cell and atomic positions

if hasattr(model,'bulk_reference'):
    bulk = model.bulk_reference.copy()
    bulk.set_calculator(model.calculator)
else:
    bulk.set_calculator(model.calculator)
    print "BOB relaxing bulk"
    bulk = relax_atoms_cell(bulk, tol=tol, traj_file=None)

print 'bulk reference len=%d cell=%r' % (len(bulk), bulk.cell)

bulk *= (N, N, N)
bulk.set_calculator(model.calculator)
print "BOB evaluating bulk"
e_bulk_per_atom = bulk.get_potential_energy()/len(bulk)
print 'e_bulk_per_atom', e_bulk_per_atom

if remove_index == 'center':
    half_cell = np.diag(bulk.cell)/2.
    remove_index = ((bulk.positions - half_cell)**2).sum(axis=1).argmin()

initial = bulk.copy()
orig_pos = initial.get_positions()[remove_index,:]

# save for initial unrelaxed energy
initial_unrelaxed = initial.copy()
orig_pos_unrelaxed = initial_unrelaxed.get_positions()[remove_index,:]

# create perturbed initial
nl = NeighborList([a0*sqrt(3.0)/4*0.6]*len(initial), self_interaction=False, bothways=True)
nl.update(initial)
indices, offsets = nl.get_neighbors(remove_index)
remove_index_f = None
initial.arrays['ind'] = array([i for i in range(len(initial)) ])
offset_factor=0.13
for i, offset in zip(indices, offsets):
   ri = initial.positions[remove_index] - (initial.positions[i] + dot(offset, initial.get_cell()))
   if remove_index_f is None:
      remove_index_f = i
   print "initial offset ", i, offset_factor, ri
   initial.positions[i] += offset_factor*ri
   offset_factor += 0.01

del initial[remove_index]
initial.set_calculator(model.calculator)

print "BOB relaxing initial"
initial = relax_atoms(initial, tol=tol, traj_file=None)
#initial = ase.io.read(os.path.join("initial-relax", "castep.castep"))
initial.set_calculator(model.calculator)

# print 'initial relaxed energy', initial.get_potential_energy()
# e_form = initial.get_potential_energy() - e_bulk_per_atom*initial.get_number_of_atoms()
# print 'relaxed vacancy formation energy', e_form

# make final
final = bulk.copy()

nl = NeighborList([a0*sqrt(3.0)/4*0.6]*len(final), self_interaction=False, bothways=True)
nl.update(final)
indices, offsets = nl.get_neighbors(remove_index_f)
offset_factor=0.13
for i, offset in zip(indices, offsets):
   ri = final.positions[remove_index_f] - (final.positions[i] + dot(offset, final.get_cell()))
   final.positions[i] += offset_factor*ri
   offset_factor += 0.01

del final[remove_index]
for (i,ii) in enumerate(initial.arrays['ind']):
   if ii == remove_index_f:
      new_remove_index_f = i
remove_index_f = new_remove_index_f
final.positions[remove_index_f] = orig_pos
final.set_calculator(model.calculator)

print "BOB relaxing final config"
final = relax_atoms(final, tol=tol, traj_file=None)
#final = ase.io.read(os.path.join("final-relax", "castep.castep"))
final.set_calculator(model.calculator)

# print 'final relaxed energy', final.get_potential_energy()
# e_form = final.get_potential_energy() - e_bulk_per_atom*final.get_number_of_atoms()
# print 'relaxed vacancy formation energy', e_form

#make interpoalted chain
images = [ Atoms(initial) ]
images += [ Atoms(initial.copy()) for i in range(n_images) ]
images += [ Atoms(final) ]
neb = NEB(images, k=0.1)
neb.interpolate()

for img in neb.images:
    img.set_calculator(model.calculator)
    # print "image energy", img.get_potential_energy()

# # make chain
# images = [ Atoms(initial) ]
# images += [ Atoms(initial.copy()) for i in range(n_images) ]
# images += [ Atoms(final) ]
# 
# neb = NEB(images, k=0.1)
# neb.interpolate()
# for img in images:
#    img.set_calculator(model.calculator)
#    ase.io.write(sys.stdout, img, 'extxyz')
# 
# # perturb intermediate images
# for img in images[1:-1]:
#    img.rattle(0.05)
# 
# # optimizer = FIRE(neb)
# optimizer = MDMin(neb, dt=0.05)
# optimizer.run(fmax=fmax)

neb_traj = open('model-'+model.name+'-vacancy-path.neb.xyz', 'w')
ase.io.write(neb_traj, neb.images, format='extxyz')
neb_traj.close()

Es = []
for img in images:
   print "BOB evaluating img"
   Es.append(img.get_potential_energy() - e_bulk_per_atom*len(img))
   ase.io.write(sys.stdout, img, 'extxyz')

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = {'vacancy_path':
                Es }
