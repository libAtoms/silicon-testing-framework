
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
from ase.units import GPa
import sys, ase.io

# set of utility routines specific this this model/testing framework
from utilities import relax_atoms_cell, evaluate
import matscipy.elasticity
from ase.optimize.precon.precon import Exp
from ase.optimize.precon.lbfgs import PreconLBFGS
import model

def do_lattice(bulk, elastic=True):
   tol = 1e-3 # max force tol for relaxation
   n_E_vs_V_steps=10

   # use one of the routines from utilities module to relax the initial
   # unit cell and atomic positions
   bulk = relax_atoms_cell(bulk, tol=tol, traj_file=None, symmetrize=True)

   print "relaxed bulk"
   ase.io.write(sys.stdout, bulk, format='extxyz')

   if elastic:
       # reset calculator to non-symmetrized one (not optimal, but would otherwise need to have optimizer used by fit_elastic_constants to reset symmetry for each relaxation):w
       bulk.set_calculator(model.calculator)
       precon = Exp(3.0)
       opt = lambda atoms, **kwargs: PreconLBFGS(atoms, precon=precon, **kwargs)
       elastic_consts = matscipy.elasticity.fit_elastic_constants(bulk, symmetry='cubic', optimizer=opt)
       c11 = elastic_consts[0][0,0]/GPa
       c12 = elastic_consts[0][0,1]/GPa
       c44 = elastic_consts[0][3,3]/GPa

   V0 = bulk.get_volume()
   dV = bulk.get_volume()*0.025
   E_vs_V=[]

   f = open("relaxed_E_vs_V_configs.xyz", "w")

   for i in range(0, -5-1, -1):
      scaled_bulk = bulk.copy()
      scaled_bulk.set_cell(scaled_bulk.get_cell()*((V0+i*dV)/V0)**(1.0/3.0), scale_atoms=True)
      scaled_bulk = relax_atoms_cell(scaled_bulk, tol=tol, traj_file=None, constant_volume=True, method='cg_n', symmetrize=True)
      # evaluate(scaled_bulk)
      ase.io.write(f, scaled_bulk, format='extxyz')
      E_vs_V.insert(0,  (scaled_bulk.get_volume()/len(scaled_bulk), scaled_bulk.get_potential_energy()/len(bulk)) )

   for i in range(1, 6+1):
      scaled_bulk = bulk.copy()
      scaled_bulk.set_cell(scaled_bulk.get_cell()*((V0+i*dV)/V0)**(1.0/3.0), scale_atoms=True)
      scaled_bulk = relax_atoms_cell(scaled_bulk, tol=tol, traj_file=None, constant_volume=True, method='cg_n', symmetrize=True)
      # evaluate(scaled_bulk)
      ase.io.write(f, scaled_bulk, format='extxyz')
      E_vs_V.append( (scaled_bulk.get_volume()/len(scaled_bulk), scaled_bulk.get_potential_energy()/len(bulk)) )

   for (V, E) in E_vs_V:
     print "EV_final ", V, E

   if elastic:
       return (c11, c12, c44, E_vs_V)
   else:
       return (E_vs_V)
