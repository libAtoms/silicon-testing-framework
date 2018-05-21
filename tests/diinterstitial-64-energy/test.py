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
import numpy as np

import ase.io, sys, os

# set of utility routines specific this this model/testing framework
from utilities import relax_atoms, relax_atoms_cell

# the current model
import model

a0 = 5.44 # initial guess at lattice constant, cell will be relaxed below
fmax = 0.01 # maximum force following relaxtion [eV/A]

if not hasattr(model, 'bulk_reference'):
    print("Did not find bulk_reference, recalculating bulk")
    # set up the a
    bulk = Diamond(symbol='Si', latticeconstant=a0)

    # specify that we will use model.calculator to compute forces, energies and stresses
    bulk.set_calculator(model.calculator)

    # use one of the routines from utilities module to relax the initial
    # unit cell and atomic positions
    bulk = relax_atoms_cell(bulk, tol=fmax, traj_file=None)
else:
    bulk = model.bulk_reference.copy()
    print("Found bulk_reference, lattice constant", bulk.get_cell()[0,0] )
    bulk.set_calculator(model.bulk_reference.calc)

a0 = bulk.cell[0,0] # get lattice constant from relaxed bulk
e0 = bulk.get_potential_energy()/8.0
print "got a0 ", a0
print "got e0 ", e0


def defect_energy(a0, e0, filename):

    
    tmp = ase.io.read(os.path.dirname(__file__)+"/"+filename, index=':', format='extxyz')
    defect = tmp[0]
    
    # adjust lattice constant
    # defect.set_cell([2*a0, 2*a0, 2*a0], scale_atoms=True)
    defect.set_calculator(model.calculator)
    
    # relax atom positions and cell
    # defect.positions += (np.random.rand((66*3))*0.01).reshape([66,3])
    defect = relax_atoms_cell(defect, tol=fmax, traj_file="model-"+model.name+"-diinterstitial-"+filename+"-relaxed.opt.xyz")
    ase.io.write(sys.stdout, defect, format='extxyz')

    edefect  = defect.get_potential_energy()
    print 'defect relaxed cell energy', edefect
    e_form = (edefect-66.0 * e0) 
    print 'defect relaxed formation energy', e_form
    return e_form

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = {'diinterstitial_formation_energy_EXT':
                defect_energy(a0, e0, "Si64_I2_EXT_1S_castep_relaxation_end.xyz") ,
'diinterstitial_formation_energy_TT':
                defect_energy(a0, e0, "Si64_I2_TT_1NNN_castep_relaxation_end.xyz") ,
'diinterstitial_formation_energy_XEX':
                defect_energy(a0, e0, "Si64_I2_XEX_0_castep_relaxation_end.xyz") ,
'diinterstitial_formation_energy_XT':
                defect_energy(a0, e0, "Si64_I2_XT_01L_castep_relaxation_end.xyz") ,
'diinterstitial_formation_energy_W':
                defect_energy(a0, e0, "Si64_I2_XT_02A_castep_relaxation_end.xyz") ,
'diinterstitial_formation_energy_XX3':
                defect_energy(a0, e0, "Si64_I2_XX3_00A_castep_relaxation_end.xyz")}

# 'diinterstitial_formation_energy_XX2':
#                 defect_energy(a0, e0, "Si64_I2_XX2_00A.xyz") ,
