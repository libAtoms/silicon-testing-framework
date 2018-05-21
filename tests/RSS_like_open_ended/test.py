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
import numpy as np
import sys
import ase, ase.io
from ase.constraints import UnitCellFilter
from utilities import relax_atoms_cell

# the current model
import model 
import os

v0 = 20.0 # initial test of volume
n_at = 8 # number of atoms in cell
N_c = 300 # number of configs to test
max_strain = 0.2 # max strain away from cube
hard_sphere_dist = 1.7 # exclude distance
tol = 0.01
rng_seed = 5

cluster_e_res = 0.01
cluster_v_res = 0.1

def cluster(configs, e_res, v_res):
    unique_structures=[]
    for at in configs:
        e = at.info['relaxed_energy_per_atom']
        v = at.get_volume()/len(at)
        # print "check ", e, v
        new = True
        for (e_i, v_i, at_i) in unique_structures:
            # print "compare to ", e_i, v_i
            if np.abs(e_i-e) < e_res and np.abs(v_i-v) < v_res:
                # print "dup", e_i, e, np.abs(e_i-e), v_i, v, np.abs(v_i-v)
                at_i.info['multiplicity'] += 1
                new=False
                break
        if new:
            # print "new"
            at.info['energy'] = e
            at.info['volume'] = v
            at.info['multiplicity'] = 1
            at_copy = at.copy()
            at_copy.set_calculator(at.get_calculator())
            unique_structures.append((e, v, at_copy))
    return unique_structures

def too_close(at, min_r):
    n = len(at)
    for i in range(n):
        for j in range(i+1,n):
            if at.get_distance(i, j, mic=True) < min_r:
                return True
    return False

np.random.seed(rng_seed)

minima = []
for config_i in range(N_c):
    a0 = (v0*n_at)**(1.0/3.0)
    at = ase.Atoms(numbers = [14] * n_at, cell=[a0,a0,a0], pbc=[1,1,1])

    # strain cell
    F = 2.0*max_strain*(np.random.rand(3,3)-0.5)
    F_diag = (F[0,0]+F[1,1]+F[2,2])/3.0
    F[0,0] += 1.0 - F_diag
    F[1,1] += 1.0 - F_diag
    F[2,2] += 1.0 - F_diag
    F = 0.5*(F+F.T)
    at.set_cell(np.dot(at.get_cell(),F))

    at.set_scaled_positions(np.random.rand(n_at,3))
    while too_close(at, hard_sphere_dist):
        at.set_scaled_positions(np.random.rand(n_at,3))
    print "initial config"
    ase.io.write(sys.stdout, at, format="extxyz")

    at.set_calculator(model.calculator)

    if os.path.basename(os.path.dirname(model.__file__)) == 'CASTEP_file':
        model.geom_method='tpsd'
        model.geom_force_tol=tol
        model.geom_stress_tol=tol
        model.geom_energy_tol=1.0e8
        model.geom_disp_tol=1.0e8
        print "BOB relax atoms"
        at = relax_atoms_cell(at, tol=tol)
        i_minim = -1
        converged = True
    else:
        config_minim = UnitCellFilter(at)
        x = config_minim.get_positions()
        converged = False
        for i_minim in range(500):
            config_minim.set_positions(x)
            grad_f = - config_minim.get_forces()
            E = config_minim.get_potential_energy()
            try:
                pred_err = predictive_error(config_minim.atoms)
            except:
                pred_err = None
            if isinstance(config_minim, UnitCellFilter):
                print config_i,"SD2: {} {} {} {}".format(i_minim, E, np.max(np.abs(grad_f[:n_at])), np.max(np.abs(grad_f[n_at:])))
            else:
                print config_i,"SD2: {} {} {}".format(i_minim, E, np.max(np.abs(grad_f[:n_at])))
            if np.max(np.abs(grad_f[:n_at])) < tol and np.max(np.abs(grad_f[n_at:])) < tol:
                converged = True
                break
            if i_minim == 0:
                alpha = 1.0e-6
            else:
                alpha = np.sum((x-x_old)*(grad_f-grad_f_old)) / np.sum((grad_f-grad_f_old)**2)
            x_old = x.copy()
            grad_f_old = grad_f.copy()
            x -= np.abs(alpha)*grad_f
    if converged:
        print "CONVERGED"
        minima.append(at.copy())
        minima[-1].info['i_minim'] = i_minim
        print "BOB get_potential_energy"
        minima[-1].info['relaxed_energy_per_atom'] = at.get_potential_energy()/len(at)
    else:
        print "FAILED TO CONVERGE"

    print "relaxed config", minima[-1].info['relaxed_energy_per_atom'], minima[-1].get_volume()/len(minima[-1])
    ase.io.write(sys.stdout, minima[-1], format="extxyz")

unique_minima = cluster(minima, e_res = cluster_e_res, v_res = cluster_v_res)
properties = { 'E' : [ x[0] for x in unique_minima ], 'V' : [ x[1] for x in unique_minima ] , 'index' : [ x[2].info['i_minim'] for x in unique_minima ] }
for (i, x) in enumerate(unique_minima):
    ase.io.write("minimum-{}.xyz".format(i), x[2], format="extxyz")
