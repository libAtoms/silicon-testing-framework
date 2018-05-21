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
import ase, numpy as np
import ase.io, sys, quippy

from ase.constraints import FixCartesian, UnitCellFilter

# set of utility routines specific this this model/testing framework
from utilities import relax_atoms, relax_atoms_cell, evaluate

# the current model
import model

debug = True

a0 = 5.44 # initial guess at lattice constant, cell will be relaxed below
tol = 1e-3 # maximum force following relaxtion [eV/A]
N = 3 # number of unit cells in each direction

do_paths=True
do_surfs=False

n_steps=7
n_steps_path_1=7
n_steps_path_2=5

if hasattr(model, 'bulk_reference'):
    a0 = model.bulk_reference.cell[0, 0]

# set up the bulk cell
bulk = ase.Atoms(['Si']*6,
   positions=[(0,     0,    0),   (a0/2, a0/2, 0),        (a0, a0/2, a0/2),
              (a0/4, a0/4, a0/4), (a0*3/4, a0*3/4, a0/4), (a0*5/4, a0*3/4, a0*3/4)],
   cell=[[a0/2, -a0/2, 0.0], [0.0, a0/2, -a0/2], [a0, a0, a0]],
   pbc=(1,1,1))

# specify that we will use model.calculator to compute forces, energies and stresses
bulk.set_calculator(model.calculator)

# use one of the routines from utilities module to relax the initial
# unit cell and atomic positions
if not hasattr(model, 'bulk_reference'):
    bulk = relax_atoms_cell(bulk, tol=tol, traj_file=None)

bulk *= (1, 1, N)
# rotate to put a1-a2 plane normal to z

bulk_energy = bulk.get_potential_energy()

orig_cell = bulk.get_cell()
new_cell = np.zeros( (3,3) )
new_cell[0,:] = (1.0, 0.0, 0.0)
ang = np.arccos(np.dot(orig_cell[0,:],orig_cell[1,:])/
                (np.linalg.norm(orig_cell[0,:])*np.linalg.norm(orig_cell[1,:])))
new_cell[1,:] = (np.cos(ang), np.sin(ang), 0.0)
new_cell[2,:] = (0.0, 0.0, 1.0)
new_cell[0,:] *= np.linalg.norm(orig_cell[0,:])
new_cell[1,:] *= np.linalg.norm(orig_cell[1,:])
new_cell[2,:] *= np.linalg.norm(orig_cell[2,:])
bulk.set_cell(new_cell, scale_atoms=True)

def sf_calc(bulk, n_steps, z_offset, p0=None, p1=None, relax_from_configs=None):
    global bulk_energy

    p0 = np.array(p0)
    p1 = np.array(p1)
    E0 = bulk_energy
    surf_unrelaxed=np.zeros( (n_steps, n_steps) , dtype=(float,3))
    surf_relaxed=np.zeros( (n_steps, n_steps) , dtype=(float,3))
    path_unrelaxed=[]
    path_relaxed=[]

    sf = bulk.copy()
    pos = sf.get_positions()
    pos[:,2] += z_offset
    sf.set_positions(pos)

    pos = sf.get_scaled_positions()
    pos -= np.floor(pos)
    sf.set_scaled_positions(pos)

    if p0 is not None and p1 is not None: # paths

        # unrelaxed
        for i in range(n_steps):
            print "start path unrelaxed step", i
            sys.stdout.flush()
            p = p0 + (p1-p0)*i/float(n_steps-1)
            sf_cur = sf.copy()

            sf_cur.arrays['move_mask_3'] = np.zeros( (len(sf_cur),3), dtype=int)
            sf_cur.arrays['move_mask_3'][:,2] = 1

            sf_cur.set_calculator(model.calculator)
            cell = sf_cur.get_cell()
            cell[2,:] += cell[0,:]*p[0] + cell[1,:]*p[1]
            sf_cur.set_cell(cell)

            E = (sf_cur.get_potential_energy()-E0)/np.linalg.norm(np.cross(cell[0,:],cell[1,:]))
            path_unrelaxed.append((i/float(n_steps-1), p[0], p[1], E))

        # relaxed
        sf_cur = sf.copy()
        sf_cur.arrays['move_mask_3'] = np.zeros( (len(sf_cur),3), dtype=int)
        sf_cur.arrays['move_mask_3'][:,2] = 1
        sf_cur.set_calculator(model.calculator)
        dp = (p1-p0)/float(n_steps-1)
        p = p0.copy()
        cell = sf_cur.get_cell()
        cell[2,:] += cell[0,:]*p[0] + cell[1,:]*p[1]
        sf_cur.set_cell(cell)
        relax_from_configs_i = 0
        for i in range(n_steps):
            print "start path relaxed step", i
            sys.stdout.flush()

            if model.name == 'CASTEP_file':
                sf_cur.calc.fix_all_cell = False
                sf_cur.set_constraint([FixCartesian(j, [1, 1, 0]) for j in range(len(sf_cur))])
                sf_cur.calc.cell_constraints = ["0 0 1", "0 0 0"]
                sf_cur.calc._old_atoms = None
            if debug:
                ase.io.write(file_unrelaxed,sf_cur,"extxyz")
                file_unrelaxed.flush()

            try:
                sf_cur.set_cell(relax_from_configs[relax_from_configs_i].get_cell(), False)
                sf_cur.set_positions(relax_from_configs[relax_from_configs_i].get_positions())
                print "got sf_cur from 'relax_from_configs'"
                sys.stdout.flush()
            except:
                pass
            sf_cur_rel = relax_atoms_cell(sf_cur, tol=tol, method='cg_n',
                           max_steps=500, traj_file=None, mask=[False, False, True, False, False, False])
            relax_from_configs_i += 1
            if debug:
                ase.io.write(file_relaxed,sf_cur,"extxyz")
                file_relaxed.flush()
            print "done relaxation"
            sys.stdout.flush()

            E = (sf_cur.get_potential_energy()-E0)/np.linalg.norm(np.cross(cell[0,:],cell[1,:]))
            path_relaxed.append((i/float(n_steps-1), p[0], p[1], E))

            cell = sf_cur.get_cell()
            cell[2,:] += cell[0,:]*dp[0] + cell[1,:]*dp[1]
            sf_cur.set_cell(cell)
            p += dp

        return (path_unrelaxed, path_relaxed)

    else: # surfaces

        # unrelaxed
        j_doffset = sf.get_cell()[1,:]*(1.0/float(n_steps-1))
        for i in range(n_steps):
            x_i = float(i)/float(n_steps-1)
            sf_cur = sf.copy()

            sf_cur.arrays['move_mask_3'] = np.zeros( (len(sf_cur),3), dtype=int)
            sf_cur.arrays['move_mask_3'][:,2] = 1

            sf_cur.set_calculator(model.calculator)
            cell = sf_cur.get_cell()
            cell[2,:] += cell[0,:]*x_i
            sf_cur.set_cell(cell)
            for j in range(n_steps-i):
                x_j = float(j)/float(n_steps-1)
                print (i, j, x_i, x_j)
                if debug:
                    ase.io.write(file_unrelaxed,sf_cur,"extxyz")
                    file_unrelaxed.flush()
                E = (sf_cur.get_potential_energy()-E0)/np.linalg.norm(np.cross(cell[0,:],cell[1,:]))
                surf_unrelaxed[i, j] = (x_i, x_j, E)
                if i+j != n_steps-1:
                    surf_unrelaxed[n_steps-j-1, n_steps-i-1] = (1.0-x_j, 1.0-x_i, E)
                cell = sf_cur.get_cell()
                cell[2,:] += j_doffset
                sf_cur.set_cell(cell, scale_atoms=False)

        # relaxed
        relax_from_configs_i = 0
        j_doffset = sf.get_cell()[1,:]*(1.0/float(n_steps-1))
        for i in range(n_steps):
            x_i = float(i)/float(n_steps-1)
            sf_cur = sf.copy()

            # set constraints
            # cs = []
            # for i_at in range(len(bulk)):
                # cs.append(FixCartesian(i_at, mask=[0,0,1]))
            # sf_cur.set_constraint(cs)
            sf_cur.arrays['move_mask_3'] = np.zeros( (len(sf_cur),3), dtype=int)
            sf_cur.arrays['move_mask_3'][:,2] = 1

            sf_cur.set_calculator(model.calculator)
            cell = sf_cur.get_cell()
            cell[2,:] += cell[0,:]*x_i
            sf_cur.set_cell(cell)
            for j in range(n_steps-i):
                x_j = float(j)/float(n_steps-1)
                print (i, j, x_i, x_j)
                try:
                    sf_cur.set_cell(relax_from_configs[relax_from_configs_i].get_cell(), False)
                    sf_cur.set_positions(relax_from_configs[relax_from_configs_i].get_positions())
                except:
                    pass
                sf_cur_rel = relax_atoms_cell(sf_cur, tol=tol, method='cg_n',
                               max_steps=500, traj_file=None, mask=[False, False, True, False, False, False])
                relax_from_configs_i += 1
                if debug:
                    ase.io.write(file_relaxed,sf_cur_rel,"extxyz")
                    file_relaxed.flush()
                E = (sf_cur_rel.get_potential_energy()-E0)/np.linalg.norm(np.cross(cell[0,:],cell[1,:]))
                surf_relaxed[i,j] = (x_i, x_j, E)
                if i+j != n_steps-1:
                    surf_relaxed[n_steps-j-1,n_steps-i-1] = (1.0-x_j, 1.0-x_i, E)
                cell = sf_cur.get_cell()
                cell[2,:] += j_doffset
                sf_cur.set_cell(cell, scale_atoms=False)

        return (surf_unrelaxed, surf_relaxed)

SFE_relaxed = None
SFE_unrelaxed = None

glide_path_unrelaxed=None
glide_path_relaxed=None
shuffle_path_unrelaxed=None
shuffle_path_relaxed=None
glide_surf_unrelaxed=None
glide_surf_relaxed=None
shuffle_surf_unrelaxed=None
shuffle_surf_relaxed=None

if do_paths:
    file_unrelaxed = open("glide_relaxation_start.xyz", "w")
    file_relaxed = open("glide_relaxed.xyz", "w")
    try:
        relax_from_configs = ase.io.read("../model-{}-test-stacking-fault-surface/glide_relaxed.xyz".format("GAP-6"),":")
    except:
        relax_from_configs = None
    (glide_path_unrelaxed, glide_path_relaxed) = sf_calc(bulk, n_steps_path_1, 0.4, (0.0, 0.0), (0.33333333, 0.6666666), relax_from_configs=relax_from_configs )
    SFE_unrelaxed = glide_path_unrelaxed[-1][3]
    SFE_relaxed = glide_path_relaxed[-1][3]
    # (glide_path_unrelaxed_2, glide_path_relaxed_2) = sf_calc(bulk, n_steps_path_2, 0.4, (0.33333333, 0.6666666), (0.0, 1.0) )
    # xo = glide_path_unrelaxed[-1][0]
    # glide_path_unrelaxed.extend([ (x[0]+xo, x[1], x[2], x[3]) for x in glide_path_unrelaxed_2] )
    # xo = glide_path_relaxed[-1][0]
    # glide_path_relaxed.extend([ (x[0]+xo, x[1], x[2], x[3]) for x in glide_path_relaxed_2] )
    file_unrelaxed = open("shuffle_relaxation_start.xyz", "w")
    file_relaxed = open("shuffle_relaxed.xyz", "w")
    try:
        relax_from_configs = ase.io.read("../model-{}-test-stacking-fault-surface/shuffle_relaxed.xyz".format("GAP-6"),":")
    except:
        relax_from_configs = None
    (shuffle_path_unrelaxed, shuffle_path_relaxed) = sf_calc(bulk, n_steps_path_1, -1.17, (0.0, 0.0), (1.0, 0.0) )


if do_surfs:
    print "do glide"
    file_unrelaxed = open("glide_relaxation_start.xyz", "w")
    file_relaxed = open("glide_relaxed.xyz", "w")
    (glide_surf_unrelaxed, glide_surf_relaxed) = sf_calc(bulk, n_steps, 0.4)
    file_unrelaxed.close()
    file_relaxed.close()

    print "do shuffle"
    file_unrelaxed = open("shuffle_relaxation_start.xyz", "w")
    file_relaxed = open("shuffle_relaxed.xyz", "w")
    (shuffle_surf_unrelaxed, shuffle_surf_relaxed) = sf_calc(bulk, n_steps, -1.17)
    file_unrelaxed.close()
    file_relaxed.close()

    if not do_paths:
        for (x_i, x_j, E) in glide_surf_unrelaxed.reshape( (-1,3) ):
            if np.abs(x_i-0.33333) < 1.0e-2 and np.abs(x_j-0.66666) < 1.0e-2:
                SFE_unrelaxed = E
        for (x_i, x_j, E) in glide_surf_relaxed.reshape( (-1,3) ):
            if np.abs(x_i-0.33333) < 1.0e-2 and np.abs(x_j-0.66666) < 1.0e-2:
                SFE_relaxed = E

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = {}
try:
    properties = { 'stacking_fault_surface_1' : [ x[0] for x in glide_surf_unrelaxed.reshape( (-1,3) ) ], 'stacking_fault_surface_2' : [ x[1] for x in glide_surf_unrelaxed.reshape( (-1,3) ) ],
                   'glide_stacking_fault_surface' : [ x[2] for x in glide_surf_relaxed.reshape( (-1,3) ) ],
                   'glide_stacking_fault_surface_unrelaxed' : [ x[2] for x in glide_surf_unrelaxed.reshape( (-1,3) ) ],
                   'shuffle_stacking_fault_surface' : [ x[2] for x in shuffle_surf_relaxed.reshape( (-1,3) ) ],
                   'shuffle_stacking_fault_surface_unrelaxed' : [ x[2] for x in shuffle_surf_unrelaxed.reshape( (-1,3) ) ],
                   'stacking_fault_energy': SFE_relaxed, 'stacking_fault_energy_unrelaxed' : SFE_unrelaxed}
except:
    pass

try:
    properties.update({ 'glide_stacking_fault_path' : glide_path_relaxed, 'glide_stacking_fault_path_unrelaxed' : glide_path_unrelaxed,
                        'shuffle_stacking_fault_path' : shuffle_path_relaxed, 'shuffle_stacking_fault_path_unrelaxed' : shuffle_path_unrelaxed,
                        'stacking_fault_energy': SFE_relaxed, 'stacking_fault_energy_unrelaxed' : SFE_unrelaxed})
except:
    pass
