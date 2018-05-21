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
import ase.units
import numpy as np
import quippy
import re
import os

import __builtin__

# the current model
import model 

a0 = 5.43 # initial guess at lattice constant, cell will be relaxed below
N = 2 # number of unit cells in each direction

# set up the a
bulk = Diamond(symbol='Si', latticeconstant=a0)
bulk *= (N, N, N)
bulk = quippy.Atoms(bulk)

def save_hook(traj, a=bulk):
   traj.append(a.copy())

quippy.system_set_random_seeds(10)

debug = True
save_traj = True
n_steps_heat=20000
n_steps_equil=10000
n_steps_liquid=5000
T_melt = 5000.0
T_equil = 2000.0

liquid_traj_filename="model-"+model.name+"-test-liquid-liquid_traj.xyz"

if __builtin__.do_io:
    print liquid_traj_filename

try:
   liquid_traj=quippy.AtomsList(liquid_traj_filename)
   if __builtin__.do_io:
       print "read %s from file" % liquid_traj_filename
except:
   if __builtin__.do_io:
       print "generating %s" % liquid_traj_filename
   if 'LIQUID_NO_RUN' in os.environ:
       sys.exit(1)

   # specify that we will use model.calculator to compute forces, energies and stresses
   bulk.set_calculator(model.calculator)
   print "doing bulk"
   print bulk.get_potential_energy()
   print "done"
   try:
       model.calculator.print_()
   except:
       pass

   ####################################################################################################
   # initial heating
   dyn_heat = quippy.Dynamics(bulk, 0.5*ase.units.fs, trajectory=None, initialtemperature=T_melt, loginterval=10, logfile="model-"+model.name+"-test-liquid-heat.log")
   dyn_heat.add_thermostat(type=quippy.THERMOSTAT_ALL_PURPOSE, T=T_melt, tau=10.0)
   dyn_heat.set_barostat(type=quippy.BAROSTAT_HOOVER_LANGEVIN, p_ext=0.0, hydrostatic_strain=True, diagonal_strain=True, finite_strain_formulation=True, tau_epsilon=200.0, T=T_melt, W_epsilon_factor=5.0)

   if debug:
      traj=quippy.AtomsList()
      dyn_heat.attach(save_hook, 10, traj=traj)

   dyn_heat.run(steps=n_steps_heat)

   if debug and __builtin__.do_io:
      traj.write("model-"+model.name+"-liquid-heat_traj.xyz")
   ####################################################################################################

   ####################################################################################################
   # equilibration heating
   dyn_equil = quippy.Dynamics(bulk, 0.25*ase.units.fs, trajectory=None, loginterval=10, logfile="model-"+model.name+"-test-liquid-equil.log")
   dyn_equil._ds.rescale_velo(T_equil)

   # quippy.system.set_random_seeds(10)
   dyn_equil.add_thermostat(type=quippy.THERMOSTAT_ALL_PURPOSE, T=T_equil, tau=50.0)
   dyn_equil.set_barostat(type=quippy.BAROSTAT_HOOVER_LANGEVIN, p_ext=0.0, hydrostatic_strain=True, diagonal_strain=True, finite_strain_formulation=True, tau_epsilon=200.0, T=T_equil, W_epsilon_factor=5.0)

   if debug and __builtin__.do_io:
      traj=quippy.AtomsList()
      dyn_equil.attach(save_hook, 10, traj=traj)

   dyn_equil.run(steps=n_steps_equil)

   if debug and __builtin__.do_io:
      traj.write("model-"+model.name+"-liquid-equil_traj.xyz")
   ####################################################################################################

   ####################################################################################################
   # collect data
   dyn_liquid = quippy.Dynamics(bulk, 0.25*ase.units.fs, trajectory=None, loginterval=10, logfile="model-"+model.name+"-test-liquid-data.log")

   # quippy.system.set_random_seeds(10)
   dyn_liquid.add_thermostat(type=quippy.THERMOSTAT_ALL_PURPOSE, T=T_equil, tau=200.0)
   dyn_liquid.set_barostat(type=quippy.BAROSTAT_HOOVER_LANGEVIN, p_ext=0.0, hydrostatic_strain=True, diagonal_strain=True, finite_strain_formulation=True, tau_epsilon=500.0, T=T_equil, W_epsilon_factor=5.0)

   liquid_traj=quippy.AtomsList()
   dyn_liquid.attach(save_hook, 10, traj=liquid_traj)

   dyn_liquid.run(steps=n_steps_liquid)

   if save_traj and __builtin__.do_io:
      liquid_traj.write(liquid_traj_filename)
   ####################################################################################################

def calc_rdf(traj, r_min, r_max, sigma, n_bins):
   rdfd = np.zeros( (n_bins) )
   rdfd_1 = np.zeros( (n_bins) )
   bin_width = (r_max-r_min)/(n_bins-1)
   r_vals = np.zeros( (n_bins) )
   for at in traj:
      quippy.rdfd_calc(rdfd_1, at, zone_center=np.zeros( (3) ), zone_atom_center=-1, bin_width=bin_width, n_bins=n_bins,
	 zone_width=-1.0, n_zones=1, gaussian_smoothing=True, gaussian_sigma=sigma, center_mask_str="", neighbour_mask_str="", 
	 bin_pos=r_vals)
      rdfd[:] += rdfd_1
   return zip(r_vals, rdfd/len(traj))

def calc_adf(traj, r_cut, n_bins):
   adfd = np.zeros( (n_bins) )
   adfd_1 = np.zeros( (n_bins) )
   a_vals = np.zeros( (n_bins) )
   for at in traj:
      quippy.adfd_calc(adfd_1, at, zone_center=np.zeros( (3) ), n_angle_bins=n_bins, dist_bin_width=r_cut, n_dist_bins=1, 
	 zone_width=-1.0, n_zones=1, center_mask_str="", neighbour_1_mask_str="", neighbour_1_max_dist=r_cut, neighbour_2_mask_str="", 
	 dist_bin_rc2=True, angle_bin_pos=a_vals)
      adfd[:] += adfd_1
   # do mean
   adfd /= len(traj)
   # normalize to 1
   adfd *= (float(n_bins)/180.0)/sum(adfd)
   return zip(a_vals, adfd)


# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = {'liquid_rdf': calc_rdf(liquid_traj, 1.0, 7.0, 0.1, 81),
              'liquid_adf': calc_adf(liquid_traj, 2.8, 51) }
