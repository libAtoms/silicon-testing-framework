# tests for silicon testing framework

* force errors on databases
  * `GAP-6_teaching_database_force_error` force error on final GAP teaching database
  * `testing_database_force_error` force error on testing database gathered from various tests
* crystal structures
  * `bulk_bc8` BC8 bulk crystal
  * `bulk_bcc` BCC bulk crystal
  * `bulk_beta-Sn` beta-Sn bulk crystal
  * `bulk_diamond` diamond bulk crystal
  * `bulk_fcc` FCC bulk crystal
  * `bulk_hcp` HCP bulk crystal
  * `bulk_hcp_short_and_fat` HCP' bulk crystal
  * `bulk_hex_diamond` hexagonal diamond bulk crystal
  * `bulk_sh` simple hexagonal bulk crystal
  * `bulk_st12` ST12 bulk crystal
* other bulk tests
  * `phonon_diamond` phonons of diamond structure
  * `qha_diamond` quasi-harmonic approximation thermal expansion of diamond structure
  * `liquid` produce a liquid and calculate its RDF and ADF
* random structure search
  * `RSS_like_open_ended` random structure search
  * `gap_6_rss_8_at_CASTEP_rerelax_trajectory_evaluate` _evaluate_ the CASTEP relaxation trajectory of two final GAP RSS minima
  * `gap_6_rss_like_open_ended_minima_rerelax_sd2` re-relax minima of final GAP RSS
* point defects
  * `vacancy-energy` unrelaxed and relaxed vacancy formation energy
  * `vacancy-path-no-relax-neb` energy of pathway linearly interpolating between relaxed vacancies
  * `vacancy-path-ref` nudged-elastic band relaxation of vacancy migration
  * `interstitial-energy` relax energy of hexagonal, tetrahedral, and dumbbell interstitial
  * `fourfold-defect-small` relaxed fourfold defect in 64 atom cell
  * `diinterstitial-64-energy` relaxed diinterstitial configurations in 64+2 atom cells
* plane defects
  * `grain_boundary_112_Sigma3_tilt` relax (112) Sigma=3 tilt grain boundary
  * `stacking-fault-surface` unrelaxed and relaxed energies on pathways along hig symmetry directions on generalized stacking fault surface
  * `surface-energy-100-relaxed` relaxed (100) buckled dimer surface energy
  * `surface-energy-110-relaxed` relaxed (111) surface energy
  * `surface-energy-111-relaxed` relaxed (111) 1x1 surface energy
  * `surface-energy-111-2n+1_x_2n+1-relaxed` relaxed (111) DAS 2n+1 x 2n+1 reconstruction surface energy
  * `surface-decohesion-100-unrelaxed` unrelaxed energy pathway for decohesion at (100) surface
  * `surface-decohesion-110-unrelaxed` unrelaxed energy pathway for decohesion at (110) surface
  * `surface-decohesion-111-unrelaxed` unrelaxed energy pathway for decohesion at (111) surface
