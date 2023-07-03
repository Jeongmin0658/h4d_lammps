/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(ewald/neq,EwaldNeq);
// clang-format on
#else

#ifndef LMP_EWALD_NEQ_H
#define LMP_EWALD_NEQ_H

#include "kspace.h"

namespace LAMMPS_NS {

class EwaldNeq : public KSpace {
 public:
  EwaldNeq(class LAMMPS *);
  virtual ~EwaldNeq();
  void init();
  void setup();
  virtual void settings(int, char **);
  virtual void compute(int, int);
  double memory_usage();

  void compute_group_group(int, int, int);

 protected:
  bigint time_origin_neq, time_final_neq, time_diff_neq;    // JK: 
  int    isNeqinsertion, isNeqconvex;                       // JK:
  double vspeed_neq, max_altitude_neq, power_neq;           // JK: for altitude dynamics in 4 dim.
  double max2_altitude_neq;
  double lambda_w,altitude_screen;                          // JK: altitude screening
  int    add_type_neq, add2_type_neq;                       // JK: for altitude dynamics in 4 dim.
  int kxmax, kymax, kzmax;
  int kcount, kmax, kmax3d, kmax_created;
  double gsqmx, volume;
  int nmax;

  double unitk[3];
  int *kxvecs, *kyvecs, *kzvecs;
  int kxmax_orig, kymax_orig, kzmax_orig;
  double *ug;
  double **eg, **vg;
  double **ek;
  double *sfacrl, *sfacim, *sfacrl_all, *sfacim_all;
  double ***cs, ***sn;

  // screening
  virtual void vertical_screening();

  // group-group interactions

  int group_allocate_flag;
  double *sfacrl_A, *sfacim_A, *sfacrl_A_all, *sfacim_A_all;
  double *sfacrl_B, *sfacim_B, *sfacrl_B_all, *sfacim_B_all;

  double rms(int, double, bigint, double);
  virtual void eik_dot_r();
  void coeffs();
  virtual void allocate();
  void deallocate();
  void slabcorr();

  // triclinic

  int triclinic;
  void eik_dot_r_triclinic();
  void coeffs_triclinic();

  // group-group interactions

  void slabcorr_groups(int, int, int);
  void allocate_groups();
  void deallocate_groups();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use Ewald with 2d simulation

The kspace style ewald cannot be used in 2d simulations.  You can use
2d Ewald in a 3d simulation; see the kspace_modify command.

E: Kspace style requires atom attribute q

The atom style defined does not have these attributes.

E: Cannot use non-periodic boundaries with Ewald

For kspace style ewald, all 3 dimensions must have periodic boundaries
unless you use the kspace_modify command to define a 2d slab with a
non-periodic z dimension.

E: Incorrect boundaries with slab Ewald

Must have periodic x,y dimensions and non-periodic z dimension to use
2d slab option with Ewald.

E: Cannot (yet) use Ewald with triclinic box and slab correction

This feature is not yet supported.

E: KSpace style is incompatible with Pair style

Setting a kspace style requires that a pair style with matching
long-range Coulombic or dispersion components be used.

E: KSpace accuracy must be > 0

The kspace accuracy designated in the input must be greater than zero.

E: Must use 'kspace_modify gewald' for uncharged system

UNDOCUMENTED

E: Cannot (yet) use K-space slab correction with compute group/group for triclinic systems

This option is not yet supported.

*/