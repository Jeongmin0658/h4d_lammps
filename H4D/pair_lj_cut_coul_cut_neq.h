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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/coul/cut/neq,PairLJCutCoulCutNeq);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_COUL_CUT_NEQ_H
#define LMP_PAIR_LJ_CUT_COUL_CUT_NEQ_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCutCoulCutNeq : public Pair {
 public:
  PairLJCutCoulCutNeq(class LAMMPS *);
  virtual ~PairLJCutCoulCutNeq();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  virtual double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

 protected:
  bigint time_origin_neq, time_final_neq, time_diff_neq;    // JK: 
  int    isNeqinsertion, isNeqconvex;                       // JK:
  double vspeed_neq, max_altitude_neq, power_neq;           // JK: for altitude dynamics in 4 dim.
  double max2_altitude_neq;
  double lambda_w,altitude_screen;                          // JK: altitude screening
  int    add_type_neq, add2_type_neq;                       // JK: for altitude dynamics in 4 dim.
  double cut_lj_global, cut_coul_global;
  double **cut_lj, **cut_ljsq;
  double **cut_coul, **cut_coulsq;
  double **epsilon, **sigma;
  double **lj1, **lj2, **lj3, **lj4, **offset;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style lj/cut/coul/cut requires atom attribute q

The atom style defined does not have this attribute.

*/
