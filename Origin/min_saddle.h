/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef MINSADDLE_H
#define MINSADDLE_H



//#include "pointers.h"
#include "min.h"

namespace LAMMPS_NS {

class MinSaddle : public Min {
 public:
  MinSaddle(class LAMMPS *);
  virtual ~MinSaddle();
  int vector_exist;
  virtual void init();
  // possible return values of iterate() method
 protected:
  //double energy_force(int);
  double energy_force(int);
  double my_energy_force(int);
};

}

#endif

/* ERROR/WARNING messages:

W: Using 'neigh_modify every 1 delay 0 check yes' setting during minimization

UNDOCUMENTED

E: Minimization could not find thermo_pe compute

This compute is created by the thermo command.  It must have been
explicitly deleted by a uncompute command.

E: Cannot use a damped dynamics min style with fix box/relax

This is a current restriction in LAMMPS.  Use another minimizer
style.

E: Cannot use a damped dynamics min style with per-atom DOF

This is a current restriction in LAMMPS.  Use another minimizer
style.

E: Cannot use hftn min style with fix box/relax

UNDOCUMENTED

E: Cannot use hftn min style with per-atom DOF

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

U: Resetting reneighboring criteria during minimization

Minimization requires that neigh_modify settings be delay = 0, every =
1, check = yes.  Since these settings were not in place, LAMMPS
changed them and will restore them to their original values after the
minimization.

U: Energy due to X extra global DOFs will be included in minimizer energies

When using fixes like box/relax, the potential energy used by the minimizer
is augmented by an additional energy provided by the fix. Thus the printed
converged energy may be different from the total potential energy.

*/
