/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Aidan Thompson (SNL)
                        improved CG and backtrack ls, added quadratic ls
   Sources: Numerical Recipes frprmn routine
            "Conjugate Gradient Method Without the Agonizing Pain" by
            JR Shewchuk, http://www-2.cs.cmu.edu/~jrs/jrspapers.html#cg
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "min_saddle.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "fix_minimize.h"
#include "compute.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "output.h"
#include "thermo.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

MinSaddle::MinSaddle(LAMMPS *lmp) : Min(lmp)
{
  vector_exist = 3;
}

/* ---------------------------------------------------------------------- */

MinSaddle::~MinSaddle()
{
}

/* ---------------------------------------------------------------------- */

void MinSaddle::init()
{
  Min::init();
  // register two atom vector in fix_minimize
  fix_minimize->add_vector(3);
  fix_minimize->add_vector(3);
}



/* ---------------------------------------------------------------------- */



/* ----------------------------------------------------------------------
   evaluate potential energy and forces
   may migrate atoms due to reneighboring
   return new energy, which should include nextra_global dof
   return negative gradient stored in atom->f
   return negative gradient for nextra_global dof in fextra
------------------------------------------------------------------------- */
double MinSaddle::my_energy_force(int resetflag)
{
  // check for reneighboring
  // always communicate since minimizer moved atoms

  int nflag = neighbor->decide();

  if (nflag == 0) {
    timer->stamp();
    comm->forward_comm();
    timer->stamp(Timer::COMM);
  } else {
    if (modify->n_min_pre_exchange) {
      timer->stamp();
      modify->min_pre_exchange();
      timer->stamp(Timer::MODIFY);
    }
    if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    if (domain->box_change) {
      domain->reset_box();
      comm->setup();
      if (neighbor->style) neighbor->setup_bins();
    }
    timer->stamp();
    comm->exchange();
    if (atom->sortfreq > 0 &&
        update->ntimestep >= atom->nextsort) atom->sort();
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    timer->stamp(Timer::COMM);
    if (modify->n_min_pre_neighbor) {
      modify->min_pre_neighbor();
      timer->stamp(Timer::MODIFY);
    }
    neighbor->build(1);
    timer->stamp(Timer::NEIGH);
    if (modify->n_min_post_neighbor) {
      modify->min_post_neighbor();
      timer->stamp(Timer::MODIFY);
    }
  }

  ev_set(update->ntimestep);
  force_clear();

  timer->stamp();

  if (modify->n_min_pre_force) {
    modify->min_pre_force(vflag);
    timer->stamp(Timer::MODIFY);
  }

  if (pair_compute_flag) {
    force->pair->compute(eflag,vflag);
    timer->stamp(Timer::PAIR);
  }

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
    timer->stamp(Timer::BOND);
  }

  if (kspace_compute_flag) {
    force->kspace->compute(eflag,vflag);
    timer->stamp(Timer::KSPACE);
  }

  if (modify->n_min_pre_reverse) {
    modify->min_pre_reverse(eflag,vflag);
    timer->stamp(Timer::MODIFY);
  }

  if (force->newton) {
    comm->reverse_comm();
    timer->stamp(Timer::COMM);
  }

  // update per-atom minimization variables stored by pair styles

  if (nextra_atom)
    for (int m = 0; m < nextra_atom; m++)
      requestor[m]->min_xf_get(m);

  // fixes that affect minimization

  if (modify->n_min_post_force) {
     timer->stamp();
     modify->min_post_force(vflag);
     timer->stamp(Timer::MODIFY);
  }

  // compute potential energy of system
  // normalize if thermo PE does

  double energy = pe_compute->compute_scalar();
  if (nextra_global) energy += modify->min_energy(fextra);
  if (output->thermo->normflag) energy /= atom->natoms;

  // if reneighbored, atoms migrated
  // if resetflag = 1, update x0 of atoms crossing PBC
  // reset vectors used by lo-level minimizer

  if (nflag) {
    if (resetflag) fix_minimize->reset_coords();
    reset_vectors();
  }

  return energy;
}


double MinSaddle::energy_force(int resetflag)
{
  // check for reneighboring
  // always communicate since minimizer moved atoms

  int nflag = neighbor->decide();

  if (nflag == 0) {
    timer->stamp();
    comm->forward_comm();
    timer->stamp(Timer::COMM);
  } else {
    if (modify->n_min_pre_exchange) {
      timer->stamp();
      modify->min_pre_exchange();
      timer->stamp(Timer::MODIFY);
    }
    if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    if (domain->box_change) {
      domain->reset_box();
      comm->setup();
      if (neighbor->style) neighbor->setup_bins();
    }
    timer->stamp();
    comm->exchange();
    if (atom->sortfreq > 0 &&
        update->ntimestep >= atom->nextsort) atom->sort();
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    timer->stamp(Timer::COMM);
    if (modify->n_min_pre_neighbor) {
      modify->min_pre_neighbor();
      timer->stamp(Timer::MODIFY);
    }
    neighbor->build(1);
    timer->stamp(Timer::NEIGH);
    if (modify->n_min_post_neighbor) {
      modify->min_post_neighbor();
      timer->stamp(Timer::MODIFY);
    }
  }

  ev_set(update->ntimestep);
  force_clear();

  timer->stamp();

  if (modify->n_min_pre_force) {
    modify->min_pre_force(vflag);
    timer->stamp(Timer::MODIFY);
  }

  if (pair_compute_flag) {
    force->pair->compute(eflag,vflag);
    timer->stamp(Timer::PAIR);
  }

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
    timer->stamp(Timer::BOND);
  }

  if (kspace_compute_flag) {
    force->kspace->compute(eflag,vflag);
    timer->stamp(Timer::KSPACE);
  }

  if (modify->n_min_pre_reverse) {
    modify->min_pre_reverse(eflag,vflag);
    timer->stamp(Timer::MODIFY);
  }

  if (force->newton) {
    comm->reverse_comm();
    timer->stamp(Timer::COMM);
  }

  // update per-atom minimization variables stored by pair styles

  if (nextra_atom)
    for (int m = 0; m < nextra_atom; m++)
      requestor[m]->min_xf_get(m);

  // fixes that affect minimization

  if (modify->n_min_post_force) {
     timer->stamp();
     modify->min_post_force(vflag);
     timer->stamp(Timer::MODIFY);
  }

  // compute potential energy of system
  // normalize if thermo PE does

  //double energy = pe_compute->compute_scalar();
  //if (nextra_global) energy += modify->min_energy(fextra);
  //if (output->thermo->normflag) energy /= atom->natoms;

  // if reneighbored, atoms migrated
  // if resetflag = 1, update x0 of atoms crossing PBC
  // reset vectors used by lo-level minimizer

  if (nflag) {
    if (resetflag) fix_minimize->reset_coords();
    reset_vectors();
  }
  // reset force and energy
  const double EPS = 1e-10;
  double iEPS = 1.0/EPS;

  // new energy
  // energy = 0.5 * |F|^2
  double fnorm2 = fnorm_sqr();
  //double ifnorm2 = 1.0 / fnorm2;
  double fnorm = sqrt(fnorm2);
  double ifnorm = 1.0/fnorm;
  double energy = fnorm2 * 0.5;

  // new force
  // f_new = H \cdot f
  
  // store  initial coordinate and force
  double *f0;
  double *x0;

  //memory->create(f0,nvec, "MINSADDLE:f0");
  f0 = fix_minimize->request_vector(vector_exist);
  x0 = fix_minimize->request_vector(vector_exist+1);
  for (int i = 0; i < nvec; ++i) {
    x0[i] = xvec[i];
    f0[i] = fvec[i];
  }
  // move along the force direction
  for (int i = 0; i < nvec; ++i) {
    xvec[i] = x0[i] + fvec[i] * ifnorm * EPS;
  }
  my_energy_force(1);
  f0 = fix_minimize->request_vector(vector_exist);
  x0 = fix_minimize->request_vector(vector_exist+1);
  for (int i = 0; i < nvec; ++i) {
    fvec[i] = (f0[i]-fvec[i])*iEPS*fnorm;
  }
  // restore to initial position
  for (int i = 0; i < nvec; ++i) {
    xvec[i] = x0[i];
  }
  //my_energy_force(1);

  //printf("Energy: %g\n", energy);

  return energy;
}


