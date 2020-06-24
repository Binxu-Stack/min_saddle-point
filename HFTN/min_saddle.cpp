/* -------------------------------------------------------------------------------------------------
 * ARTn: Activation Relaxation Technique nouveau
 * Bin Xu, xubinrun@gmail.com; Lingti Kong, konglt@gmail.com
------------------------------------------------------------------------------------------------- */
#include "min_saddle.h"

//#define DEBUG
#include <cmath>
#include <cstdlib>
#include <cstring>
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



/* -------------------------------------------------------------------------------------------------
 * lapack or MKL-lapack is used to evaluate the lowest eigenvalue of the matrix in Lanczos.
------------------------------------------------------------------------------------------------- */
//#ifdef MKL
//#include "mkl.h"
//#define dstev_  dstev
//#else
//extern "C" {
//extern void dstev_(char *, int*, double *, double *, double *, int *, double *, int *);
//};
//#endif
//
#define EPS_ENERGY 1.e-8
//
//enum{MAXITER,MAXEVAL,ETOL,FTOL,DOWNHILL,ZEROALPHA,ZEROFORCE,ZEROQUAD};



/* -------------------------------------------------------------------------------------------------
 * Constructor of ARTn
------------------------------------------------------------------------------------------------- */
MinSaddle::MinSaddle(LAMMPS *lmp): MinHFTN(lmp) {
}

/* -------------------------------------------------------------------------------------------------
 * Reconstruct the energy force from min class
------------------------------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   evaluate potential energy and forces
   may migrate atoms due to reneighboring
   return new energy, which should include nextra_global dof
   return negative gradient stored in atom->f
   return negative gradient for nextra_global dof in fextra
------------------------------------------------------------------------- */

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
  
  // reset force and energy
  const double EPS = 1e-12;
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
  //double *f0;
  for (int i = 0; i < nvec; ++i) {
    x0[i] = xvec[i];
    f0[i] = fvec[i];
  }
  // move along the force direction
  for (int i = 0; i < nvec; ++i) {
    xvec[i] = x0[i] + fvec[i] * ifnorm * EPS;
  }
  my_energy_force(1);
  for (int i = 0; i < nvec; ++i) {
    fvec[i] = (f0[i]-fvec[i])*iEPS*fnorm;
  }
  // restore to initial position
  for (int i = 0; i < nvec; ++i) {
    xvec[i] = x0[i];
  }
  //my_energy_force(1);



  //double energy = pe_compute->compute_scalar();
  //if (nextra_global) energy += modify->min_energy(fextra);
  //if (output->thermo->normflag) energy /= atom->natoms;

  // if reneighbored, atoms migrated
  // if resetflag = 1, update x0 of atoms crossing PBC
  // reset vectors used by lo-level minimizer

  //if (nflag) {
  //  if (resetflag) fix_minimize->reset_coords();
  //  reset_vectors();
  //}


  return energy;
}

double MinSaddle::my_energy_force(int resetflag){
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

void MinSaddle::setup_style()
{
  MinHFTN::setup_style();

  fix_minimize->add_vector(3); // x0
  fix_minimize->add_vector(3); // f0
  return;
}

void MinSaddle::reset_vectors()
{
  MinHFTN::reset_vectors();
  x0 = fix_minimize->request_vector(NUM_HFTN_ATOM_BASED_VECTORS);
  f0 = fix_minimize->request_vector(NUM_HFTN_ATOM_BASED_VECTORS+1);

  return;
}


