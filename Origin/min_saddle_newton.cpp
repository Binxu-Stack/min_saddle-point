/* -------------------------------------------------------------------------------------------------
 * ARTn: Activation Relaxation Technique nouveau
 * Bin Xu, xubinrun@gmail.com; Lingti Kong, konglt@gmail.com
------------------------------------------------------------------------------------------------- */
#include "min_saddle_newton.h"

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
MinSaddleNewton::MinSaddleNewton(LAMMPS *lmp): MinSaddleSD(lmp) {
}


int MinSaddleNewton::iterate(int maxiter)
{
  int i,m,n,fail,ntimestep;
  double fdotf;
  double *fatom,*hatom;

  // test
  double step = 0.1;

  // initialize working vectors
  ecurrent = energy_force(0);

  for (i = 0; i < nvec; i++) h[i] = fvec[i];
  if (nextra_atom)
    for (m = 0; m < nextra_atom; m++) {
      fatom = fextra_atom[m];
      hatom = hextra_atom[m];
      n = extra_nlen[m];
      for (i = 0; i < n; i++) hatom[i] = fatom[i];
    }
  if (nextra_global)
    for (i = 0; i < nextra_global; i++) hextra[i] = fextra[i];

  // add new vector
  fix_minimize->add_vector(3);
  double * x0tmp;
  for (int iter = 0; iter < maxiter; iter++) {

    if (timer->check_timeout(niter))
      return TIMEOUT;

    ntimestep = ++update->ntimestep;
    niter++;

    // line minimization along h from current position x
    // h = downhill gradient direction

    eprevious = ecurrent;

    // test
    //fail = (this->*linemin)(ecurrent,alpha_final);
    //if (fail) return fail;
    
    double alpha, alpha_max;
    double hme = 0.0;
    double hmaxall, alphamax;
    hme = 0.0;
    for (i = 0; i < nvec; i++)
      hme = MAX(hme,fabs(h[i]));

    MPI_Allreduce(&hme,&hmaxall,1,MPI_DOUBLE,MPI_MAX,world);
    alpha_max = dmax/hmaxall;
    alpha = MIN(alpha_max,step);
    x0tmp = fix_minimize->request_vector(5);
    //for( i = 0; i < nvec; ++i) {
    //  if(atom->tag[i]==4001){
    //    printf("Index: %i\n", i);
    //  }
    //}
    for (i = 0; i < nvec; ++i) {
      x0tmp[i] = xvec[i];
    }
    for (i = 0; i < nvec; i++) xvec[i] += alpha*h[i];
    //ecurrent = alpha_step(alpha, 1);
    ecurrent = energy_force(0);
    int refuse_flag;
    refuse_flag = 0;
    if (ecurrent < eprevious) {
      step *= 1.2;
    } else {
      x0tmp = fix_minimize->request_vector(5);
      for (i = 0; i < nvec; ++i) xvec[i] = x0tmp[i];
      step *= 0.5;
      ecurrent = energy_force(0);
      refuse_flag = 1;
    }
    
    

    // function evaluation criterion

    if (neval >= update->max_eval) return MAXEVAL;

    // energy tolerance criterion

    if (refuse_flag != 1 && fabs(ecurrent-eprevious) <
        update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
      return ETOL;

    // force tolerance criterion

    fdotf = fnorm_sqr();
    if (fdotf < update->ftol*update->ftol) return FTOL;

    // set new search direction h to f = -Grad(x)

    for (i = 0; i < nvec; i++) h[i] = fvec[i];
    if (nextra_atom)
      for (m = 0; m < nextra_atom; m++) {
        fatom = fextra_atom[m];
        hatom = hextra_atom[m];
        n = extra_nlen[m];
        for (i = 0; i < n; i++) hatom[i] = fatom[i];
      }
    if (nextra_global)
      for (i = 0; i < nextra_global; i++) hextra[i] = fextra[i];

    // output for thermo, dump, restart files

    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }

  return MAXITER;
}


