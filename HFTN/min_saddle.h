#ifdef MINIMIZE_CLASS

MinimizeStyle(saddle,MinSaddle)

#else

#ifndef MINSADDLE_H
#define MINSADDLE_H

#include "min_hftn.h"

using namespace std;

namespace LAMMPS_NS {

class MinSaddle: public MinHFTN {
public:
  MinSaddle(class LAMMPS *);
  void setup_style();
  void reset_vectors();

private:
  enum {
    VEC_XK=0,   //-- ATOM POSITIONS AT SUBITER START
    VEC_CG_P,   //-- STEP p IN CG SUBITER
    VEC_CG_D,   //-- DIRECTION d IN CG SUBITER
    VEC_CG_HD,  //-- HESSIAN-VECTOR PRODUCT Hd
    VEC_CG_R,   //-- RESIDUAL r IN CG SUBITER
    VEC_DIF1,   //-- FOR FINITE DIFFERENCING
    VEC_DIF2,   //-- FOR FINITE DIFFERENCING
    NUM_HFTN_ATOM_BASED_VECTORS
  };


  double energy_force(int);
  double * f0;
  double * x0;
  double my_energy_force(int);
};
}
#endif
#endif
