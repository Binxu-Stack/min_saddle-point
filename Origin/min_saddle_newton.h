#ifdef MINIMIZE_CLASS

MinimizeStyle(saddle_newton,MinSaddleNewton)

#else

#ifndef MINSADDLE_NEWTON_H
#define MINSADDLE_NEWTON_H

#include "min_saddle_sd.h"

using namespace std;

namespace LAMMPS_NS {

class MinSaddleNewton: public MinSaddleSD {
public:
  MinSaddleNewton(class LAMMPS *);
  int iterate(int);

};
}
#endif
#endif
