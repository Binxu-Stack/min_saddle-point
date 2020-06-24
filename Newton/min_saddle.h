#ifdef MINIMIZE_CLASS

MinimizeStyle(saddle,MinSaddle)

#else

#ifndef MINSADDLE_H
#define MINSADDLE_H

#include "min_sd.h"

using namespace std;

namespace LAMMPS_NS {

class MinSaddle: public MinSD {
public:
  MinSaddle(class LAMMPS *);
  int iterate(int);

private:
  double energy_force(int);
  double my_energy_force(int);
};
}
#endif
#endif
