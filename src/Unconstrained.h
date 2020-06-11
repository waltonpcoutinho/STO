#ifndef UNCONSTRAINED_H
#define UNCONSTRAINED_H

#include<iomanip>
#include<iostream>
#include<fstream>
#include<vector>
#include<cstdlib>
#include<cstring>
#include<climits>
#include<cmath>
//3rd party libraries
//#include <cmaes.h>

#include "Util.h"
#include "Dynamics.h"
#include "CPLEXmodel.h"

using namespace std;
//using namespace libcmaes;

//global variables
static CPLEXmodel* modelGlobal;
static Dynamics::dynamics legGlobal;
static int statusGlobal;

class Unconstrained{

   public:
      Unconstrained();
      ~Unconstrained();

      //void callCMAES(CPLEXmodel*, Dynamics::dynamics, double, double&);
      IloAlgorithm::Status bisection(CPLEXmodel*, Dynamics::dynamics, double, double, double&, double&, double&);
      IloAlgorithm::Status minFlightTime(CPLEXmodel*, Dynamics::dynamics, double, double, double&, double&, double&);
      void plotObjEval(Data*, CPLEXmodel*, Dynamics::dynamics, double, double);

   private:

};
#endif
