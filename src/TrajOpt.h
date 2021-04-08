#ifndef TRAJOPT_H
#define TRAJOPT_H

#include<iomanip>
#include<iostream>
#include<fstream>
#include<vector>
#include<list>
#include<cstdlib>
#include<cstring>
#include<climits>
#include<utility>
#include<algorithm>// std::sort
#include<cmath>
#include<cassert>

#include "Util.h"
#include "Dynamics.h"
#include "Data.h"
#include "Unconstrained.h"
#include "CPLEXmodel.h"
#include "AMPLmodel.h"

using namespace std;

class TrajOpt{

   public:
      TrajOpt(Data*, Dynamics*, CPLEXmodel*, AMPLmodel*, string);
      ~TrajOpt();
       //IloAlgorithm::Status findTraj(int*, int, vector<vector<double>>&, double*, double*, double*);
       IloAlgorithm::Status findTraj(vector<int>&, int,
             vector<vector<double>>&, double*, double*, double*, 
             int&, double**, double**, double**);

   private:
      Dynamics* dynPtr;
      Data* data;
      CPLEXmodel* model;
      AMPLmodel* nlp;
      string method;
};
#endif
