#ifndef SOLUTION_H
#define SOLUTION_H

#include <cstdlib>
#include <cmath>
#include <list>
#include <vector>
#include <algorithm>
#include <iterator>
#include <limits>

#include "Data.h"

using namespace std;

struct glider{
   int maxRouteSize;
   int infOrig;
   bool isFeasible;
   double routeCost;
   vector<int> route;
   double* flightTimes;
   double* stepSizes;
   double* errors; //errors per arc
   vector<int> collisionList;
   vector<vector<double>> trajectory;
   double** normEpsilons; //errors per time step for all arcs
   double** normTaylor1st; //first term of Taylor's approx.
   double** relEpsilons; //relative errors per time step and arcs
};

class Solution {

   public:
   Solution();
   ~Solution();

   //global and local solutions
   vector<glider> globalSol;
   vector<glider> localSol;

   //global and local solution costs
   double ILScost;
   double globalSolCost;
   double localSolCost;
   
   //get route cost
   double computeRoutesCost(vector<glider>&, double** flightTimes);
   
   //get solution cost
   double getTotalCost(vector<glider>&);

   //total number of LS iterations
   vector<double> computingTime;
   
   //return version
   bool isMakespan();
   int getMkspRoute(vector<glider>&);
   
   private:
   bool makespan;
};


#endif /* SOLUTION_H_ */
