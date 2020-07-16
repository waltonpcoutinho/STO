#ifndef DATA_H
#define DATA_H

#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cfloat>

using namespace std;

class Data{

   friend class Solution;
   public:
      Data(int, char**);
      ~Data();

      int T;//discretisation size
      int maxFleetSize;
      int fleetSize;
      double t0;
      double tf;
      double PI;
      int noWaypoints;
      int noLandPoints;
      string fileName;
      string instType;
      //instance data
      double* depot;
      double* xBar;
      double* yBar;
      double* zBar;
      double* radiusBar;
      double* minZBar;
      double* maxZBar;
      double* xTilde;
      double* yTilde;
      double* zTilde;
      double* radiusTilde;
      //distance matrices
      double* distFromDepot;
      double** distMatrix;
      double** distToLanding;
      //bounds over variables
      double xLb;
      double yLb;
      double zLb;
      double xUb;
      double yUb;
      double zUb;

   private:
      void computeConstants();
      void getInstName(char* instPath);
};
#endif
