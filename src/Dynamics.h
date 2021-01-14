#ifndef DYNAMICS_H
#define DYNAMICS_H

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
#include<cfloat>
#include<cmath>
#include<adolc/adolc.h>
#include<ilcplex/ilocplex.h>

#include "Data.h"
#include "Util.h"

using namespace std;

class Dynamics{

   public:
      Dynamics(Data*);
      ~Dynamics();

      //dimensions
      static const int stateSize = 6;//state vector size
      static const int controlSize = 2;//control vector size
      //sytem data
      int legArrayNumRows;
      int legArrayNumCols;
      class dynamics{
         friend class Dynamics;
         private:
         bool status;
         double* Yeq;
         double* Ueq;
         double** stateMatrix;
         double** controlMatrix;
         double stateMatrixMaxNorm;
         double ctrlMatrixMaxNorm;
         public:
         int orig;
         int dest;
         double* Yo;
         double* Uo;
         double* Yf;
         double* Uf;
         double* vo;
         double fRadius;
         double fzMin;
         double fzMax;
         double deltaTa;
         double deltaTb;
         double get_Yeq(int) const;
         double get_Ueq(int) const;
         double get_A(int,int) const;
         double get_B(int,int) const;
         double getMaxNormA() const;
         double getMaxNormB() const;
         void get_leg(int&,int&,bool&) const;
         bool getStatus() const;
         bool isLanding;
      };
      friend class dynamics;
      dynamics** legDyn;
      double** flightTimes;
      double** flightDist;
      double** I;
      //glider data
      double CD0;//coefficient of drag at zero-lift (dimensionless)
      double K;//aerodynamic coefficient (dimensionless)
      double gravity;//(m/s^2)
      double mass;//(kg)
      double S;//(m^2)
      double alpha;//camera openning angle (degrees->rad)
      double vLb;
      double gammaLb;//(degrees->rad)
      double phiLb;//(degrees->rad)
      double CLlb;
      double muLb;//(degrees->rad)
      double vUb;
      double gammaUb;//(degrees->rad)
      double phiUb;//(degrees->rad)
      double CLub;
      double muUb;//(degrees->rad)
      //upper bounds on derivatives of states and control
      double xDot_ub;
      double yDot_ub;
      double hDot_ub;
      double vDot_ub;
      double gammaDot_ub;
      double phiDot_ub;
      double CLDot_ub;
      double muDot_ub;//(degrees->rad)
      double yDotUbMaxNorm;
      double uDotUbMaxNorm;
      //wind data
      double rho;//density of the air at sea level (kg/m3)
      double beta;// wind strength (1/s);
      //other data
      double PI;
      short int tag; //parameter for ADOL-C
      double safetyDist;

      //public functions
      double windSpeed(double);
      dynamics getLeg(int,int);
      bool checkRouteStatus(vector<int>&,int);
      void printFlightMatrices();

   private:
      Data* dataPtr;
      void steadyState(dynamics*);
      void gliderODE_ADOL(dynamics);
      void ODE_Jac(dynamics*);
      void analyticJac(dynamics);
      void computeDynamics();
      void computeFlightTimes();
      void randomFlightTimes();
      void computeFlightDist();
};
#endif










