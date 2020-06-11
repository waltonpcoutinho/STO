#ifndef CPLEXMODEL_H
#define CPLEXMODEL_H

#include<cmath>
#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<list>
#include<cstdlib>
#include<ctime>
#include<climits>
#include<algorithm>
#include<stdio.h>
#include<ilcplex/ilocplex.h>
#include<assert.h>

#include "Data.h"
#include "Dynamics.h"
#include "Util.h"

using namespace std;

class CPLEXmodel{

   public:
      CPLEXmodel(Data*,Dynamics*,int); 
      ~CPLEXmodel();
      
      IloAlgorithm::Status solveModel(const Dynamics::dynamics, double, double);
      IloCplex getCplexModel();

      double** solMatrix;
      double h;
      double objVal;
      double time;
      double error;
      int tree;
      int numIter;
      
   private:
      IloEnv env;
      IloCplex SOCP;
      IloNumVarArray x;
      IloNumVarArray y;
      IloNumVarArray z;
      IloNumVarArray vel;
      IloNumVarArray gamma;
      IloNumVarArray phi;
      IloNumVarArray CL;
      IloNumVarArray mu;
      IloNumVar epsilon;
      
      //arrays to store solution
      IloNumArray xt;
      IloNumArray yt;
      IloNumArray zt;
      IloNumArray vt;
      IloNumArray gammat;
      IloNumArray phit;
      IloNumArray CLt;
      IloNumArray mut;

      //problem dimensions
      int T;

      Data *dataPtr;
      Dynamics *gliderPtr;

      void createModel(IloModel, const Dynamics::dynamics leg, double, double);
      void finishCplex();
      void popSolMatrix(IloNumArray, IloNumArray, IloNumArray, IloNumArray, 
                        IloNumArray, IloNumArray, IloNumArray, IloNumArray);

};
#endif
