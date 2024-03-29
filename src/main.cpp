//stl libraries
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <cstdlib>
#include <cstdio>
#include <cfloat>
#include <cmath>
#include <algorithm>

//internal libraries
#include "Data.h"
#include "Util.h"
#include "Dynamics.h"
#include "CPLEXmodel.h"
#include "AMPLmodel.h"
#include "Unconstrained.h"
#include "TrajOpt.h"

//#include "ILS/ILSRVND.h"
#include "Solution.h"

using namespace std;

void callTrajOpt(glider& route, Data* data, Dynamics* dyn, CPLEXmodel* model, AMPLmodel* nlp, string method);

double getRoutelb(glider uav, Dynamics* dynPtr){
   double lb = 0;
   int routeSize = uav.route.size();
   for(int k = 0; k < routeSize - 1; k++){
      const Dynamics::dynamics leg = dynPtr->getLeg(uav.route[k],uav.route[k+1]);
      lb += leg.deltaTa;
   }
   return lb;
}

int main(int argc, char** argv) 
{
   //set random seed (user function)
   setSeed();

   //load instance and glider data
   Data* dataPtr = new Data(argc, argv);
   Dynamics* gliderPtr = new Dynamics(dataPtr);
   
   //create cplex environment
   const int T = dataPtr->T;
   CPLEXmodel* modelPtr = new CPLEXmodel(dataPtr,gliderPtr,T);

   //initialise AMPL api
   AMPLmodel* nlpPtr = new AMPLmodel(dataPtr,gliderPtr,T);

   //initialise solution object
   Solution* sol = new Solution();
   if(sol->isMakespan()){
      cout << "Running MAKESPAN version!" << endl;
   }

   //initialise a glider object (single route)
   glider drone;
   int noWay = dataPtr->noWaypoints;
   int noLand = dataPtr->noLandPoints;
   //Random insertion
   //insert depot
   drone.route.push_back(0);
   //insert waypoints
   for(int i = 1; i <= noWay; i++){
      drone.route.push_back(i);
   }
   random_shuffle(drone.route.begin() + 1, drone.route.end());
   //insert landing site
   int lowerLimit = noWay + 1;
   int upperLimit = noWay + noLand;
   int landingSpot = lowerLimit + (rand() % (upperLimit - lowerLimit + 1));
   drone.route.push_back(landingSpot);

   //initialise computing time variables
   double startingTime = 0;

   //call STO-iterative
   startingTime = wallTime();
   callTrajOpt(drone, dataPtr, gliderPtr, modelPtr, nlpPtr, "STO");
   sol->computingTime.push_back(wallTime() - startingTime);
   cout << "STO ojb. val. = " << drone.routeCost << endl;
   cout << "STO time(s) " << sol->computingTime.back() << endl;
   sol->globalSol.push_back(drone);
   
   //call STO-NLP
   nlpPtr->setLocalSolver("ipopt");
   startingTime = wallTime();
   callTrajOpt(drone, dataPtr, gliderPtr, modelPtr, nlpPtr, "STO-NLP");
   sol->computingTime.push_back(wallTime() - startingTime);
   cout << "STO-NLP ojb. val. = " << drone.routeCost << endl;
   cout << "STO-NLP time(s) " << sol->computingTime.back() << endl;
   sol->globalSol.push_back(drone);

   // //call STO-NLP
   // nlpPtr->setLocalSolver("octeract-engine");
   // startingTime = wallTime();
   // callTrajOpt(drone, dataPtr, gliderPtr, modelPtr, nlpPtr, "STO-NLP");
   // sol->computingTime.push_back(wallTime() - startingTime);
   // cout << "STO-NLP ojb. val. = " << drone.routeCost << endl;
   // cout << "STO-NLP time(s) " << sol->computingTime.back() << endl;
   // sol->globalSol.push_back(drone);

   // //call STO-NLP
   // nlpPtr->setLocalSolver("bonmin");
   // startingTime = wallTime();
   // callTrajOpt(drone, dataPtr, gliderPtr, modelPtr, nlpPtr, "STO-NLP");
   // sol->computingTime.push_back(wallTime() - startingTime);
   // cout << "STO-NLP ojb. val. = " << drone.routeCost << endl;
   // cout << "STO-NLP time(s) " << sol->computingTime.back() << endl;
   // sol->globalSol.push_back(drone);

   // //call STO-NLP
   // nlpPtr->setLocalSolver("lgo");
   // startingTime = wallTime();
   // callTrajOpt(drone, dataPtr, gliderPtr, modelPtr, nlpPtr, "STO-NLP");
   // sol->computingTime.push_back(wallTime() - startingTime);
   // cout << "STO-NLP ojb. val. = " << drone.routeCost << endl;
   // cout << "STO-NLP time(s) " << sol->computingTime.back() << endl;
   // sol->globalSol.push_back(drone);
   
   double routeLb = getRoutelb(drone, gliderPtr);
   cout << "\n\n Lower bound on the flight time = " << routeLb << endl;

   //write solution to file
   writeSol2File(dataPtr,sol->globalSol);
   
   //write results to file
   writeRes2File(dataPtr,routeLb,sol,dataPtr->instType);

   //delete pointers
   delete sol;
   delete modelPtr;
   delete gliderPtr;
   delete dataPtr;
   return 0;
}

void callTrajOpt(glider& drone, Data* data, Dynamics* dyn, CPLEXmodel* model, AMPLmodel* nlp, string method)
{
   double const beta = 0.5;
   //Initialise aux vars
   int stateSize = dyn->stateSize;
   int controlSize = dyn->controlSize;
   const int T = data->T;

   //Initialise cplex's status variable
   IloAlgorithm::Status status = IloAlgorithm::Infeasible;

   //Initialise trajectory optimisation object
   TrajOpt* trajectory = new TrajOpt(data, dyn, model, nlp, method);

   //Compute trajectories and update real costs
   int rSize = drone.route.size();

   //Initialise glider variables
   drone.routeCost = 0;
   drone.isFeasible = false;
   drone.flightTimes = new double[rSize - 1];
   drone.stepSizes = new double[rSize - 1];
   drone.errors = new double[rSize - 1];
   for(int j = 0; j < rSize - 1; j++){
      drone.flightTimes[j] = 0;
      drone.stepSizes[j] = 0;
      drone.errors[j] = 0;
   }
   drone.normEpsilons = new double* [rSize -1];
   for(int j = 0; j < rSize - 1; j++){
      drone.normEpsilons[j] = new double[T];
   }
   drone.normTaylor1st = new double* [rSize -1];
   for(int j = 0; j < rSize - 1; j++){
      drone.normTaylor1st[j] = new double[T];
   }
   drone.relEpsilons = new double* [rSize -1];
   for(int j = 0; j < rSize - 1; j++){
      drone.relEpsilons[j] = new double[T];
   }
   
   //Initialise solution array
   drone.trajectory.resize((rSize - 1)*T,vector<double>(stateSize + controlSize));

   //Optimise trajectory for the given route
   status = trajectory->findTraj(drone.route, drone.route.size(),
         drone.trajectory, drone.flightTimes,
         drone.stepSizes, drone.errors, drone.infOrig, 
         drone.normEpsilons, drone.normTaylor1st, drone.relEpsilons);

   //Save solution info
   if(status == IloAlgorithm::Optimal){
      //Update route cost
      drone.isFeasible = true;
      for(int j = 0; j < rSize - 1; j++){
         //Computation of the real objective function value
         //arc cost
         double arcCost = drone.flightTimes[j] + drone.errors[j];
         drone.routeCost += arcCost;
      }     
   }else{
      cout << "Failed to find feasible trajectory." << endl;
      drone.routeCost = DBL_MAX;
   }

}










