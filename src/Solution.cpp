#include "Solution.h"

Solution::Solution()
{
   //set version
   makespan = true;
   //set timers
   stoTime = 0;
   stoNlpTime = 0;
}

Solution::~Solution()
{
}

bool Solution::isMakespan()
{
   return makespan;
}

double Solution::computeRoutesCost(vector<glider>& gliders, double** flightTimes)
{
   for(int i = 0; i < gliders.size(); i++){
      //reset route cost
      gliders[i].routeCost = 0;
      //compute new cost
      int size = gliders[i].route.size();
      for(int j = 0; j < size - 1; j++){
         int ii = gliders[i].route[j];
         int jj = gliders[i].route[j+1];
         double cij = flightTimes[ii][jj];
         //update routes cost
         gliders[i].routeCost += cij;
      }
   }
}

double Solution::getTotalCost(vector<glider>& gliders)
{
   double flightTime = 0;
   double totalError = 0;

   for(int i = 0; i < gliders.size(); i++){
      if(gliders[i].isFeasible){
         int rSize = gliders[i].route.size();
         for(int j = 0; j < rSize - 1; j++){
            totalError += gliders[i].errors[j];
         }
      }
   }

   //compute makespan or total flight time
   if(makespan){
      for(int i = 0; i < gliders.size(); i++){
         if(gliders[i].routeCost > flightTime + 0.001){
            flightTime = gliders[i].routeCost;
         }
      }
   }else{
      for(int i = 0; i < gliders.size(); i++){
         flightTime += gliders[i].routeCost;
      }
   }

   return flightTime + totalError;
}


int Solution::getMkspRoute(vector<glider>& gliders)
{
   double cost = 0;
   int index = -1;

   //compute total routing cost
   if(makespan){
      for(int i = 0; i < gliders.size(); i++){
         if(gliders[i].routeCost > cost + 0.001){
            cost = gliders[i].routeCost;
            index = i;
         }
      }
   }else{
      cout << "Conflicting options! Solution::getMkspRoute()" << endl;
      exit(1);
   }

   return index;
}

