#include "Util.h"

double fleetSize_util = 0;
double availGlider_utils = 0;

void setSeed()
{
   // initialise random seed
   int seed = time(NULL);
   srand(seed);
   cout << "\n\nRandom seed: " << seed << endl;
}

bool randomBool()
{
   static auto gen = bind(uniform_int_distribution<>(0,1),default_random_engine());
   return gen();
}

double fRand(double fMin, double fMax)
{
   double f = (double)rand() / RAND_MAX;
   return fMin + f * (fMax - fMin);
}

double cpuTime(){
    return (double)clock() / CLOCKS_PER_SEC;
}

double wallTime(){
   struct timeval time;
   if (gettimeofday(&time,NULL)){
      //  Handle error
      return 0;
   }
   return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

void randomSeq(int* route, int seqSize)
{
   const int AMOUNT = seqSize; //amount of random numbers that need to be generated
   const int MAX = seqSize; //maximum value (of course, this must be at least the same as AMOUNT;
   int* value = route;

   //generate random numbers:
   for(int i = 0; i < AMOUNT; i++){
      bool check; //variable to check or number is already used
      int n; //variable to store the number in
      do{
         n = rand()%MAX;
         //check or number is already used:
         check = true;
         for(int j = 0; j < i; j++){
            if(n == value[j]){
               check = false; //set check to false
               break;
            }
         }
      }while(!check); //loop until new, unique number is found
      value[i] = n; //store the generated number in the array
      if(i == 0){
         value[i] = 0;
      }
   }
   //at this point in the program we have an array value[] with a serie of unique random numbers
}

double l2_norm(const double x, const double y, const double z) 
{
   return sqrt(x*x + y*y + z*z);
}
string genRandStr(const int len){

   char key[len];

   static const char alphanum[] =
      "0123456789"
      "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      "abcdefghijklmnopqrstuvwxyz";
   
   default_random_engine generator;
   uniform_int_distribution<int> uniform(0,sizeof(alphanum)-1);

   for (int i = 0; i < len; i++) {
      key[i] = alphanum[uniform(generator)];
   }

   return key;
}

template<typename Type> void printMatrix(Type mat, int M, int N)
{
    cout<<"\nPrinting Matrix : " << endl;
    for(int i = 0 ; i < M; ++i){
        for(int j = 0 ; j < N; ++j){
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void printVec(double const* vec, int length, string text)
{
   cout << text << "(";
   for(int i = 0; i < length -1; i++){
      cout << vec[i] << ", ";
   }
   cout << vec[length-1] << ") " << endl;
}

void printVec(int const* vec, int length, string text)
{
   cout << text << "(";
   for(int i = 0; i < length -1; i++){
      cout << vec[i] << ", ";
   }
   cout << vec[length-1] << ") " << endl;
}

void printVec(vector<int> vec, string text)
{
   cout << text << "(";
   for(int i = 0; i < vec.size()-1; i++){
      cout << vec[i] << ", ";
   }
   cout << vec[vec.size()-1] << ") " << endl;
}

void printRoutes(vector<glider> gliders)
{
   for(int i = 0; i < gliders.size(); i++){
      cout << "glider[" << i << "] cost= " << gliders[i].routeCost << " | ";
      printVec(gliders[i].route, "route: ");
      std::cout.flush();
   }
}

double routeTotalFlightTime(glider uav){
   double totalTime = 0;
   int routeSize = uav.route.size();
   for(int i = 0; i < routeSize ; i++){
      totalTime += uav.flightTimes[i];
   }
   return totalTime;
}

double routeTotalError(glider uav){
   double totalError = 0;
   int routeSize = uav.route.size();
   for(int i = 0; i < routeSize ; i++){
      totalError += uav.errors[i];
   }
   return totalError;
}

double routeAvgStep(glider uav){
   double totalStep = 0;
   int routeSize = uav.route.size();
   for(int i = 0; i < routeSize ; i++){
      totalStep += uav.stepSizes[i];
   }
   return totalStep/routeSize;
}

void writeRes2File(Data const* dataPtr, Solution const* sol)
{
   vector<glider> glob = sol->globalSol;

   double stoFlTime = routeTotalFlightTime(glob[0]);
   double stoNlpFlTime = routeTotalFlightTime(glob[1]);

   double stoError = routeTotalError(glob[0]);
   double stoNlpError = routeTotalError(glob[1]);

   double stoStep = routeAvgStep(glob[0]);
   double stoNlpStep = routeAvgStep(glob[1]);

   double stoTime = sol->stoTime;
   double stoNlpTime = sol->stoNlpTime;

   //write results to file
   ofstream writeRes("Results/compResults_" + genRandStr(12) + ".out", ios::app);
   writeRes << dataPtr->fileName << " "
            << stoFlTime << " "
            << stoNlpFlTime << " "
            << stoError << " " 
            << stoNlpError << " " 
            << stoStep << " " 
            << stoNlpStep << " " 
            << stoTime << " "
            << stoNlpTime << endl;
   writeRes.close();
}

void writeInfList(Data const* dataPtr)
{
   //write results to file
   ofstream file("Results/infeasibleInstances.out", ios::app);
   file << dataPtr->fileName << endl;
   file.close();
}

void writeSol2File(Data const* data,vector<glider> gliders)
{
   cout << "\nPrinting results to file." << endl;
   cout << "Instance: " << data->fileName << endl;
   //number of available gliders 
   cout << "Number of available gliders = " << data->maxFleetSize << endl;
   cout << " Deleting empty routes." << endl;
   vector<glider> temp;
   //printRoutes(gliders);
   //delete empty routes
   for(int i = 0; i < gliders.size(); i++){
      if(gliders[i].route.size() > 2){
         temp.push_back(gliders[i]);
      }
   }
   //copy non-empty routes back 
   gliders.clear();
   gliders = temp;
   //printRoutes(gliders);
   cout << " #Gliders = " << gliders.size() << endl;
   cout << "Done!" << endl;
   
   //get file destination folder
   char destination[100];
   strcpy(destination,"Results/Solutions/solution_");
   strcat(destination,data->fileName.c_str());
   strcat(destination,".out");
   
   //delete file with solutions iff there is one
   if (FILE *file = fopen(destination, "r")) {
         if( remove(destination) != 0 ){
            perror("Error deleting file");
         }
      fclose(file);
   }

   ofstream writeSol(destination, ios::out);

   //write instance name
   writeSol << data->fileName << endl;
   //write discretisation size
   writeSol << data->T << endl;
   //write number of routes/gliders
   writeSol << data->maxFleetSize << " " << gliders.size() << endl;
   //write route sizes
   for(int i = 0; i < gliders.size(); i++){
      writeSol << gliders[i].route.size() << " ";
   }
   writeSol << endl;
   //print routes
   for(int i = 0; i < gliders.size(); i++){
      int routeSize = gliders[i].route.size();
      for(int j = 0; j < routeSize; j++){
         writeSol << gliders[i].route[j] << " ";
      }
      writeSol << endl;
   }
   //print flight times of each leg and each glider
   for(int i = 0; i < gliders.size(); i++){
      int routeSize = gliders[i].route.size();
      for(int j = 0; j < routeSize - 1; j++){
         writeSol << gliders[i].flightTimes[j] << " ";
      }
      writeSol << endl;
   }
   //print step sizes of each leg and each glider
   for(int i = 0; i < gliders.size(); i++){
      int routeSize = gliders[i].route.size();
      for(int j = 0; j < routeSize - 1; j++){
         writeSol << gliders[i].stepSizes[j] << " ";
      }
      writeSol << endl;
   }
   //print errors of each leg and each glider
   for(int i = 0; i < gliders.size(); i++){
      int routeSize = gliders[i].route.size();
      for(int j = 0; j < routeSize - 1; j++){
         writeSol << gliders[i].errors[j] << " ";
      }
      writeSol << endl;
   }
   //print trajectories of each glider
   for(int i = 0; i < gliders.size(); i++){
      int routeSize = gliders[i].route.size();
      int nLines = (routeSize-1)*data->T;
      for(int j = 0; j < nLines; j++){
         int nCols =  gliders[i].trajectory[j].size();
         for(int k = 0; k < nCols; k++){
            writeSol << gliders[i].trajectory[j][k] << " ";
         }
         writeSol << endl;
      }
   }
}

















