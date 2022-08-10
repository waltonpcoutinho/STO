#include "Util.h"

void setSeed()
{
   // initialise random seed
   int seed = time(NULL);
   seed = 1637296301;
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

   static const char alphanum[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
   
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
   for(int i = 0; i < routeSize - 1; i++){
      totalTime += uav.flightTimes[i];
   }
   return totalTime;
}

double routeTotalError(glider uav){
   double totalError = 0;
   int routeSize = uav.route.size();
   for(int i = 0; i < routeSize - 1; i++){
      totalError += uav.errors[i];
   }
   return totalError;
}

double routeAvgStep(glider uav){
   double totalStep = 0;
   int routeSize = uav.route.size();
   for(int i = 0; i < routeSize - 1; i++){
      totalStep += uav.stepSizes[i];
   }
   return totalStep/routeSize;
}

bool checkIfEmpty(fstream& pFile)
{
   //test if file is empty by trying to read from it 
   bool test = pFile.peek() == std::ifstream::traits_type::eof();

   //the test above sets failbit itself, since it reached EOF 
   //which would disallow the rest of the functions to work
   //so we need to reset the state of the stream to goodbit
   pFile.clear();
   pFile.seekg(0, ios::beg);

   //return the result of the test
   return test;
}

void writeRes2File(Data const* dataPtr, double lb, Solution const* sol, string instType)
{
   //auxiliary vectors
   vector<double> gliderFlightTime;
   vector<double> gliderError;
   vector<double> gliderStep;

   //compute results for each glider in solution
   for(int i = 0; i < sol->globalSol.size(); i++){
      gliderFlightTime.push_back(routeTotalFlightTime(sol->globalSol[i]));
      gliderError.push_back(routeTotalError(sol->globalSol[i]));
      gliderStep.push_back(routeAvgStep(sol->globalSol[i]));
   }

   //write results to file
   string resultsPath = "Results/compResults_" + instType + ".out";
   fstream writeRes(resultsPath, ios::out | ios::in | ios::app);

   //check if file is empty and write header if so
   if(checkIfEmpty(writeRes)){
      writeRes << "Instance " << "Lower_bound ";
      //for each glider (solver or any other alternative) write the following:
      for(int i = 0; i < sol->globalSol.size();i++){
         writeRes << "flight_time_" << i << " ";
      }
      for(int i = 0; i < sol->globalSol.size();i++){
         writeRes << "error_" << i << " ";
      }
      for(int i = 0; i < sol->globalSol.size();i++){
         writeRes << "step_" << i << " ";
      }
      for(int i = 0; i < sol->globalSol.size();i++){
         writeRes << "time_" << i << " ";
      }
      writeRes << endl;
   }

   writeRes << dataPtr->fileName << " " << lb << " ";

   for(int i = 0; i < sol->globalSol.size(); i++){
      writeRes << gliderFlightTime[i] << " ";
   }
   for(int i = 0; i < sol->globalSol.size(); i++){
      writeRes << gliderError[i] << " ";
   }
   for(int i = 0; i < sol->globalSol.size(); i++){
      writeRes << gliderStep[i] << " ";
   }
   for(int i = 0; i < sol->computingTime.size(); i++){
      writeRes << sol->computingTime[i] << " ";
   }
   writeRes << endl;
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
   //append discretisation size to file name
   char aux[10];
   sprintf(aux,"_N_%d_",data->T);
   strcat(destination,aux);
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
   //print the max norm of epsilon for each route and time step
   writeSol << "max norm of epsilon" << endl;
   for(int i = 0; i < gliders.size(); i++){
      int routeSize = gliders[i].route.size();
      for(int j = 0; j < routeSize - 1; j++){
         for(int t = 0; t < data->T; t++){
            writeSol << gliders[i].normEpsilons[j][t] << " ";
         }
         writeSol << endl;
      }
   }
   //print maxnorm of 1st term of Taylor approx for each route and time step
   writeSol << "max norm of 1st Taylor term" << endl;
   for(int i = 0; i < gliders.size(); i++){
      int routeSize = gliders[i].route.size();
      for(int j = 0; j < routeSize - 1; j++){
         for(int t = 0; t < data->T; t++){
            writeSol << gliders[i].normTaylor1st[j][t] << " ";
         }
         writeSol << endl;
      }
   }
   //print the relative errors for each route and time step
   writeSol << "relative errors" << endl;
   for(int i = 0; i < gliders.size(); i++){
      int routeSize = gliders[i].route.size();
      for(int j = 0; j < routeSize - 1; j++){
         for(int t = 0; t < data->T; t++){
            writeSol << gliders[i].relEpsilons[j][t] << " ";
         }
         writeSol << endl;
      }
   }
}

















