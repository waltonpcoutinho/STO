#include "Data.h"

Data::Data(int argc, char** argv)
{
   //load input 
   if(argc < 5){
      cerr << "Not enough input arguments! " << endl;
      cerr << "Calling sequence: ./exeGTO [INSTANCE PATH] [NUMBER OF STEPS] [MAX FLEET SIZE] [INSTANCE TYPE]" << endl;
      exit(1);
   }
   if(argc > 5){//
      cerr << "Too many input arguments! " << endl;
      cerr << "Calling sequence: ./exeGTO [INSTANCE PATH] [NUMBER OF STEPS] [MAX FLEET SIZE] [INSTANCE TYPE]" << endl;
      exit(1);
   }
   //extract instance name from file   
   char* instPath = argv[1];   
   getInstName(instPath);
   
   //constants
   T = atoi(argv[2]);
   PI = 4*atan(1);   
   maxFleetSize = atoi(argv[3]);
   printf("\nDiscretisation size: %d\n",T);

   //get instance type
   instType = argv[4];

   //read instance
   ifstream readInst(instPath, ios::in);
   if(!readInst){  
      cerr << "Failed to open instance file in Instances/" << endl;
      exit(1);
   }
   
   //read depot coordinates
   depot = new double [3];
   readInst >> depot[0] >> depot[1] >> depot[2];  
   
   //read waypoints data from instance file
   noWaypoints = 0;
   readInst >> noWaypoints;
   xBar = new double [noWaypoints];
   yBar = new double [noWaypoints];
   zBar = new double [noWaypoints];
   minZBar = new double [noWaypoints];
   maxZBar = new double [noWaypoints];
   radiusBar = new double [noWaypoints];   
   for(int i = 0; i < noWaypoints; i++){
      readInst >> xBar[i];
      readInst >> yBar[i];
      readInst >> zBar[i];
      readInst >> radiusBar[i];
      readInst >> minZBar[i];
      readInst >> maxZBar[i];
   }

   //read landing sites data from instance file
   noLandPoints = 0;
   readInst >> noLandPoints;
   xTilde = new double [noLandPoints];
   yTilde = new double [noLandPoints];
   zTilde = new double [noLandPoints];
   radiusTilde = new double [noLandPoints];   
   for(int i = 0; i < noLandPoints; i++){
      readInst >> xTilde[i];
      readInst >> yTilde[i];
      readInst >> zTilde[i];
      readInst >> radiusTilde[i];
   }

   //initialise distance matrices
   //distance from depot to waypoints
   distFromDepot = new double [noWaypoints];
   //distance between waypoints
   distMatrix = new double* [noWaypoints];
   for(int i = 0; i < noWaypoints; i++){
      distMatrix[i] = new double[noWaypoints];
   }
   //distance between waypoints and landing sites
   distToLanding = new double* [noWaypoints];
   for(int i = 0; i < noWaypoints; i++){
      distToLanding[i] = new double[noLandPoints];
   }
   
   //compute dist. matrices and other constants
   computeConstants();

   //close file and end constructor
   readInst.close();
   
   //compute fleet size
   //fleetSize = floor((double) noWaypoints/2.0);
   fleetSize = maxFleetSize;
   if(fleetSize < 1){
      fleetSize = 1;
   }
   if(fleetSize >= maxFleetSize){
      fleetSize = maxFleetSize;
   }
   cout << "Fleet size: " << fleetSize << endl;
}

// destructor
Data::~Data()
{
   delete[] xBar;
   delete[] yBar;
   delete[] zBar;
   delete[] minZBar;
   delete[] maxZBar;
   delete[] radiusBar;

   delete[] xTilde;
   delete[] yTilde;
   delete[] zTilde;
   delete[] radiusTilde;
   
   for(int i = 0; i < noWaypoints; i++){
      delete[] distMatrix[i];
   }
   delete[] distMatrix;
   
   for(int i = 0; i < noLandPoints; i++){
      delete[] distToLanding[i];
   }
   delete[] distToLanding;

   xBar = NULL;
   yBar = NULL;
   zBar = NULL;
   minZBar = NULL;
   maxZBar = NULL;
   radiusBar = NULL;

   xTilde = NULL;
   yTilde = NULL;
   zTilde = NULL;
   radiusTilde = NULL;

   distMatrix = NULL;
   distToLanding = NULL;
}

void Data::computeConstants()
{
   //waypoints' coordinates
   vector<double> x(noWaypoints);
   vector<double> y(noWaypoints);
   vector<double> z(noWaypoints);
   for(int i = 0; i < noWaypoints; i++){
      x[i] = xBar[i];
      y[i] = yBar[i];
      z[i] = zBar[i];
   }  
   //calculating distance from depot to each waypoint
   for(int i = 0; i < noWaypoints; i++){
      distFromDepot[i] = sqrt((x[i]-depot[0])*(x[i]-depot[0])+
                              (y[i]-depot[1])*(y[i]-depot[1])+
                              (z[i]-depot[2])*(z[i]-depot[2]));

   }
   //calculating distances between waypoints
   for(int i = 0; i < noWaypoints; i++){
      for(int j = 0; j < noWaypoints; j++){
         distMatrix[i][j] = sqrt((x[i]-x[j])*(x[i]-x[j])+
                                 (y[i]-y[j])*(y[i]-y[j])+
                                 (z[i]-z[j])*(z[i]-z[j]));
      }
   }
   //landing sites coordinates
   vector<double> xLan(noLandPoints);
   vector<double> yLan(noLandPoints);
   vector<double> zLan(noLandPoints);
   for(int i = 0; i < noLandPoints; i++){
      xLan[i] = xTilde[i];
      yLan[i] = yTilde[i];
      zLan[i] = zTilde[i];
   }   
   //calculating distances from waypoints to landing sites
   for(int i = 0; i < noWaypoints; i++){
      for(int j = 0; j < noLandPoints; j++){
         distToLanding[i][j] = sqrt((x[i]-xLan[j])*(x[i]-xLan[j])+
                                    (y[i]-yLan[j])*(y[i]-yLan[j])+
                                    (z[i]-zLan[j])*(z[i]-zLan[j]));
      }
   }
/*   
   //print matrices
   cout << "\nEuclidean Distance Matrices: " << endl;
   cout << "Distance from depot" << endl;
   for(int i = 0; i < noWaypoints; i++){
      cout << distFromDepot[i] << endl;

   }
   cout << endl;
   cout << "Distance between waypoints" << endl;
   for(int i = 0; i < noWaypoints; i++){
      for(int j = 0; j < noWaypoints; j++){
         cout << distMatrix[i][j] << "\t";
      }
      cout << endl;
   }
   cout << "\nDistance between waypoints and landing sites" << endl;
   for(int i = 0; i < noWaypoints; i++){
      for(int j = 0; j < noLandPoints; j++){
         cout << distToLanding[i][j] << "\t";
      }
      cout << endl;
   }
*/
   //create new auxiliary vectors and insert 
   //waypoints and landing sites coordinates
   vector<double> X(x.size() + xLan.size());
   vector<double> Y(y.size() + yLan.size());
   vector<double> Z(z.size() + zLan.size());

   X.insert(X.end(),x.begin(),x.end());
   X.insert(X.end(),xLan.begin(),xLan.end());
   Y.insert(Y.end(),y.begin(),y.end());
   Y.insert(Y.end(),yLan.begin(),yLan.end());
   Z.insert(Z.end(),z.begin(),z.end());
   Z.insert(Z.end(),zLan.begin(),zLan.end());
   
   //sort vectors
   sort(X.begin(),X.end());
   sort(Y.begin(),Y.end());
   sort(Z.begin(),Z.end());

   double deltaX = X.back() - X.front();
   double deltaY = Y.back() - Y.front();
   double deltaZ = Z.back() - Z.front();

   double spaceDiagonal = sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);

   //compute x,y,z bounds
   xLb = X.front() - spaceDiagonal;
   yLb = Y.front() - spaceDiagonal;
   zLb = 0.001;
   xUb = X.back() + spaceDiagonal;
   yUb = Y.back() + spaceDiagonal;
   zUb = Z.back() + spaceDiagonal;
}

void Data::getInstName(char* instPath)
{
   string str = instPath;
   string DirectoryName;		
   string::size_type found = str.find_first_of("/");
   DirectoryName.append(instPath,found);
   string::size_type loc = str.find_last_of(".", str.size());
   string::size_type loc2 = str.find_last_of("/", str.size());
   if (loc != string::npos) {
      fileName.append(instPath, loc2+1, loc-loc2-1 );
   }
   else {
      fileName.append(instPath, loc2+1, str.size() );
   }
   cout << "Instance name: " << fileName;
}





