#include "Dynamics.h"

using namespace std;

Dynamics::Dynamics(Data* data)
:dataPtr(data)
{
   //constants
   PI = 4*atan(1);

   //safety distance in metres
   safetyDist = 5;

   //load system, glider and wind data
   const char* gliderWindPath = "Data/GliderData.dat";
   ifstream readGliderDat(gliderWindPath, ios::in);

   // Load data from file
   if(!readGliderDat){  
      cerr << "Failed to open data files in Data/" << endl;
      exit(1);
   }
   string tempString;   
   //load glider and wind data
   readGliderDat >> rho;
   readGliderDat >> CD0;
   readGliderDat >> K;
   readGliderDat >> gravity;
   readGliderDat >> mass;
   readGliderDat >> S;
   readGliderDat >> alpha;
   readGliderDat >> beta;
   
   //read bounds from file
   const char* boundsPath = "Data/boundsOnVars.dat";   
   ifstream readBounds(boundsPath, ios::in);
   if(!readBounds){  
      cerr << "Failed to open data file in Data/" << endl;
      exit(1);
   }
   //bounds on variables
   readBounds >> vLb >> vUb;
   readBounds >> gammaLb >> gammaUb;
   readBounds >> phiLb >> phiUb;
   readBounds >> CLlb >> CLub;
   readBounds >> muLb >> muUb;
   //bounds on derivatives
   readBounds >> xDot_ub;
   readBounds >> yDot_ub;
   readBounds >> hDot_ub;
   readBounds >> vDot_ub;
   readBounds >> gammaDot_ub;
   readBounds >> phiDot_ub;
   readBounds >> CLDot_ub;
   readBounds >> muDot_ub;

   double stateDotUb[] = {xDot_ub, yDot_ub, hDot_ub, vDot_ub, gammaDot_ub, phiDot_ub};
   double ctrlDotUb[] = {CLDot_ub, muDot_ub};

   yDotUbMaxNorm = *max_element(stateDotUb, stateDotUb + 6);
   uDotUbMaxNorm = *max_element(ctrlDotUb, ctrlDotUb + 2);

   //compute A, B and the stead state conditions
   //for all pairs of vertices
   computeDynamics();

   //compute flight time matrices
   computeFlightTimes();
   
   //compute flight distance matrices
   computeFlightDist();

   //erase file readers
   readBounds.close();
   readGliderDat.close();
}

// destrutor
Dynamics::~Dynamics()
{
   //delete all data from the legs
   for(int i = 0; i < legArrayNumRows; i++){
      for(int j = 0; j < legArrayNumCols; j++){
         if(!(j == 0 || i == j || (i == 0 && j > dataPtr->noWaypoints))){
            delete[] legDyn[i][j].Yeq;
            delete[] legDyn[i][j].Ueq;
            for(int k = 0; k < stateSize; k++){
               delete[] legDyn[i][j].stateMatrix[k];
            }
            delete[] legDyn[i][j].stateMatrix;
            for(int k = 0; k < stateSize; k++){
               delete[] legDyn[i][j].controlMatrix[k];
            }
            delete[] legDyn[i][j].controlMatrix;
         }
      }
   }
   //delete legs themselves
   for(int i = 0; i < legArrayNumRows; i++){
      delete[] legDyn[i];
   }
   delete[] legDyn;
   //delete identity matrix
   for(int i = 0; i < stateSize; i++){
      delete[] I[i];
   }
   delete[] I;

   legDyn = NULL;
   I = NULL;
}

void Dynamics::gliderODE_ADOL(dynamics leg)
{
   //active independent vars
   adouble z,v,gamma,phi;
   adouble CL, mu;
   //active dependent vars
   adouble CD, L, D, U, dU;
   adouble* dY = new adouble [stateSize];
   //passive vars
   double* pdY = new double [stateSize];
   
   //tape array and/or tape file specifier
   tag = 0;                        
   //start tracing
   trace_on(tag);                          
   
   //state vars  
   z <<= leg.get_Yeq(2);
   v <<= leg.get_Yeq(3);
   gamma <<= leg.get_Yeq(4);
   phi <<= leg.get_Yeq(5);
   //control vars
   CL <<= leg.get_Ueq(0);
   mu <<= leg.get_Ueq(1);

   //dependent vars
   CD = CD0 + K*CL*CL;
   L = 0.5*rho*S*CL*v*v;
   D = 0.5*rho*S*CD*v*v;
   U = beta*z;
   dU = beta*v*sin(gamma);

   //system of ODEs
   dY[0] = v*cos(gamma)*sin(phi) + U; 
   dY[1] = v*cos(gamma)*cos(phi);
   dY[2] = v*sin(gamma);
   dY[3] = -D/mass -gravity*sin(gamma) -dU*cos(gamma)*sin(phi);
   dY[4] = L*cos(mu)/(mass*v) -gravity*cos(gamma)/v + dU*sin(gamma)*sin(phi)/v;
   dY[5] = L*sin(mu)/(mass*v*cos(gamma)) - dU*cos(phi)/(v*cos(gamma));

   //mark and pass dependents
   for(int i = 0; i < stateSize; i++){
      dY[i] >>= pdY[i]; 
   }

   delete[] pdY;
   delete[] dY;   
   //finish tracing
   trace_off();
}


void Dynamics::ODE_Jac(dynamics* leg){

   //Allocate memory for the Jacobian
   double** J = new double* [stateSize];
   for(int i = 0; i < stateSize; i++){
      J[i] = new double [stateSize];
   }

   //# of relevant variables for the jacobian
   int numSteady = stateSize - 2 + controlSize;

   //steady-state vector 
   //ignoring x and y columns
   double steady[] = {leg->get_Yeq(2),leg->get_Yeq(3),leg->get_Yeq(4),leg->get_Yeq(5),
                      leg->get_Ueq(0),leg->get_Ueq(1)};

   //compute jacobian
   int status = jacobian(tag,numSteady,numSteady,steady,J);

   //pass the values from the jacobian to the system matrices
   for(int i = 0; i < stateSize; i++){
      for(int j = 0; j < stateSize - 2; j++){
         leg->stateMatrix[i][j+2] = J[i][j];
      }
      for(int j = 4; j < 4 + controlSize; j++){
         leg->controlMatrix[i][j-4] = J[i][j];
      }
   }

#if DEBUG > 0
   cout << "Steady flight values" << endl;
   for(int i = 0; i < 6; i++){
      cout << steady[i] << " ";
   }
   cout << endl;
   cout << "State matrix: " << endl;
   for(int i = 0; i < stateSize; i++){
      for(int j = 0; j < stateSize - 2; j++){
         cout << setw(10) << leg->stateMatrix[i][j+2] << "\t";
      }
      cout << endl;
   }
   cout << "Control matrix: " << endl;
   for(int i = 0; i < stateSize; i++){
      for(int j = 4; j < 4 + controlSize; j++){
         cout << setw(10) << leg->controlMatrix[i][j-4] << "\t";
      }
      cout << endl;
   }
   PAUSE;
#endif

   //delete auxiliary jacobian matrix
   for(int i = 0; i < numSteady; i++){
      delete[] J[i];
   }
   delete[] J;
}

void Dynamics::steadyState(dynamics* leg)
{
   //steady state variables
   double CL_0  = 0;
   double v_0 = 0;
   double gamma_0 = 0;
   double phi_0 = 0;
   double mu_0 = 0;
   
   double posi[] = {leg->Yo[0],leg->Yo[1],leg->Yo[2]};
   double posf[] = {leg->Yf[0],leg->Yf[1],leg->Yf[2]};

   double heightDiff = posf[2] - posi[2];
   if(heightDiff > 10){
      leg->status = false;//Infeasible
      //printf("Target %d cannot be reached from origin %d (%d)!\n",leg->dest,leg->orig,leg->status);
   }else{
      leg->status = true;//feasible
      if(heightDiff >= 0 - epsl && heightDiff <= 10){
         CL_0 = 0.37;
         v_0 = 12.48;
         gamma_0 = -0.02;
         phi_0 = -1.57;
         mu_0 = 0;
      }else{
         //heightDiff is negative
         //set steady-descent equilibrium points
         //compute gamma and check intervals
         double length = sqrt(pow((posf[0] - posi[0]),2) + pow((posf[1] - posi[1]),2));
         //deltah is always negative, therefore gamma is always
         //negative,i.e., V always points down 
         gamma_0 = -atan(heightDiff/length);
         //max range compute from z to the ground
         if(gamma_0 <= gammaLb){
            gamma_0 = gammaLb;
         }
         if(gamma_0 >= gammaUb){
            gamma_0 = gammaUb;
         }
         //compute CL and check intervals
         double b = -1/cos(gamma_0) + 1 + sqrt(2*CD0*K);
         double Delta = b*b - 4*K*CD0;
         //ignore the complex part of Delta
         double sqrtDelta = Delta > 0 ? sqrt(Delta) : 0;
         double CL1 = (-b + sqrtDelta)/(2*K);
         double CL2 = (-b - sqrtDelta)/(2*K);
         if(abs(CL1) <= abs(CL2)){
            CL_0 = CL1;
         }else{
            CL_0 = CL2;
         }
         if(CL_0 <= CLlb){
            CL_0 = CLlb;
         }
         if(CL_0 >= CLub){
            CL_0 = CLub;
         }
         //compute v and check intervals
         double CD = CD0 + K*(CL_0*CL_0);
         double exp = double(3.0)/double(4.0);
         //v_0 = sqrt((2*mass*gravity*CD)/(rho*S*pow((CL_0*CL_0 + CD*CD),exp)));
         v_0 = CL_0 > 0 ? sqrt((2*mass*gravity*cos(gamma_0))/(CL_0*rho*S)) : 0;
         if(v_0 <= vLb){
            v_0 = vLb;
         }
         if(v_0 >= vUb){
            v_0 = vUb;
         }
      }
      leg->Yeq[3] = v_0;
      leg->Yeq[4] = gamma_0;
      leg->Yeq[5] = phi_0;
      leg->Ueq[0] = CL_0;
      leg->Ueq[1] = mu_0;
   }
}

void Dynamics::computeDynamics()
{
   //allocate memory for the array of data structures
   //each data struct contains A,B,Yeq,Ueq for
   //each pair of points in the instance
   int nWaypoints = dataPtr->noWaypoints;
   int nLandPoints = dataPtr->noLandPoints;
   legArrayNumRows = 1 + nWaypoints;
   legArrayNumCols = 1 + nWaypoints + nLandPoints;
   legDyn = new dynamics* [legArrayNumRows];
   for(int i = 0; i < legArrayNumRows; i++){
      legDyn[i] = new dynamics [legArrayNumCols];
   }

   //fill matrix
   for(int i = 0; i < legArrayNumRows; i++){
      for(int j = 0; j < legArrayNumCols; j++){
         if(j == 0 || i == j || (i == 0 && j > nWaypoints)){
            //cout << "x" << " ";
            legDyn[i][j].orig = -1;
            legDyn[i][j].dest = -1;
            legDyn[i][j].status = false;
            legDyn[i][j].deltaTa = 0;
            legDyn[i][j].deltaTb = 0;
            legDyn[i][j].stateMatrix = NULL;
            legDyn[i][j].stateMatrix = NULL;
            legDyn[i][j].Yeq = NULL;
            legDyn[i][j].Ueq = NULL;
            legDyn[i][j].Yo = NULL;
            legDyn[i][j].Uo = NULL;
            legDyn[i][j].Yf = NULL;
            legDyn[i][j].Uf = NULL;
         }else{
            //legDyn[i][j].orig = i - 1;
            //legDyn[i][j].dest = j - 1;
            legDyn[i][j].orig = i;
            legDyn[i][j].dest = j;
            legDyn[i][j].isLanding = false;
            if(j > nWaypoints){
               legDyn[i][j].isLanding = true;
            }
            //cout << legDyn[i][j].orig << "," << legDyn[i][j].dest << " ";

            legDyn[i][j].status = false;
            //allocate memory to system matrices
            legDyn[i][j].stateMatrix = new double* [stateSize];
            for(int k = 0; k < stateSize; k++){
               legDyn[i][j].stateMatrix[k] = new double [stateSize];
            }
            legDyn[i][j].controlMatrix = new double* [stateSize];
            for(int k = 0; k < stateSize; k++){
               legDyn[i][j].controlMatrix[k] = new double [controlSize];
            }
            //initialise state and control matrices with zeros
            for(int r = 0; r < stateSize; r++){
               for(int c = 0; c < stateSize; c++){
                  legDyn[i][j].stateMatrix[r][c] = 0;
               }
               for(int c = 0; c < controlSize; c++){
                  legDyn[i][j].controlMatrix[r][c] = 0;
               }
            }
            //allocate space for eq. vectors and fill with 0 
            //Y\setminus{x,y,z}
            legDyn[i][j].Yeq = new double [stateSize];
            legDyn[i][j].Ueq = new double [controlSize];
            for(int k = 0; k < stateSize; k++){
               legDyn[i][j].Yeq[k] = 0;
            }
            for(int k = 0; k < controlSize; k++){
               legDyn[i][j].Ueq[k] = 0;
            }
            //allocate space for initial and final states
            legDyn[i][j].Yo = new double [stateSize];
            legDyn[i][j].Uo = new double [controlSize];
            if(i == 0){
                  legDyn[i][j].Yo[0] = dataPtr->depot[0];
                  legDyn[i][j].Yo[1] = dataPtr->depot[1];
                  legDyn[i][j].Yo[2] = dataPtr->depot[2];
            }else{
                  legDyn[i][j].Yo[0] = dataPtr->xBar[i-1];
                  legDyn[i][j].Yo[1] = dataPtr->yBar[i-1];
                  legDyn[i][j].Yo[2] = dataPtr->zBar[i-1];
            }
            for(int k = 0; k < controlSize; k++){
               legDyn[i][j].Uo[k] = 0;
            }
            legDyn[i][j].Yf = new double [stateSize];
            legDyn[i][j].Uf = new double [controlSize];
            if(j > nWaypoints){
                  legDyn[i][j].Yf[0] = dataPtr->xTilde[j-1-nWaypoints];
                  legDyn[i][j].Yf[1] = dataPtr->yTilde[j-1-nWaypoints];
                  legDyn[i][j].Yf[2] = dataPtr->zTilde[j-1-nWaypoints];
                  legDyn[i][j].fRadius = dataPtr->radiusTilde[j-1-nWaypoints];
                  legDyn[i][j].fzMin = 0.0;
                  legDyn[i][j].fzMax = 0.0;
            }else{
                  legDyn[i][j].Yf[0] = dataPtr->xBar[j-1];
                  legDyn[i][j].Yf[1] = dataPtr->yBar[j-1];
                  legDyn[i][j].Yf[2] = dataPtr->zBar[j-1];
                  legDyn[i][j].fRadius = dataPtr->radiusBar[j-1];
                  legDyn[i][j].fzMin = dataPtr->minZBar[j-1];
                  legDyn[i][j].fzMax = dataPtr->maxZBar[j-1];
            }
            for(int k = 0; k < controlSize; k++){
               legDyn[i][j].Uf[k] = 0;
            }
            legDyn[i][j].vo = new double [3];
            for(int k = 0; k < 3; k++){
               legDyn[i][j].vo[k] = 0;
            }
            //compute lower bound and estimate upper bound on flight time
            double deltaX = legDyn[i][j].Yf[0] - legDyn[i][j].Yo[0];
            double deltaY = legDyn[i][j].Yf[1] - legDyn[i][j].Yo[1];
            double deltaZ = legDyn[i][j].Yf[2] - legDyn[i][j].Yo[2];
            double deltaS = sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);
            legDyn[i][j].deltaTa = deltaS/vUb;
            legDyn[i][j].deltaTb = deltaS/vLb;
            
            //compute the respective steady state values
            //for the given leg
            steadyState(&legDyn[i][j]);
            gliderODE_ADOL(legDyn[i][j]);
            ODE_Jac(&legDyn[i][j]);

            //compute max norms of A and B
            double maxA = 0;
            for(int s1 = 0; s1 < stateSize; s1++){
               for(int s2 = 0; s2 < stateSize; s2++){
                  double aux = abs(legDyn[i][j].stateMatrix[s1][s2]);
                  if(aux > maxA){
                     maxA = aux;
                  }
               }
            }
            legDyn[i][j].stateMatrixMaxNorm = maxA;
            
            double maxB = 0;
            for(int s = 0; s < stateSize; s++){
               for(int c = 0; c < controlSize; c++){
                  double aux = abs(legDyn[i][j].controlMatrix[s][c]);
                  if(aux > maxB){
                     maxB = aux;
                  }
               }
            }
            legDyn[i][j].ctrlMatrixMaxNorm = maxB;
         }
      }
   }

   //declare identity matrix
   I = new double* [stateSize];
   for(int i = 0; i < stateSize; i++){
      I[i] = new double [stateSize];
   }
   //fill Identity matrix
   for(int i = 0; i < stateSize; i++){
      for(int j = 0; j < stateSize; j++){
         if(i == j){
            I[i][j] = 1;
         }
         else{
            I[i][j] = 0;
         }
      }
   }
}

void Dynamics::computeFlightTimes()
{
   int nWaypoints = dataPtr->noWaypoints;
   int nLandPoints = dataPtr->noLandPoints;
   int numRows = legArrayNumRows;
   int numCols = legArrayNumCols;
   
   flightTimes = new double* [numRows];
   for(int i = 0; i < numRows; i++){
      flightTimes[i] = new double [numCols];
   }
   
   //fill in flight time matrix
   for(int i = 0; i < numRows; i++){
      for(int j = 0; j < numCols; j++){
         if(i == 0 && j > nWaypoints){
            flightTimes[i][j] = 0;
         }
         if(j == 0 || i == j){
            flightTimes[i][j] = DBL_MAX;
         }else{
            if(i == 0){
               if(j <= nWaypoints){
                  int jj = j - 1;
                  flightTimes[i][j] = dataPtr->distFromDepot[jj]/legDyn[i][j].Yeq[3];
               }
            }else if(j > nWaypoints){
               int ii = i - 1;
               int jj = j - (nWaypoints + 1);
               flightTimes[i][j] = dataPtr->distToLanding[ii][jj]/legDyn[i][j].Yeq[3];
            }else{
               int ii = i - 1;
               int jj = j - 1;
               flightTimes[i][j] = dataPtr->distMatrix[ii][jj]/legDyn[i][j].Yeq[3];
            }
         }
      }
   }
}

void Dynamics::randomFlightTimes()
{
   int nWaypoints = dataPtr->noWaypoints;
   int nLandPoints = dataPtr->noLandPoints;
   int numRows = legArrayNumRows;
   int numCols = legArrayNumCols;
   
   flightTimes = new double* [numRows];
   for(int i = 0; i < numRows; i++){
      flightTimes[i] = new double [numCols];
   }
   
   //fill in flight time matrix
   for(int i = 0; i < numRows; i++){
      for(int j = 0; j < numCols; j++){
         if(i == 0 && j > nWaypoints){
            flightTimes[i][j] = 0;
         }
         if(j == 0 || i == j){
            flightTimes[i][j] = DBL_MAX;
         }else{
            if(i == 0){
               if(j <= nWaypoints){
                  int jj = j - 1;
                  flightTimes[i][j] = dataPtr->distFromDepot[jj]/fRand(0.1,vUb);
               }
            }else if(j > nWaypoints){
               int ii = i - 1;
               int jj = j - (nWaypoints + 1);
               flightTimes[i][j] = dataPtr->distToLanding[ii][jj]/fRand(0.1,vUb);
            }else{
               int ii = i - 1;
               int jj = j - 1;
               flightTimes[i][j] = dataPtr->distMatrix[ii][jj]/fRand(0.1,vUb);
            }
         }
      }
   }
}

void Dynamics::computeFlightDist()
{
   int nWaypoints = dataPtr->noWaypoints;
   int nLandPoints = dataPtr->noLandPoints;
   int numRows = legArrayNumRows;
   int numCols = legArrayNumCols;
   
   flightDist = new double* [numRows];
   for(int i = 0; i < numRows; i++){
      flightDist[i] = new double [numCols];
   }

   for(int i = 0; i < numRows; i++){
      for(int j = 0; j < numCols; j++){
         flightDist[i][j] = 0;
      }
   }
   
   //fill in flight time matrix
   for(int i = 0; i < numRows; i++){
      for(int j = 0; j < numCols; j++){
         if(i == 0 && j > nWaypoints){
            flightDist[i][j] = 0;
         }
         if(j == 0 || i == j){
            flightDist[i][j] = DBL_MAX;
         }else{
            if(i == 0){
               if(j <= nWaypoints){
                  int jj = j - 1;
                  flightDist[i][j] = dataPtr->distFromDepot[jj];
               }
            }else if(j > nWaypoints){
               int ii = i - 1;
               int jj = j - (nWaypoints + 1);
               flightDist[i][j] = dataPtr->distToLanding[ii][jj];
            }else{
               int ii = i - 1;
               int jj = j - 1;
               flightDist[i][j] = dataPtr->distMatrix[ii][jj];
            }
         }
      }
   } 
}

void Dynamics::printFlightMatrices()
{
   int numRows = legArrayNumRows;
   int numCols = legArrayNumCols;
   int nWaypoints = dataPtr->noWaypoints;

   //print log for debbuging
   cout << "\nEquilibrium velocities" << endl;
   for(int i = 0; i < numRows; i++){
      for(int j = 0; j < numCols; j++){
         if(j == 0 || i == j || (i == 0 && j > nWaypoints)){
            cout << "-\t";
         }else{
            cout << legDyn[i][j].Yeq[3] << "\t";
         }
      }
      cout << endl;
   }

   cout << "\nTime matrix" << endl;
   for(int i = 0; i < numRows; i++){
      for(int j = 0; j < numCols; j++){
         if(flightTimes[i][j] > 99999){
            cout << "-\t";
         }else{
            printf("%.2f\t",flightTimes[i][j]);
         }
      }
      cout << endl;
   }
 
   cout << "\nFlight dist matrix" << endl;
   for(int i = 0; i < numRows; i++){
      for(int j = 0; j < numCols; j++){
         if(flightDist[i][j] > 99999){
            cout << "-\t";
         }else{
            printf("%.2f\t",flightDist[i][j]);
         }
      }
      cout << endl;
   }
   cout.flush();
   //PAUSE;
}


Dynamics::dynamics Dynamics::getLeg(int i, int j)
{
   if(i >= legArrayNumRows || j >= legArrayNumCols){
      cout << "i or j out of bounds" << endl;
      cout << "EXITTING PROGRAMME!" << endl;
      exit(1);
   }
   if(i == j){
      cout << "Please, provide i != j!" << endl;
      exit(1);
   }else{
      return legDyn[i][j];
   }
}

double Dynamics::windSpeed(double height)
{
   return beta*height;
}

void Dynamics::analyticJac(dynamics leg)
{
   cout << "\nSteady-state values Analytic Jac: " << endl;
   cout << "v_eq = " << leg.Yeq[3] << endl;
   cout << "gamma_eq = " << leg.Yeq[4] << endl;
   cout << "phi_eq = " << leg.Yeq[5] << endl;
   cout << "Cl_eq = " << leg.Ueq[0] << endl;
   cout << "mu_eq = " << leg.Ueq[1] << endl;
   cout << "\n" << endl;

   double v_0 = leg.Yeq[3];
   double gamma_0 = leg.Yeq[4];
   double phi_0 = leg.Yeq[5];
   double CL_0 = leg.Ueq[0];
   double mu_0 = leg.Ueq[1];

   vector <vector<double>> A(stateSize,vector<double>(stateSize));
   A =//A[0][..]
      {{0, 0, beta, 
          cos(gamma_0)*cos(phi_0), (-v_0)*sin(gamma_0)*cos(phi_0), 
          -v_0*cos(gamma_0)*sin(phi_0)},
      //A[1][..]
      {0, 0, 0, 
         cos(gamma_0)*sin(phi_0), -v_0*sin(phi_0)*sin(gamma_0), 
         (v_0)*cos(gamma_0)*cos(phi_0)}, 
      //A[2][..]
      {0, 0, 0, 
         -sin(gamma_0), -v_0*cos(gamma_0), 0},
      //A[3][..]
      {0, 0, 0, 
         -(((CD0 + CL_0*CL_0*K)*rho*S*v_0)/(double)mass) 
            -beta*cos(gamma_0)*sin(gamma_0)*cos(phi_0), 
         (-gravity)*cos(gamma_0) - beta*v_0*cos(gamma_0)*cos(gamma_0)*cos(phi_0) 
            + beta*v_0*sin(gamma_0)*sin(gamma_0)*cos(phi_0), 
         (-beta)*v_0*cos(gamma_0)*sin(phi_0)*sin(gamma_0)},
      //A[4][..]
      {0, 0, 0, 
         (gravity*cos(gamma_0))/(double)(v_0*v_0) + (CL_0*rho*S*cos(mu_0))/(double)(2*mass), 
         (gravity*sin(gamma_0))/(double)v_0 + 2*beta*cos(gamma_0)*sin(gamma_0)*cos(phi_0), 
         (beta*sin(phi_0)*sin(gamma_0)*sin(gamma_0))},
      //A[5][..]
      {0, 0, 0, 
         ((CL_0*rho*S*(1/(double)cos(gamma_0))*sin(mu_0))/(double)(2*mass)), 
         (-beta)*sin(phi_0)*(1/(double)cos(gamma_0))*(1/(double)cos(gamma_0)) 
            -(CL_0*rho*S*v_0*(1/(double)cos(gamma_0))*sin(mu_0)*tan(gamma_0))/(double)(2*mass),
         (beta*cos(phi_0)*tan(gamma_0))}};

   //define derivatives of control matrix
   vector <vector<double>> B(stateSize,vector<double>(controlSize));
   B = {{0, 0},
        {0, 0},
        {0, 0},
        //B[3][..]
        {-((CL_0*K*rho*S*v_0*v_0)/(double)mass),0}, 
        //B[4][..]
        {(rho*S*v_0*cos(mu_0))/(double)(2*mass), 
           -((CL_0*rho*S*v_0*sin(mu_0))/(double)(2*mass))}, 
        //B[5][..]
        {((rho*S*v_0*(1/(double)cos(gamma_0))*sin(mu_0))/(double)(2*mass)), 
           ((CL_0*rho*S*v_0*cos(mu_0)*(1/(double)cos(gamma_0)))/(double)(2*mass))}};

   ios::fmtflags oldSettings = cout.flags();
   cout << fixed << setiosflags(ios::showpoint) << setprecision(12);
   cout << "\nAnalytic Jac: " << endl;
   for(int i = 0; i < stateSize; i++){
      for(int j = 2; j < stateSize; j++){
         cout << A[i][j] << " ";
      }
      for(int j = 0; j < controlSize; j++){
         cout << B[i][j] << " ";
      }
      cout << endl;
   }
   cout.flags(oldSettings);
}

//#################################################################################
//#################### Nested class declarations #################################

//nested class's get functions
double Dynamics::dynamics::get_Yeq(int i) const
{
   if(i >= 0 && i <= Dynamics::stateSize){
      return dynamics::Yeq[i];
   }else{
      cout << "Index out of bounds!\n EXITING PROGRAM!" << endl;
      exit(1);
   }
}

double Dynamics::dynamics::get_Ueq(int i) const
{
   if(i >= 0 && i <= Dynamics::stateSize){
      return dynamics::Ueq[i];
   }else{
      cout << "Index out of bounds!\n EXITING PROGRAM!" << endl;
      exit(1);
   }
}

double Dynamics::dynamics::get_A(int i, int j) const
{
   if(i >= 0 && i <= Dynamics::stateSize 
         && j >= 0 && j <= Dynamics::stateSize){
      return dynamics::stateMatrix[i][j];
   }else{
      cout << "Indices out of bounds!\n EXITING PROGRAM!" << endl;
      exit(1);
   }

}

double Dynamics::dynamics::get_B(int i, int j) const
{
   if(i >= 0 && i <= Dynamics::stateSize 
         && j >= 0 && j <= Dynamics::stateSize){
      return dynamics::controlMatrix[i][j];
   }else{
      cout << "Indices out of bounds!\n EXITING PROGRAM!" << endl;
      exit(1);
   }
}

bool Dynamics::dynamics::getStatus() const
{
   return dynamics::status;
}

bool Dynamics::checkRouteStatus(vector<int>& route, int routeSize)
{
   bool feasible = false;
   for(int i = 0; i < routeSize-1; i++){
      dynamics leg = getLeg(route[i],route[i+1]);
      feasible = leg.status;
      if(!feasible){
         cout << "Arc (" << route[i] <<","<< route[i+1] <<") is forbidden!" << endl;
         return false;
      }
   }
   return true;

}

double Dynamics::dynamics::getMaxNormA() const
{
   return stateMatrixMaxNorm;
}

double Dynamics::dynamics::getMaxNormB() const
{
   return ctrlMatrixMaxNorm;
}







