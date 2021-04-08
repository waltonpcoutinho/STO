#include "CPLEXmodel.h"

ILOSTLBEGIN // import namespace std

CPLEXmodel::CPLEXmodel(Data* data, Dynamics* glider, int T)
   :dataPtr(data), gliderPtr(glider)
{
   this->T = dataPtr->T;
}	
CPLEXmodel::~CPLEXmodel()
{
   env.end();
}

double CPLEXmodel::errorNorm(IloArray<IloNumArray> epst)
{
   double max = 0;
   
   int state = gliderPtr->stateSize;

   for(int t = 0; t < T; t++){
      for(int s = 0; s < state; s++){
         if(epst[t][s] > max){
            max = epst[t][s];
         }
      }
   }

   return max;
}

IloAlgorithm::Status CPLEXmodel::solveModel(const Dynamics::dynamics leg, double t0, double tf)
{
   static pthread_mutex_t cs_mutex = PTHREAD_MUTEX_INITIALIZER;
   IloAlgorithm::Status status;  

   this->env = IloEnv();

   IloModel model(env);
   this->SOCP = IloCplex(model);

   this->x = IloNumVarArray(env,T);
   this->y = IloNumVarArray(env,T);
   this->z = IloNumVarArray(env,T);
   this->vel = IloNumVarArray(env,T);
   this->gamma = IloNumVarArray(env,T);
   this->phi = IloNumVarArray(env,T);
   this->CL = IloNumVarArray(env,T);
   this->mu = IloNumVarArray(env,T);
   
   this->xt = IloNumArray(env,T);
   this->yt = IloNumArray(env,T);
   this->zt = IloNumArray(env,T);
   this->vt = IloNumArray(env,T);
   this->gammat = IloNumArray(env,T);
   this->phit = IloNumArray(env,T);
   this->CLt = IloNumArray(env,T);
   this->mut = IloNumArray(env,T);

   //Declaring the error variables
   const int stateSize = gliderPtr->stateSize;
   this->epsilon = IloArray<IloNumVarArray>(env, T);
   this->epst = IloArray<IloNumArray>(env, T);
   for(int t = 0; t < T; t++){
      this->epst[t] = IloNumArray(env, stateSize); 
   }
   this->auxZ = IloNumVarArray(env,T);

   try{
      //create model
      createModel(model,leg,t0,tf);
      //remove all display of barrier algorithm
      SOCP.setOut(env.getNullStream());      
      //set number of threads and cores
      SOCP.setParam(IloCplex::TiLim, 3600);
      SOCP.setParam(IloCplex::Threads, 1);
      SOCP.setParam(IloCplex::ParallelMode, 1);
      SOCP.exportModel("Results/SOCPmodel.lp");
      SOCP.solve();
      status = SOCP.getStatus();      
      //if feasible get solution
      if(status == IloAlgorithm::Optimal){
         //generate output
         SOCP.getValues(xt,x);
         SOCP.getValues(yt,y);
         SOCP.getValues(zt,z);
         SOCP.getValues(vt,vel);
         SOCP.getValues(gammat,gamma);
         SOCP.getValues(phit,phi);
         SOCP.getValues(CLt,CL);
         SOCP.getValues(mut,mu);
         for(int t = 0; t < T; t++){
            SOCP.getValues(epst[t],epsilon[t]);
         }
         //populate solution matrix
         popSolMatrix(xt,yt,zt,vt,gammat,phit,CLt,mut);
         //get solution information
         objVal = SOCP.getObjValue();
         error = errorNorm(epst);
         time = SOCP.getTime();
         tree = SOCP.getNnodes64();
         numIter = SOCP.getNiterations64(); 
      }else{
         cerr << status << endl;
         objVal = IloInfinity;
         time = SOCP.getTime();
         tree = SOCP.getNnodes64();
         numIter = SOCP.getNiterations64();
      }
      this->SOCP.clearModel();
      this->SOCP.end();
      model.end();
      this->env.end();
   }
   catch (IloException& exception) {
      cerr << "Concert exception caught: " << exception << endl;
   }
   catch (...) {
      cerr << "Unknown exception caught" << endl;
   }

   pthread_mutex_unlock(&cs_mutex);
   return status;
}

void CPLEXmodel::finishCplex()
{
   SOCP.end();
}

IloCplex CPLEXmodel::getCplexModel()
{
   return SOCP;
}

void CPLEXmodel::createModel(IloModel model, const Dynamics::dynamics leg, double t0, double tf)
{
   //load data to this routine (angles in rad)
   const double xLb = dataPtr->xLb;
   const double xUb = dataPtr->xUb;
   const double yLb = dataPtr->yLb;
   const double yUb = dataPtr->yUb;
   const double zLb = dataPtr->zLb;
   const double zUb = dataPtr->zUb;
   const double vLb = gliderPtr->vLb;
   const double gammaLb = gliderPtr->gammaLb;
   const double phiLb = gliderPtr->phiLb;
   const double CLlb = gliderPtr->CLlb;
   const double muLb = gliderPtr->muLb;
   const double vUb = gliderPtr->vUb;
   const double gammaUb = gliderPtr->gammaUb;
   const double phiUb = gliderPtr->phiUb;
   const double CLub = gliderPtr->CLub;
   const double muUb = gliderPtr->muUb;
   //make a copy of the system's data
   const int stateSize = gliderPtr->stateSize;
   const int controlSize = gliderPtr->controlSize;
   double const* const* ID = gliderPtr->I;
   const double grav = gliderPtr->gravity;
   const double mass = gliderPtr->mass;   
   const double alpha = gliderPtr->alpha;   
   //initial and final conditions
   double const* Yo = leg.Yo;
   double const* Yf = leg.Yf;
   double const* Uo = leg.Uo;
   double const* vo = leg.vo;
   int const orig = leg.orig;
   
   //update x,y,z components of equilibrium conditions
   const double Yeq[] = {Yo[0],Yo[1],Yo[2],leg.get_Yeq(3),leg.get_Yeq(4),leg.get_Yeq(5)};   
   const double Ueq[] = {leg.get_Ueq(0),leg.get_Ueq(1)};   
   
   //norms of derivatives and matrices
   const double normuDotUb = gliderPtr->uDotUbMaxNorm;
   const double normyDotUb = gliderPtr->yDotUbMaxNorm;
   const double normA = leg.getMaxNormA();
   const double normB = leg.getMaxNormB();

   //step size
   h = (tf - t0)/(double)(T-1);

   //start modelling
   char name[100];

   //create variables corresponding to coordinates x_{t}, y_{t}
   //and z_{t} \forall t \in T
   //================================================================
   //IloNumVarArray x(env,T,xLb,xUb);
   for(int t = 0; t < T; t++){
      IloNumVar temp(env,xLb,xUb,ILOFLOAT);
      x[t] = temp;
      sprintf(name, "x(%d)", t);
      x[t].setName(name);
   }
   //IloNumVarArray y(env,T,yLb,yUb);
   for(int t = 0; t < T; t++){
      IloNumVar temp(env,yLb,yUb,ILOFLOAT);
      y[t] = temp;
      sprintf(name, "y(%d)", t);
      y[t].setName(name);
   }
   //IloNumVarArray z(env,T,zLb,zUb);
   for(int t = 0; t < T; t++){
      IloNumVar temp(env,zLb,zUb,ILOFLOAT);
      z[t] = temp;
      sprintf(name, "z(%d)", t);
      z[t].setName(name);
   } 
   //create state variables vel_{t}, gamma_{t}, phi_{t} and control
   //variables CL_{t} and \mu_{t} \forall t \in T
   //================================================================
   for(int t = 0; t < T; t++){
      IloNumVar temp(env,vLb,vUb,ILOFLOAT);
      vel[t] = temp;
      sprintf(name, "vel(%d)", t);
      vel[t].setName(name);
   }
   for(int t = 0; t < T; t++){
      IloNumVar temp(env,gammaLb,gammaUb,ILOFLOAT);
      gamma[t] = temp;
      sprintf(name, "gamma(%d)", t);
      gamma[t].setName(name);
   }
   for(int t = 0; t < T; t++){
      IloNumVar temp(env,phiLb,phiUb,ILOFLOAT);
      phi[t] = temp;
      sprintf(name, "phi(%d)", t);
      phi[t].setName(name);
   }
   for(int t = 0; t < T; t++){
      IloNumVar temp(env,CLlb,CLub,ILOFLOAT);
      CL[t] = temp;
      sprintf(name, "CL(%d)", t);
      CL[t].setName(name);
      model.add(CL[t]);
   }
   for(int t = 0; t < T; t++){
      IloNumVar temp(env,muLb,muUb,ILOFLOAT);
      mu[t] = temp;
      sprintf(name, "mu(%d)", t);
      mu[t].setName(name);
      model.add(mu[t]);
   }
   
   for(int t = 0; t < T; t++){
      IloNumVar temp(env,0,IloInfinity,ILOFLOAT);
      auxZ[t] = temp;
      sprintf(name, "z_aux(%d)", t);
      auxZ[t].setName(name);
      model.add(auxZ[t]);
   }
   
   //create the epsilon and aux variables for computing the
   //discretisation errors
   for(int t = 0; t < T; t++){
      epsilon[t] = IloNumVarArray(env, stateSize); 
      for(int s = 0; s < stateSize; s++){
         epsilon[t][s] = IloNumVar(env, -IloInfinity, IloInfinity, ILOFLOAT);
         sprintf(name, "_eps(%d,%d)", t, s);
         epsilon[t][s].setName(name);
         model.add(epsilon[t][s]);
      }
   }

   //create objective function
   IloExpr objFunction(env);
   for(int t = 0; t < T; t++){
      objFunction += auxZ[t];
   }
   model.add(IloMinimize(env,objFunction));

   //error bound constraints
   //===============================================
   for(int t = 0; t < T-1; t++){
      for(int s = 0; s < stateSize; s++){
         IloRange cons = (auxZ[t] - epsilon[t][s] >= 0);
         sprintf(name,"auxZ_mod1_(%d,%d)",t,s);
         cons.setName(name);
         model.add(cons);
      }
   } 
   for(int t = 0; t < T-1; t++){
      for(int s = 0; s < stateSize; s++){
         IloRange cons = (auxZ[t] + epsilon[t][s] >= 0);
         sprintf(name,"auxZ_mod2_(%d,%d)",t,s);
         cons.setName(name);
         model.add(cons);
      }
   } 
   for(int t = 0; t < T-1; t++){
      for(int s = 0; s < stateSize; s++){
         IloRange cons = (auxZ[t] <= 0.5*h*h*(normA*normyDotUb + normB*normuDotUb));
         sprintf(name,"auxZ_mod2_(%d,%d)",t,s);
         cons.setName(name);
         model.add(cons);
      }
   } 
 
   //set of system dynamics constraints
   //===============================================
   for(int t = 0; t < T-1; t++){
      vector <IloNumVar> Yt = {x[t],y[t],z[t],vel[t],gamma[t],phi[t]};
      vector <IloNumVar> Ytp1 = {x[t+1],y[t+1],z[t+1],vel[t+1],gamma[t+1],phi[t+1]};
      vector <IloNumVar> Ut = {CL[t],mu[t]};
      for(int i = 0; i < stateSize; i++){
         IloExpr exp(env);
         exp -= Ytp1[i];
         for(int j = 0; j < stateSize; j++){
               exp += (h*leg.get_A(i,j) + ID[i][j])*Yt[j] - h*leg.get_A(i,j)*Yeq[j];
         }
         for(int j = 0; j < controlSize; j++){
            exp += h*leg.get_B(i,j)*Ut[j] - h*leg.get_B(i,j)*leg.get_Ueq(j);
         }
         exp -= epsilon[t][i];
         IloRange cons = (exp <= 0);
         sprintf(name,"sys1(A%d,)_t%d",i,t);
         cons.setName(name);
         model.add(cons);
      }
   }   

   for(int t = 0; t < T-1; t++){
      vector <IloNumVar> Yt = {x[t],y[t],z[t],vel[t],gamma[t],phi[t]};
      vector <IloNumVar> Ytp1 = {x[t+1],y[t+1],z[t+1],vel[t+1],gamma[t+1],phi[t+1]};
      vector <IloNumVar> Ut = {CL[t],mu[t]};
      for(int i = 0; i < stateSize; i++){
         IloExpr exp(env);
         exp -= Ytp1[i];
         for(int j = 0; j < stateSize; j++){
               exp += (h*leg.get_A(i,j) + ID[i][j])*Yt[j] - h*leg.get_A(i,j)*Yeq[j];
         }
         for(int j = 0; j < controlSize; j++){
            exp += h*leg.get_B(i,j)*Ut[j] - h*leg.get_B(i,j)*leg.get_Ueq(j);
         }
         exp += epsilon[t][i];
         IloRange cons = (exp >= 0);
         sprintf(name,"sys2(A%d,)_t%d",i,t);
         cons.setName(name);
         model.add(cons);
      }
   }   
   //end dynamics constraints
   //===============================================

   ////===============================================
   ////boundary conditions
   ////define initial state contraints and add to model
   x[0].setBounds(Yo[0],Yo[0]);
   y[0].setBounds(Yo[1],Yo[1]);
   z[0].setBounds(Yo[2],Yo[2]);
   vel[0].setBounds(Yo[3],Yo[3]);
   gamma[0].setBounds(Yo[4],Yo[4]);
   phi[0].setBounds(Yo[5],Yo[5]);
   //define initial control contraints and add to model
   CL[0].setBounds(Uo[0],Uo[0]);
   mu[0].setBounds(Uo[1],Uo[1]);
   //end boundary conditions
   //===============================================
   
   //add volume/visiting constraints
   ////===============================================
   IloNumVar aux1(env, -IloInfinity, IloInfinity, ILOFLOAT);
   IloNumVar aux2(env, -IloInfinity, IloInfinity, ILOFLOAT);
   IloNumVar aux3(env, -IloInfinity, IloInfinity, ILOFLOAT);
   IloNumVar rad(env, 0, IloInfinity, ILOFLOAT);
   model.add(x[T-1] - Yf[0] == aux1);
   model.add(y[T-1] - Yf[1] == aux2);
   model.add(z[T-1] - Yf[2] == aux3);
   if(!leg.isLanding){
      //add cone constraints
      model.add(rad == (z[T-1] + leg.fRadius)*tan(alpha));
      IloConstraint vol = (aux1*aux1 + aux2*aux2 - rad*rad <= 0);
      z[T-1].setBounds(leg.fzMin,leg.fzMax);
      sprintf(name,"final_cone");
      vol.setName(name);
      model.add(vol);
      //define level flight (photographing) constraints
      gamma[T-1].setBounds(-0.09,0.09);
      mu[T-1].setBounds(-0.09,0.09);
   }else{
      IloConstraint vol = (aux1*aux1 + aux2*aux2 + aux3*aux3 <= leg.fRadius);
      sprintf(name,"landing_region");
      vol.setName(name);
      model.add(vol);
   }
   //end visiting/volume constraits
   //===============================================
   
   //add small perturbation constraints
   ////===============================================
   double c_eps = 0.001;
   for(int t = 0; t < T; t++){
      double term1 = t/(double)(T-1);
      double term2 = 1 - t/(double)(T-1);
      double glb_CL = term1*(Ueq[0] - c_eps) + term2*CLlb;
      double gub_CL = term1*(Ueq[0] + c_eps) + term2*CLub;
      CL[t].setBounds(glb_CL, gub_CL);
      
      double glb_mu = term1*(Ueq[1] - c_eps) + term2*muLb;
      double gub_mu = term1*(Ueq[1] + c_eps) + term2*muUb;
      mu[t].setBounds(glb_mu, gub_mu);
   }
   ////===============================================

   //finished modelling
}

void CPLEXmodel::popSolMatrix(IloNumArray xt, IloNumArray yt, IloNumArray zt, IloNumArray vt, IloNumArray gammat, IloNumArray phit, IloNumArray CLt, IloNumArray mut)
{
   //allocate memory to the solution matrix
   int state = gliderPtr->stateSize;
   int ctrl = gliderPtr->controlSize;
   solMatrix = new double* [T];
   for(int t = 0; t < T; t++){
      solMatrix[t] = new double [state+ctrl];
   }

   for(int t = 0; t < T; t++){
      solMatrix[t][0] = xt[t];
      solMatrix[t][1] = yt[t];
      solMatrix[t][2] = zt[t];
      solMatrix[t][3] = vt[t];
      solMatrix[t][4] = gammat[t];
      solMatrix[t][5] = phit[t];
      solMatrix[t][6] = CLt[t];
      solMatrix[t][7] = mut[t];
   }
}

