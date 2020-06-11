#include "CPLEXmodel.h"

ILOSTLBEGIN // import namespace std

CPLEXmodel::CPLEXmodel(Data* data, Dynamics* glider, int T)
   :dataPtr(data), gliderPtr(glider)
   //xt(env), yt(env), zt(env), x(env,T), y(env,T), z(env,T),
   //vt(env), gammat(env), phit(env), CLt(env), mut(env),
   //vel(env,T), gamma(env,T), phi(env,T), CL(env,T), mu(env,T),
   //epsilon(env, 0, IloInfinity, ILOFLOAT)
{
   this->T = dataPtr->T;
}	
CPLEXmodel::~CPLEXmodel()
{
   env.end();
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

   this->epsilon = IloNumVar(env, 0, IloInfinity, ILOFLOAT, "_eps");
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
         error = SOCP.getValue(epsilon);
         //populate solution matrix
         popSolMatrix(xt,yt,zt,vt,gammat,phit,CLt,mut);
         //get solution information
         objVal = SOCP.getObjValue();
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

   //create objective function
   IloExpr objFunction(env);
   objFunction == 0;
   objFunction += epsilon;
   model.add(IloMinimize(env,objFunction));
 
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
         exp -= epsilon;
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
         exp += epsilon;
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
      double term1 = t/(double)T;
      double term2 = (T-t)/(double)T;
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


   //double* vf = new double[3];

   //double xf = solMatrix[T-1][0];
   //double yf = solMatrix[T-1][1];
   //double zf = solMatrix[T-1][2];
   //double xi = solMatrix[T-2][0];
   //double yi = solMatrix[T-2][1];
   //double zi = solMatrix[T-2][2];
   //vf[0] = (xf - xi)/h; 
   //vf[1] = (yf - yi)/h; 
   //vf[2] = (zf - zi)/h; 

   //cout << "FINAL VEL VECTOR= [" << vf[0] << "," << vf[1] << "," << vf[2] << "]" << endl;
   //getchar();






   ////avoid gamma oscilation
   //IloNumVarArray deltaGamma(env,T-1,-IloInfinity,IloInfinity,ILOFLOAT);   
   //for(int t = 0; t < T-1; t++){
      //sprintf(name,"dGamma(%d,%d)",t+1,t);
      //deltaGamma[t].setName(name);
   //}   
   //for(int t = 0; t < T-1; t++){
      //IloRange cons = (deltaGamma[t] -gamma[t+1] + gamma[t] == 0);
      //sprintf(name,"diff_Gamma(%d)",t);
      //cons.setName(name);
      //model.add(cons);
   //}
   //for(int t = 0; t < T-1; t++){
      //IloRange cons = (deltaGamma[t]*deltaGamma[t] <= greekXi*greekXi);
      //sprintf(name,"dlt_Gamma(%d)",t);
      //cons.setName(name);
      //model.add(cons);
   //}
   
   ////avoid phi oscilation
   //IloNumVarArray deltaPhi(env,T-1,-IloInfinity,IloInfinity,ILOFLOAT);   
   //for(int t = 0; t < T-1; t++){
      //sprintf(name,"dPhi(%d,%d)",t+1,t);
      //deltaGamma[t].setName(name);
   //}   
   //for(int t = 0; t < T-1; t++){
      //IloRange cons = (deltaPhi[t] -phi[t+1] + phi[t] == 0);
      //sprintf(name,"diff_Phi(%d)",t);
      //cons.setName(name);
      //model.add(cons);
   //}
   //for(int t = 0; t < T-1; t++){
      //IloRange cons = (deltaPhi[t]*deltaPhi[t] <= greekXi*greekXi);
      //sprintf(name,"dlt_Phi(%d)",t);
      //cons.setName(name);
      //model.add(cons);
   //}

   ////avoid V oscilation
   //IloNumVarArray deltaV(env,T-1,-IloInfinity,IloInfinity,ILOFLOAT);   
   //for(int t = 0; t < T-1; t++){
      //sprintf(name,"dV(%d,%d)",t+1,t);
      //deltaV[t].setName(name);
   //}   
   //for(int t = 0; t < T-1; t++){
      //IloRange cons = (deltaV[t] -vel[t+1] + vel[t] == 0);
      //sprintf(name,"diff_V(%d)",t);
      //cons.setName(name);
      //model.add(cons);
   //}
   //for(int t = 0; t < T-1; t++){
      //IloRange cons = (deltaV[t]*deltaV[t] <= greekXi*greekXi);
      //sprintf(name,"dlt_V(%d)",t);
      //cons.setName(name);
      //model.add(cons);
   //}
   
   //minimise distance from central axis of target
   //add distance SOCP constraint
   ////===============================================
   //IloConstraint distance = (aux1*aux1 + aux2*aux2 + -fDist*fDist <= 0);
   //sprintf(name,"final_distance");
   //distance.setName(name);
   //model.add(distance);
   ////===============================================

   ////===============================================
   //stabilisation constraints
   //first order smoothness conditions
   //IloExpr p1Exp(env);
   //for(int t = 0; t < T-1; t++){
      //p1Exp += ((CL[t+1]-CL[t])*(CL[t+1]-CL[t]))/h;
   //}
   //IloRange p1 = (pi1 - p1Exp >= 0);
   //sprintf(name, "pi1");
   //p1.setName(name);
   //model.add(p1);    

   //IloExpr p2Exp(env);
   //for(int t = 0; t < T-1; t++){
      //p2Exp += ((mu[t+1]-mu[t])*(mu[t+1]-mu[t]))/h;
   //}
   //IloRange p2 = (pi2 - p2Exp >= 0);
   //sprintf(name, "pi2");
   //p2.setName(name);
   //model.add(p2);

   //IloExpr p3Exp(env);
   //for(int t = 0; t < T-1; t++){
      //p3Exp += ((phi[t+1]-phi[t])*(phi[t+1]-phi[t]))/h;
   //}
   //IloRange p3 = (pi3 - p3Exp >= 0);
   //sprintf(name, "pi3");
   //p3.setName(name);
 
   //model.add(p3);
   //IloExpr p3Exp(env);
   //for(int t = 0; t < T-1; t++){
      //p3Exp += h*h*(((phi[t+1]+phi[t])/2)*((phi[t+1]+phi[t])/2));
   //}
   //IloRange p3 = (pi3 - p3Exp >= 0);
   //sprintf(name, "pi3");
   //p3.setName(name);
   //model.add(p3);

   //avoid CL oscilation
/*
   IloNumVarArray deltaCL(env,T-1,-IloInfinity,IloInfinity,ILOFLOAT);   
   for(int t = 0; t < T-1; t++){
      sprintf(name,"dCL(%d,%d)",t+1,t);
      deltaCL[t].setName(name);
   }   
   for(int t = 0; t < T-1; t++){
      IloRange cons = (deltaCL[t] -CL[t+1] + CL[t] == 0);
      sprintf(name,"diff_CL(%d)",t);
      cons.setName(name);
      model.add(cons);
   }
   for(int t = 0; t < T-1; t++){
      IloRange cons = (deltaCL[t]*deltaCL[t] <= greekXi*greekXi);
      sprintf(name,"dlt_CL(%d)",t);
      cons.setName(name);
      model.add(cons);
   }
   
   //avoid mu oscilation
   IloNumVarArray deltaMu(env,T-1,-IloInfinity,IloInfinity,ILOFLOAT);   
   for(int t = 0; t < T-1; t++){
      sprintf(name,"dMu(%d,%d)",t+1,t);
      deltaMu[t].setName(name);
   }   
   for(int t = 0; t < T-1; t++){
      IloRange cons = (deltaMu[t] -mu[t+1] + mu[t] == 0);
      sprintf(name,"diff_mu(%d)",t);
      cons.setName(name);
      model.add(cons);
   }
   for(int t = 0; t < T-1; t++){
      IloRange cons = (deltaMu[t]*deltaMu[t] <= greekXi*greekXi);
      sprintf(name,"dlt_mu(%d)",t);
      cons.setName(name);
      model.add(cons);
   }
*/
   //end stabilisation constraints
   //===============================================

   //=============================================
