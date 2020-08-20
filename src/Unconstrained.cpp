#include "Unconstrained.h"


Unconstrained::Unconstrained()
{

}
Unconstrained::~Unconstrained()
{
}

IloAlgorithm::Status 
Unconstrained::minFlightTime(CPLEXmodel* modelPtr, Dynamics::dynamics leg, double tf_lb, double tf_ub, double& tf, double& step, double& error)
{
   //=========================================================
#if DEBUG > 0
   cout << "\nFinding minimum flight time from LB";
#endif
   IloAlgorithm::Status status = IloAlgorithm::Infeasible;
   int iter = 0;
   int preIter = 0;

   double t0 = 0;
   double timeCnt = 0;
   double ta = tf_lb;
   double tb = tf_ub;
   double tOpt = 0;
   
   double stepa = 0;
   double stepb = 0;
   double stepOpt = 0;

   double errora = 0;
   double errorb = 0;
   double errorOpt = 0;

   double fa = INF;
   double fb = INF;
   double fOpt = INF;

   //check time interval [ta,tb]
#if DEBUG > 0
   printf("\nInitial time intervals [%.1f %.1f] and [%.1f %.1f]\n",t0,ta,t0,tb);
#endif
   if(tb-ta <= 1){
      cout << "Time interval too small!" << endl;
      return IloAlgorithm::Infeasible;
   }

   //check feasibility of tb
   status = modelPtr->solveModel(leg,t0,tb);
   preIter++;
   if(status == IloAlgorithm::Optimal){
      timeCnt += modelPtr->time;
      fb = penal*tb + modelPtr->objVal;
      stepb = modelPtr->h;
      errorb = modelPtr->error;
   }

   //terminate if tb is infeasible
   if(status != IloAlgorithm::Optimal){
      printf("ta \t tb \t fa \t\t fb\n");
      printf("%.1f \t %.1f \t %f \t %f \n",ta,tb,fa,fb);
      cout << "CPLEX infeasible solution for tb!" << endl;
      return status;
   }
   
   //reset status
   status = IloAlgorithm::Infeasible;
   
   //find a feasible ta
   while(fa >= INF && ta <= tb){
      status = modelPtr->solveModel(leg,t0,ta);
      preIter++;
      if(status == IloAlgorithm::Optimal){
         timeCnt += modelPtr->time;
         fa = penal*ta + modelPtr->objVal;
         stepa = modelPtr->h;
         errora = modelPtr->error;
      }else{
         fa = INF;
         ta++;
      }      
   }
   
   //copy optimal solution and respective step size
   if(status == IloAlgorithm::Optimal){
      fOpt = fa;
      tf = ta;
      step = stepa;
      error = errora;
#if DEBUG > 0
      cout << "=====================================================" << endl;
      cout << "OPTIMAL! "<< endl;
      cout << "Solution (t*): " << tf  << " | Step Size= " << step << " | error= " << error << endl;
      cout << "Obj. Value (f*): " << fOpt << endl;
      cout << "CPLEX CPU time (s): " << timeCnt << endl;
#endif
      return IloAlgorithm::Optimal;
   }else{
      //return infeasible if solution is not optimal
#if DEBUG > 0
      cout << "=====================================================" << endl;
      cout << "INFEASIBLE! "<< endl;
      cout << "Solution (t*): " << tOpt << endl;
#endif
      status = IloAlgorithm::Infeasible;
      return IloAlgorithm::Infeasible;
   }
}

IloAlgorithm::Status 
Unconstrained::bisection(CPLEXmodel* modelPtr, Dynamics::dynamics leg, double tf_lb, double tf_ub, double& tf, double& step, double& error)
{
   //=========================================================
   //bisection initialisation
   cout << "\nBissection method fligt time optimisation";
   //int status = 0;
   IloAlgorithm::Status status = IloAlgorithm::Infeasible;
   int iter = 0;
   int preIter = 0;

   double t0 = 0;
   double timeCnt = 0;
   double ta = tf_lb;
   double tb = tf_ub;
   double tm = 0.5*(ta+tb);
   double tl = 0;
   double tr = 0;
   double tOpt = 0;
   
   double stepa = 0;
   double stepb = 0;
   double stepm = 0;
   double stepl = 0;
   double stepr = 0;
   double stepOpt = 0;

   double errora = 0;
   double errorb = 0;
   double errorm = 0;
   double errorl = 0;
   double errorr = 0;
   double errorOpt = 0;

   double fa = INF;
   double fb = INF;
   double fm = INF;
   double fl = INF;
   double fr = INF;
   double fOpt = INF;

   printf("\nInitial time intervals [%.1f %.1f] and [%.1f %.1f]\n",t0,ta,t0,tb);
   if(tb-ta <= 1){
      cout << "Time interval too small!" << endl;
      return IloAlgorithm::Infeasible;
   }

   //status = modelPtr->solveModel(leg,t0,ta);
   //preIter++;
   //if(status == IloAlgorithm::Optimal){
      //timeCnt += modelPtr->time;
      //fa = penal*ta + modelPtr->objVal;
   //}   
   status = modelPtr->solveModel(leg,t0,tb);
   preIter++;
   if(status == IloAlgorithm::Optimal){
      timeCnt += modelPtr->time;
      fb = penal*tb + modelPtr->objVal;
      stepb = modelPtr->h;
      errorb = modelPtr->error;
   }

   if(status != IloAlgorithm::Optimal){
      printf("ta \t tb \t fa \t\t fb\n");
      printf("%.1f \t %.1f \t %f \t %f \n",ta,tb,fa,fb);
      cout << "CPLEX infeasible solution for both t bounds!" << endl;
      return status;
   }
   
   while(fa >= INF){
      ta++;
      status = modelPtr->solveModel(leg,t0,ta);
      preIter++;
      if(status == IloAlgorithm::Optimal){
         timeCnt += modelPtr->time;
         fa = penal*ta + modelPtr->objVal;
         stepa = modelPtr->h;
         errora = modelPtr->error;
      }else fa = INF;
   }
   //while(fb >= INF){
      //tb++;
      //status = modelPtr->solveModel(leg,t0,tb);
      //preIter++;
      //if(status == IloAlgorithm::Optimal){
         //timeCnt += modelPtr->time;
         //fb = penal*tb + modelPtr->objVal;
      //}else fb = INF;
   //}
   tm = 0.5*(ta+tb);
   status = modelPtr->solveModel(leg,t0,tm);
   preIter++;
   if(status == IloAlgorithm::Optimal){
      timeCnt += modelPtr->time;
      fm = penal*tm + modelPtr->objVal;
      stepm = modelPtr->h;
      errorm = modelPtr->error;
   }else{
      ta = INF;
      tb = INF;
      cout << "Infeasible bounds" << endl;
   }

   cout << "\n\n-----------------------------------------------------" << endl;
   cout << "Pre-processing phase - New bounds:" << endl;
   printf("ta \t tb \t fa \t\t fb \t\t iter\n");
   printf("%.1f \t %.1f \t %f \t %f \t %d \n",ta,tb,fa<INF?fa:-1.0,fb<INF?fb:-1.0,preIter);
   cout << "-----------------------------------------------------" << endl;
   cout << "CPU (s): " << timeCnt << endl;
   cout << endl;

   if(abs(tb-ta) <= 1){
      if(fa < fb){
         fOpt = fa;
         tOpt = ta;
         stepOpt = stepa;
         errorOpt = errora;
      }else{
         fOpt = fb;
         tOpt = tb;
         stepOpt = stepb;
         errorOpt = errorb;
      }
   }

   //bisection main loop
   cout << "Bisection-CPLEX" << endl;
   cout << "=====================================================" << endl;
   printf("ta \t tb \t fa \t\t fb \t\t iter\n");
   cout << "=====================================================" << endl;
   while(abs(tb-ta) > 1 && tb < maxtb){
      iter++;
      printf("%.1f \t %.1f \t %f \t %f \t %d \n",ta,tb,fa<INF?fa:-1.0,fb<INF?fb:-1.0,iter);
      tl = 0.5*(ta+tm);
      tr = 0.5*(tb+tm);
      status = modelPtr->solveModel(leg,t0,tl);
      if(status == IloAlgorithm::Optimal){
         timeCnt += modelPtr->time;
         fl = penal*tl + modelPtr->objVal;
         stepl = modelPtr->h;
         errorl = modelPtr->error;
      }else fl = INF;
      status = modelPtr->solveModel(leg,t0,tr);
      if(status == IloAlgorithm::Optimal){
         timeCnt += modelPtr->time;
         fr = penal*tr + modelPtr->objVal;
         stepr = modelPtr->h;
         errorr = modelPtr->error;
      }else fr = INF;
      double fAll[] = {fa,fb,fm,fl,fr};
      fOpt = *min_element(fAll,fAll+5);
      if(fOpt == fa||fOpt == fl){
         if(fOpt == fa){
            tOpt = ta;
            stepOpt = stepa;
            errorOpt = errora;
         }
         if(fOpt == fl){
            tOpt = tl;
            stepOpt = stepl;
            errorOpt = errorl;
         }
         tb = tm;
         tm = tl;
         fb = fm;
         fm = fl;
      }else{
         if(fOpt == fm){
            tOpt = tm;
            stepOpt = stepm;
            errorOpt = errorm;
            ta = tl;
            tb = tr;
            fa = fl;
            fb = fr;
         }else{
            if(fOpt == fr || fOpt == fb){
               if(fOpt == fr){
                  tOpt = tr;
                  stepOpt = stepr;
                  errorOpt = errorr;
                  ta = tm;
                  tm = tr;
                  fa = fm;
                  fm = fr;
               }
               if(fOpt == fb){
                  tOpt = tb;
                  stepOpt = stepb;
                  errorOpt = errorb;
                  tb += tm;
                  ta = tm;
                  tm = tb;
                  fa = fm;
                  fm = fb;
                  //fb = INF;
               }
            }
         }
      }
      //cout << tOpt << " | " << fOpt << endl;
   }
   //copy optimal solution and respective step size
   tf = tOpt;
   step = stepOpt;
   error = errorOpt;

   if(fb < 1000){
      cout << "=====================================================" << endl;
      cout << "OPTIMAL! "<< endl;
      cout << "Solution (t*): " << tf  << " | Step Size= " << step << endl;
      cout << "Obj. Value (f*): " << fOpt << endl;
      cout << "CPLEX CPU time (s): " << timeCnt << endl;
   }else{
      cout << "=====================================================" << endl;
      cout << "INFEASIBLE! "<< endl;
      cout << "Solution (t*): " << tOpt << endl;
      status = IloAlgorithm::Infeasible;
      exit(1);
   }
   //end of bissection method
   //=========================================================
   //getchar();
   return status;
}

void Unconstrained::
plotObjEval(Data* dataPtr, CPLEXmodel* modelPtr, Dynamics::dynamics leg, double tf_lb, double tf_ub)
{
   cout << "\n\nPlotting Objective Function Shape!" << endl;
   //int instSize = 1 + dataPtr->noWaypoints + dataPtr->noLandPoints;
   //if(instSize > 2){
      //cout << "Do not call this function for |V| > 2!" << endl;
      //exit(1);
   //}
   double t0 = 0;
   IloAlgorithm::Status status;
   CPLEXmodel model = *modelPtr;
   ofstream file("Results/objFunc.out", ios::out);
   file << dataPtr->fileName << endl;
   for(double t = tf_lb; t <= tf_ub; t++){
      status = model.solveModel(leg,t0,t);
      double cplex = model.objVal;
      double obj = penal*t + cplex;
      if(cplex > DBL_MAX -1) cplex = -1;
      if(obj > DBL_MAX -1) obj = -1;
      file << t <<" "<< obj <<" "<< cplex << endl;
      cout << t <<" "<< obj <<" "<< cplex << endl;
   }
   file.close();
   system("python ../Python/plotObjFunc.py");
}

//libcmaes::FitFunc objFunction = [](const double* deltaT, const int N = 1)
//{
   //static pthread_mutex_t cs_mutex = PTHREAD_MUTEX_INITIALIZER;
   
   //int t0 = 0;
   //double tf = t0 + *deltaT;
   //statusGlobal = modelGlobal->solveModel(legGlobal,t0,tf);
   //return *deltaT + modelGlobal->objVal;
   
   //pthread_mutex_unlock(&cs_mutex);
//};

//ProgressFunc<CMAParameters<GenoPheno<pwqBoundStrategy>>,CMASolutions> 
//select_time = [](const CMAParameters<GenoPheno<pwqBoundStrategy>> &cmaparams, const CMASolutions &cmasols)
//{
   //if(cmasols.niter() == 0){
      //cout << "\nCMA-ES flight time optimisation";
      //printf("\n---------------------------------------");
      //printf("\n    obj       sol    iter  funEv  time");
      //printf("\n---------------------------------------\n");
   //}else{
      //int printStep = 10;
      //if(cmasols.niter() % printStep == 0 || cmasols.niter() == 1){
         //double objVal = fabs(cmasols.best_candidate().get_fvalue());
         //int nIter = cmasols.niter();
         //int nFuncEvals = cmasols.fevals();
         //double time = printStep*(double)cmasols.elapsed_last_iter()/1000.0; 
         //vector <double> sol = cmasols.get_best_seen_candidate().get_x();
         //int sol_size = cmasols.best_candidate().get_x_size();
         //printf("%10.6f  %6.4f  %.4d  %.4d  %6.4f \n",objVal,sol[0],nIter,nFuncEvals,time);
      //}
   //}
  //return 0;
//};

//void Unconstrained::callCMAES(CPLEXmodel* model, Dynamics::dynamics leg, double deltaT0, double& tOpt)
//{
   //static pthread_mutex_t cs_mutex = PTHREAD_MUTEX_INITIALIZER;
   ////####### starting CMAES ##################
   //double t0 = 0;
   ////copy model data to global variable
   //modelGlobal = model;
   //legGlobal = leg;

   //int dim = 1; // problem dimensions.
   //vector<double> flightTime(dim,deltaT0);
   //// arrays for lower and upper parameter bounds, respectively         
   //double sigma = 1;
   //double lbounds[dim];
   //double ubounds[dim];
   //for(int i = 0; i < dim; i++){
      //lbounds[i] = 1;
      //ubounds[i] = DBL_MAX;
    //}
   ////set lb and ub objects
   //// genotype / phenotype transform associated to bounds.
   //GenoPheno<pwqBoundStrategy> gp(lbounds,ubounds,dim);  
   //// -1 for automatically\decided lambda, 0 is for random seeding of the internal generator.
   //CMAParameters<GenoPheno<pwqBoundStrategy>> cmaparams(flightTime,sigma,-1,0,gp);
   //cmaparams.set_max_fevals(150);
   //cmaparams.set_ftolerance(0.01);
   //cmaparams.set_algo(aCMAES);
   ////call solver, sort candidates and print opt. sol
   //CMASolutions cmasols = cmaes<GenoPheno<pwqBoundStrategy>>(objFunction,cmaparams,select_time);
   //cmasols.sort_candidates();
   //cout << "\n\nbest solution: " << cmasols << endl;
   //cout << "optimization took " << cmasols.elapsed_time() / 1000.0 << " seconds\n";

   //double tf = t0 + (cmasols.best_candidate().get_fvalue() - modelGlobal->objVal);
   //cout << "Minimum flight time: " << tf << endl;
   //tOpt = tf;
   //int status = statusGlobal;
   ////###### finishing CMAES ##################
   //pthread_mutex_unlock(&cs_mutex);
//}






















