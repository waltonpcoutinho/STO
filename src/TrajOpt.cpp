#include "TrajOpt.h"

TrajOpt::TrajOpt(Data* dataPtr, Dynamics* gliderPtr, CPLEXmodel* modelPtr, AMPLmodel* nlpPtr, string usrMethod)
:dynPtr(gliderPtr), data(dataPtr), model(modelPtr), nlp(nlpPtr), method(usrMethod) 
{
   if(method != "STO" and method != "STO-NLP" and method != "STO+STO-NLP"){
      cout << method <<" is an invalid option!" << endl;
      cout << "Possible options are 'STO', 'STO-NLP' or STO+STO-NLP" << endl;
      exit(0);
   }

   if(method == "STO"){
      cout << "\n\nRunning STO!" << endl;
   }
   if(method == "STO-NLP"){
      cout << "\n\nRunning STO-NLP!" << endl;
   }
   if(method == "STO+STO-NLP"){
      cout << "\n\nRunning STO + STO-NLP!" << endl;
   }
}

// destrutor
TrajOpt::~TrajOpt()
{
}

IloAlgorithm::Status TrajOpt::findTraj(vector<int>& route, int seqSize, vector<vector<double>>& solution, double* flightTimes, double* stepSizes, double* errors, int& infsblArc, double** normEpsilons, 
double** normTaylor1st, double** relEpsilons)
{
   //##################################################
   //check route feasibility and compute optimal trajectories
   IloAlgorithm::Status status = IloAlgorithm::Infeasible;
   const int T = data->T;
   const int stateSize = dynPtr->stateSize;
   const int controlSize = dynPtr->controlSize;
   bool feasible = dynPtr->checkRouteStatus(route,seqSize);
   if(feasible){
      //initialisation
      const Dynamics::dynamics leg0 = dynPtr->getLeg(route[0],route[1]);
      double Ybar[stateSize] = {leg0.Yo[0], leg0.Yo[1], leg0.Yo[2],
                                leg0.get_Yeq(3),leg0.get_Yeq(4),leg0.get_Yeq(5)};
      double Ubar[controlSize] = {leg0.get_Ueq(0),leg0.get_Ueq(1)};      
      double vfAux[3] = {0, 0, 0};
      
      int counter = 0;

      //main loop of Traj. Opt.
      for(int k = 0; k < seqSize - 1; k++){
         //Print some log
         cout << "Iteration: " << k
            << " | Arc: (" << route[k] << "," << route[k+1] << ")" << endl;

         //get first first leg
         const Dynamics::dynamics leg = dynPtr->getLeg(route[k],route[k+1]);
#if DEBUG > 0
         cout << "Origin: "  << leg.orig << endl;
         cout << "Dest: "  << leg.dest << endl;
#endif
         //double targetVel = 0;
         //update boundary conditions
         //states
         double Yo_aux[] = {Ybar[0],Ybar[1],Ybar[2],Ybar[3],Ybar[4],Ybar[5]};
         for(int ii = 0; ii < stateSize; ii++){
            leg.Yo[ii] = Yo_aux[ii];
         }
         for(int ii = 3; ii < stateSize; ii++){
            leg.Yf[ii] = 0;
         }
         //controls
         leg.Uo[0] = Ubar[0]; 
         leg.Uo[1] = Ubar[1];
         //velocites
         leg.vo[0] = vfAux[0];
         leg.vo[1] = vfAux[1];
         leg.vo[2] = vfAux[2];
         
         //print boundary conditions
         const double Yeq[] = {leg.Yo[0],leg.Yo[1],leg.Yo[2],
                               leg.get_Yeq(3),leg.get_Yeq(4),leg.get_Yeq(5)};
#if DEBUG > 0
         cout << "\nBoundary conditions: " << endl;
         printVec(leg.vo,3,"vo: ");
         printVec(Yeq,stateSize,"Yeq: ");
         printVec(leg.Yo,stateSize,"Yo: ");
         printVec(leg.Uo,controlSize,"Uo: ");
         printVec(leg.Yf,stateSize,"Yf: ");
#endif
         //define some aux variables
         double tOpt = 0;
         double stepOpt = 0;
         double errorOpt;
         status = IloAlgorithm::Infeasible;

         if(method == "STO"){
            /*
             * START STO
             */
            //run the unconstrained optimization
            Unconstrained* uncOpt = new Unconstrained();
            double tfLb = leg.deltaTa;
            double tfUb = leg.deltaTb;
            double tfGuess = (tfUb + tfLb)/2;
            double deltaT = tfUb - tfLb;

            //check Ub for tf against maxtb
            if(!(tfUb < maxtb)){
               cout << "tfUb= " << tfUb << " , maxtb= " << maxtb << endl;
               cout << "maxtb too small?" << endl;               
            }

            //compute optimal trajectory between waypoint k and k+1
            while(status != IloAlgorithm::Optimal && tfUb < maxtb){
               //call optimiser
               status = uncOpt->minFlightTime(model,leg,tfLb,tfUb,tOpt,stepOpt,errorOpt);
               //shift time interval if status != optimal
               if(status != IloAlgorithm::Optimal){
                  tfLb += deltaT;
                  tfUb += deltaT;
               }
            }
            /*
             * FINISH STO
             */
         }else if(method == "STO-NLP"){
            /*
             * START STO-NLP
             */
            //call NLP solver
            int statusNLP = nlp->solveNLP(leg, tOpt, stepOpt, errorOpt);

            if(statusNLP == 1){
               status = IloAlgorithm::Optimal;
            }else{
               status = IloAlgorithm::Infeasible;
            }
            /*
             * FINISH STO-NLP
             */
         }else if(method == "STO+STO-NLP"){
            /*
             * START STO+STO-NLP
             */
            //run the unconstrained optimization
            Unconstrained* uncOpt = new Unconstrained();
            double tfLb = leg.deltaTa;
            double tfUb = leg.deltaTb;
            double tfGuess = (tfUb + tfLb)/2;
            double deltaT = tfUb - tfLb;

            //check Ub for tf against maxtb
            if(!(tfUb < maxtb)){
               cout << "tfUb= " << tfUb << " , maxtb= " << maxtb << endl;
               cout << "maxtb too small?" << endl;               
               //exit(1);
            }

            //compute optimal trajectory between waypoint k and k+1
            while(status != IloAlgorithm::Optimal && tfUb < maxtb){
               //call optimiser
               status = uncOpt->minFlightTime(model,leg,tfLb,tfUb,tOpt,stepOpt,errorOpt);
               //shift time interval if status != optimal
               if(status != IloAlgorithm::Optimal){
                  tfLb += deltaT;
                  tfUb += deltaT;
               }
            /*
             * FINISH STO+STO-NLP
             */
            }
            
            if(status != IloAlgorithm::Optimal){
               status = IloAlgorithm::Infeasible;
               break;
            }

            //set up STO solution as initial guess for STO-NLP
            for(int ii = 0; ii < T; ii++){
               for(int jj = 0; jj < stateSize + controlSize; jj++){
                  nlp->solMatrix[ii][jj] = model->solMatrix[ii][jj];
               }
            }
            
            //call NLP solver
            int statusNLP = nlp->solveNLP(leg, tOpt, stepOpt, errorOpt);

            if(statusNLP == 1){
               status = IloAlgorithm::Optimal;
            }
         }

         //save solution if optimal
         //return infeasible otherwise
         if(status == IloAlgorithm::Optimal){
            if(method == "STO"){
               for(int ii = 0; ii < T; ii++){
                  for(int jj = 0; jj < stateSize + controlSize; jj++){
                     solution[counter][jj] = model->solMatrix[ii][jj];
                  }
                  counter++;
               }
            }
            if(method == "STO-NLP" or method == "STO+STO-NLP"){
               for(int ii = 0; ii < T; ii++){
                  for(int jj = 0; jj < stateSize + controlSize; jj++){
                     solution[counter][jj] = nlp->solMatrix[ii][jj];
                  }
                  counter++;
               }
            }
         }else{
            infsblArc = k;
            cout << "Arc (" << route[k] << "," << route[k+1] << ") infeasible!" << endl;
            flightTimes[k] = DBL_MAX;
            stepSizes[k] = DBL_MAX;
            errors[k] = DBL_MAX;
            return IloAlgorithm::Infeasible;
         }
         //---------------------------------------------------------------------

         flightTimes[k] = tOpt;
         stepSizes[k] = stepOpt;
         errors[k] = errorOpt;
         
         for(int t = 0; t < T; t++){
            normEpsilons[k][t] = nlp->normErrors[t];
            normTaylor1st[k][t] = nlp->normTaylor1st[t];
            relEpsilons[k][t] = nlp->relativeErrors[t];
         }

         //update initial conditions for next leg
         Ybar[0] = solution[counter-1][0];
         Ybar[1] = solution[counter-1][1];
         Ybar[2] = solution[counter-1][2];
         Ybar[3] = solution[counter-1][3];
         Ybar[4] = solution[counter-1][4];
         Ybar[5] = solution[counter-1][5];
         Ubar[0] = solution[counter-1][6];
         Ubar[1] = solution[counter-1][7];
         
         double xf = solution[counter-1][0];
         double yf = solution[counter-1][1];
         double zf = solution[counter-1][2];
         double xi = solution[counter-2][0];
         double yi = solution[counter-2][1];
         double zi = solution[counter-2][2];
         vfAux[0] = (xf - xi)/stepOpt; 
         vfAux[1] = (yf - yi)/stepOpt; 
         vfAux[2] = (zf - zi)/stepOpt;
      }
   }else{
      cout << "The status of this route is infeasible!" << endl;
      return IloAlgorithm::Infeasible;
   };
   //##################################################
   //end of check route feasibility and compute optimal trajectories
   return status;
}


//IloAlgorithm::Status TrajOpt::findTraj(int* route, int seqSize, vector<vector<double>>& solution, double* flightTimes, double* stepSizes, double* errors)
//{
   ////##################################################
   ////check route feasibility and compute optimal trajectories
   //IloAlgorithm::Status status = IloAlgorithm::Infeasible;
   //const int T = data->T;
   //const int stateSize = dynPtr->stateSize;
   //const int controlSize = dynPtr->controlSize;
   //bool feasible = dynPtr->checkRouteStatus(route,seqSize);
   //if(feasible){
      ////initialisation
      //const Dynamics::dynamics leg0 = dynPtr->getLeg(route[0],route[1]);
      //double Ybar[stateSize] = {leg0.Yo[0], leg0.Yo[1], leg0.Yo[2],
                                //leg0.get_Yeq(3),leg0.get_Yeq(4),leg0.get_Yeq(5)};
      //double Ubar[controlSize] = {leg0.get_Ueq(0),leg0.get_Ueq(1)};      
      //double vfAux[3] = {0, 0, 0};
      
      //int counter = 0;

      ////main loop of Traj. Opt.
      //for(int k = 0; k < seqSize - 1; k++){
         ////get first first leg
         //const Dynamics::dynamics leg = dynPtr->getLeg(route[k],route[k+1]);
         //cout << "Origin: "  << leg.orig << endl;
         //cout << "Dest: "  << leg.dest << endl;
         ////double targetVel = 0;
         ////update boundary conditions
         ////states
         //double Yo_aux[] = {Ybar[0],Ybar[1],Ybar[2],Ybar[3],Ybar[4],Ybar[5]};
         //for(int ii = 0; ii < stateSize; ii++){
            //leg.Yo[ii] = Yo_aux[ii];
         //}
         //for(int ii = 3; ii < stateSize; ii++){
            //leg.Yf[ii] = 0;
         //}
         ////controls
         //leg.Uo[0] = Ubar[0]; 
         //leg.Uo[1] = Ubar[1];
         ////velocites
         //leg.vo[0] = vfAux[0];
         //leg.vo[1] = vfAux[1];
         //leg.vo[2] = vfAux[2];
         
         ////print boundary conditions
         //const double Yeq[] = {leg.Yo[0],leg.Yo[1],leg.Yo[2],
                               //leg.get_Yeq(3),leg.get_Yeq(4),leg.get_Yeq(5)};
         
         //cout << "\nBoundary conditions: " << endl;
         //printVec(leg.vo,3,"vo: ");
         //printVec(Yeq,stateSize,"Yeq: ");
         //printVec(leg.Yo,stateSize,"Yo: ");
         //printVec(leg.Uo,controlSize,"Uo: ");
         //printVec(leg.Yf,stateSize,"Yf: ");
         
         ////---------------------------------------------------------------------
         ////run the unconstrained optimization
         //Unconstrained* uncOpt = new Unconstrained();
         //double tfLb = leg.deltaTa;
         //double tfUb = leg.deltaTb;
         //double tfGuess = (tfUb + tfLb)/2;
         //double deltaT = tfUb - tfLb;
         //double tOpt = 0;
         //double stepOpt = 0;
         //double errorOpt;
         ////uncOpt->callCMAES(model,leg,tfGuess,tOpt);
         //status = IloAlgorithm::Infeasible;
         //while(status != IloAlgorithm::Optimal && tfUb < maxtb){
            ////status = uncOpt->bisection(model,leg,tfLb,tfUb,tOpt,stepOpt,errorOpt);
            //status = uncOpt->minFlightTime(model,leg,tfLb,tfUb,tOpt,stepOpt,errorOpt);
            ////getchar();
            //if(status == IloAlgorithm::Infeasible){
               //cout << "NOT OPTIMAL!" << endl;               
               ////exit(1);
            //}
            ////uncOpt->plotObjEval(data,model,leg,tfLb,tfUb);
            ////getchar();
            //tfLb += deltaT;
            //tfUb += deltaT;
         //}
         ////save solution
         //if(status == IloAlgorithm::Optimal){
            //for(int ii = 0; ii < T; ii++){
               //for(int jj = 0; jj < stateSize + controlSize; jj++){
                  //solution[counter][jj] = model->solMatrix[ii][jj];
                  ////cout << solutionFinal[counter][j] << " ";
               //}
               ////cout << endl;
               //counter++;
            //}
         //}else{
            ////cout << "Infeasible leg! (" << i <<", "<< j <<")"<< endl;
            //return IloAlgorithm::Infeasible;
         //}
         ////---------------------------------------------------------------------

         ////save flight time for this leg
         //flightTimes[k] = tOpt;
         //stepSizes[k] = stepOpt;
         //errors[k] = errorOpt;
         ////cout << "Flight time= " << flightTimes[k] << " | Step Size= " << stepSizes[k] << endl;
         ////getchar();

         ////update initial conditions for next leg
         //Ybar[0] = model->solMatrix[T-1][0];
         //Ybar[1] = model->solMatrix[T-1][1];
         //Ybar[2] = model->solMatrix[T-1][2];
         //Ybar[3] = model->solMatrix[T-1][3];
         //Ybar[4] = model->solMatrix[T-1][4];
         //Ybar[5] = model->solMatrix[T-1][5];
         //Ubar[0] = model->solMatrix[T-1][6];
         //Ubar[1] = model->solMatrix[T-1][7];
         
         //double xf = model->solMatrix[T-1][0];
         //double yf = model->solMatrix[T-1][1];
         //double zf = model->solMatrix[T-1][2];
         //double xi = model->solMatrix[T-2][0];
         //double yi = model->solMatrix[T-2][1];
         //double zi = model->solMatrix[T-2][2];
         //vfAux[0] = (xf - xi)/stepOpt; 
         //vfAux[1] = (yf - yi)/stepOpt; 
         //vfAux[2] = (zf - zi)/stepOpt;
         ////cout << "FINAL VEL VECTOR2= [" << vfAux[0] << "," << vfAux[1] << "," << vfAux[2] << "]" << endl;
         ////getchar();
      //}
   //}else{
      //cout << "The status of this route is infesible!" << endl;
      ////writeInfList(data);
      //return IloAlgorithm::Infeasible;
   //};
   ////##################################################
   ////end of check route feasibility and compute optimal trajectories
   //return status;
//}












