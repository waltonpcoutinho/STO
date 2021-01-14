#include "AMPLmodel.h"

AMPLmodel::AMPLmodel(Data* data, Dynamics* glider, int T)
   :dataPtr(data), gliderPtr(glider)
{
   //Set discretisation size
   dscrtSize = T;
   //allocate memory to the solution matrix
   int state = gliderPtr->stateSize;
   int ctrl = gliderPtr->controlSize;
   solMatrix = new double* [dscrtSize];
   for(int t = 0; t < dscrtSize; t++){
      solMatrix[t] = new double [state+ctrl];
   }
   for(int t = 0; t < dscrtSize; t++){
      for(int j = 0; j < state+ctrl; j++){
         solMatrix[t][j] = 0.0;
      }
   }

   //default local solver
   localSolver = "worhp_ampl";
}	

AMPLmodel::~AMPLmodel()
{
   //comment
}
void AMPLmodel::setLocalSolver(string solver)
{
   localSolver = solver;
}

double AMPLmodel::errorNorm(ampl::AMPL* ampl)
{
   double sum = 0;
   
   int state = gliderPtr->stateSize;
   int ctrl = gliderPtr->controlSize;

   double** errorMatrix;
   errorMatrix = new double* [dscrtSize];
   for(int t = 0; t < dscrtSize; t++){
      errorMatrix[t] = new double [state];
   }

   // Get the values of the error variables
   ampl::DataFrame errors = ampl->getVariable("epsilon").getValues();
   int rowIdx = 0;
   for(int t = 0; t < errors.getNumRows(); t++){
      ampl::DataFrame::Row row = errors.getRowByIndex(t); 
      int col = row[1].dbl() - 1; //AMPL cols start by 1
      errorMatrix[rowIdx][col] = row[2].dbl();
      if(col == 5){
         rowIdx++;
      }
   }

   for(int t = 0; t < dscrtSize; t++){
      double sumSquare = 0;
      for(int i = 0; i < state; i++){
         sumSquare += errorMatrix[t][i]*errorMatrix[t][i];
      }
      sum += sqrt(sumSquare);
   }

   cout << "Error norm = " << sum << endl;

   return sum;
}

int AMPLmodel::solveNLP(Dynamics::dynamics arc, double& flightTime, double& step, double& error)
{
   //state lower and upper bounds
   const double Ylb[] = {dataPtr->xLb, dataPtr->yLb, dataPtr->zLb, gliderPtr->vLb, gliderPtr->gammaLb, gliderPtr->phiLb};
   const double Yub[] = {dataPtr->xUb, dataPtr->yUb, dataPtr->zUb, gliderPtr->vUb, gliderPtr->gammaUb, gliderPtr->phiUb};

   //control lower and upper bound
   const double Ulb[] = {gliderPtr->CLlb, gliderPtr->muLb};
   const double Uub[] = {gliderPtr->CLub, gliderPtr->muUb};
   
   //norms of derivatives and matrices
   const double normuDotUb = gliderPtr->uDotUbMaxNorm;
   const double normyDotUb = gliderPtr->yDotUbMaxNorm;
   const double normA = arc.getMaxNormA();
   const double normB = arc.getMaxNormB();
   
   //make a copy of the system's data
   const int stateSize = gliderPtr->stateSize;
   const int controlSize = gliderPtr->controlSize;
   
   //initial and final conditions
   double const* Yo = arc.Yo;
   double const* Yf = arc.Yf;
   double const* Uo = arc.Uo;
   double const* vo = arc.vo;
   int const orig = arc.orig;
   
   //update x,y,z components of equilibrium conditions
   const double Yeq[] = {Yo[0],Yo[1],Yo[2],arc.get_Yeq(3),arc.get_Yeq(4),arc.get_Yeq(5)};   
   const double Ueq[] = {arc.get_Ueq(0),arc.get_Ueq(1)};   

   /*
    * Start AMPL API
    * */

   int status = 0;

   ampl::AMPL ampl;

   //make sure it is reset
   ampl.eval("reset;");

   // Interpret the two files
   ampl.read("src/AMPLfiles/grtop.mod");
   ampl.readData("src/AMPLfiles/grtop.dat");

   // Set discretisation size
   ampl::Parameter N = ampl.getParameter("N");
   N.set(dscrtSize);

   /*
    * setting the values of the parameters
    */

   //set the values of A
   ampl::Parameter A = ampl.getParameter("A");

   double stateVals[36];
   double stateRows[] = {1,2,3,4,5,6};
   double stateCols[] = {1,2,3,4,5,6};

   int count = 0;
   for(int i = 0; i < stateSize; i++){
      for(int j = 0; j < stateSize; j++){
         stateVals[count] = arc.get_A(i,j);
         count++;
      }
   }

   A.setValues(6, stateRows, 6, stateCols, stateVals, false);
   
   //set the values of B
   ampl::Parameter B = ampl.getParameter("B");

   double controlVals[12];
   double controlRows[] = {1,2,3,4,5,6};
   double controlCols[] = {1,2};

   count = 0;
   for(int i = 0; i < stateSize; i++){
      for(int j = 0; j < controlSize; j++){
         controlVals[count] = arc.get_B(i,j);
         count++;
      }
   }

   B.setValues(6, controlRows, 2, controlCols, controlVals, false);
   
   //set the values of Yeq, Ueq
   ampl::Parameter AMPLYeq = ampl.getParameter("Yeq");
   ampl::Parameter AMPLUeq = ampl.getParameter("Ueq");

   AMPLYeq.setValues(Yeq, 6);
   AMPLUeq.setValues(Ueq, 2);

   //set the values of Ylb, Yub, Ulb, Uub
   ampl::Parameter AMPLYlb = ampl.getParameter("Ylb");
   ampl::Parameter AMPLYub = ampl.getParameter("Yub");

   AMPLYlb.setValues(Ylb, 6);
   AMPLYub.setValues(Yub, 6);
   
   ampl::Parameter AMPLUlb = ampl.getParameter("Ulb");
   ampl::Parameter AMPLUub = ampl.getParameter("Uub");

   AMPLUlb.setValues(Ulb, 2);
   AMPLUub.setValues(Uub, 2);
   
   //set the values of the upper bunds on norms
   ampl::Parameter AMPLnorm_yDotUb = ampl.getParameter("norm_yDotub");
   ampl::Parameter AMPLnorm_uDotUb = ampl.getParameter("norm_uDotub");
   AMPLnorm_yDotUb.set(normyDotUb);
   AMPLnorm_uDotUb.set(normuDotUb);
   
   ampl::Parameter AMPLnormA = ampl.getParameter("normA");
   ampl::Parameter AMPLnormB = ampl.getParameter("normB");
   AMPLnormA.set(normA);
   AMPLnormB.set(normB);
   
   //set the values of Y0, U0
   ampl::Parameter AMPLY0 = ampl.getParameter("Y0");
   ampl::Parameter AMPLU0 = ampl.getParameter("U0");

   AMPLY0.setValues(Yo, 6);
   AMPLU0.setValues(Uo, 2);

   //set the values of scalar params
   ampl::Parameter UbarH = ampl.getParameter("UbarH");
   ampl::Parameter barH = ampl.getParameter("barH");
   ampl::Parameter hatGamma = ampl.getParameter("hatGamma");
   ampl::Parameter hatMu = ampl.getParameter("hatMu");

   ampl::Parameter xBar = ampl.getParameter("xBar");
   ampl::Parameter yBar = ampl.getParameter("yBar");
   ampl::Parameter rBar = ampl.getParameter("rBar");

   ampl::Parameter xTilde = ampl.getParameter("xTilde");
   ampl::Parameter yTilde = ampl.getParameter("yTilde");
   ampl::Parameter rTilde = ampl.getParameter("rTilde");

   if(!arc.isLanding){
      //add visiting constraints
      ampl.read("src/AMPLfiles/visiting.mod");
      UbarH.set(arc.fzMin);
      barH.set(arc.fzMax);
      hatGamma.set(0.09);
      hatMu.set(0.09);
      xBar.set(Yf[0]);
      yBar.set(Yf[1]);
      rBar.set(arc.fRadius);
   }else{
      //add landing constraints
      ampl.read("src/AMPLfiles/landing.mod");
      xTilde.set(Yf[0]);
      yTilde.set(Yf[1]);
      rTilde.set(arc.fRadius);
   }
   
   //setting the value of \tau_f for avoiding 
   //solutions with \tau_f = 0
   ampl::Parameter AMPLtflb = ampl.getParameter("tflb");
   cout << "minimum flight time = " << arc.deltaTa << endl;
   AMPLtflb.set(arc.deltaTa);
   
   /*
    * end of setting the values of the parameters
    */

   //export model
   ampl.eval("expand >Results/NLPmodel.lp;");
   ampl.eval("write bResults/NLPmodel;");
   
   /*
    * Call optmimiser
    */

   // Set the output handler to accumulate the output messages
   MyOutputHandler output;
   ampl.setOutputHandler(&output);
   
   //set solver
   if(localSolver == "worhp_ampl"){
      ampl.setOption("solver", "worhp_ampl");
   }else if(localSolver == "ipopt"){
      ampl.setOption("solver", "ipopt");
      //ampl.setOption("ipopt_options", "iprint=1");
   }

   // Solve
   ampl.solve();

   //Display solver's message
   ampl.eval("display solve_message;");
   string message = output.getStatus();
   if(localSolver == "worhp_ampl"){
      cout << "WORHP-log:" << endl;
      cout << message.substr(0,string::npos);
      cout << flush;
   }else if(localSolver == "ipopt"){
      cout << "IPOPT-log:" << endl;
      cout << message.substr(19,52);
      cout << flush;
   }

   //Trigger the output handler
   ampl.eval("display solve_result;");
   //Get output message
   if(output.getStatus().find("solved") != string::npos){
      status = 1;
      //cout << "Optimal Solution found by NLP solver" << endl;
      
      //Get solution value
      //error = ampl.getValue("epsilon").dbl();
      error = errorNorm(&ampl);
      flightTime = ampl.getValue("tf").dbl();
      step = ampl.getValue("h").dbl();

      cout << "optimal flight time = " << flightTime << endl;
      cout << "optimal step size = " << step << endl;
      cout << endl;

      // Get the values of the variables Y and U in dataframes
      ampl::DataFrame state = ampl.getVariable("Y").getValues();
      int rowIdx = 0;
      for(int t = 0; t < state.getNumRows(); t++){
         ampl::DataFrame::Row row = state.getRowByIndex(t); 
         int col = row[1].dbl() - 1; //AMPL cols start by 1
         solMatrix[rowIdx][col] = row[2].dbl();
         if(col == 5){
            rowIdx++;
         }
      }
      
      ampl::DataFrame control = ampl.getVariable("U").getValues();
      rowIdx = 0;
      for(int t = 0; t < control.getNumRows(); t++){
         ampl::DataFrame::Row row = control.getRowByIndex(t); 
         int col = row[1].dbl() + 5; //AMPL cols start by 1
         solMatrix[rowIdx][col] = row[2].dbl();
         if(col == 7){
            rowIdx++;
         }
      } 
   }else{
      cout << "Failed to solve NLP!" << endl;
   }

   //reset and close AMPL env
   ampl.eval("reset;");
   ampl.close();

   return status;
}
