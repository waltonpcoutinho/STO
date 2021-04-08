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
   //allocate memory for the error Matrix
   errorMatrix = new double* [dscrtSize];
   for(int t = 0; t < dscrtSize; t++){
      errorMatrix[t] = new double [state];
   }
   for(int t = 0; t < dscrtSize; t++){
      for(int j = 0; j < state; j++){
         errorMatrix[t][j] = 0.0;
      }
   }

   normErrors = new double [dscrtSize];
   for(int t = 0; t < dscrtSize; t++){
      normErrors[t] = 0.0;
   }
   normTaylor1st = new double [dscrtSize];
   for(int t = 0; t < dscrtSize; t++){
      normTaylor1st[t] = 0.0;
   }
   relativeErrors = new double [dscrtSize];
   for(int t = 0; t < dscrtSize; t++){
      relativeErrors[t] = 0.0;
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
   double max = 0;
   
   int state = gliderPtr->stateSize;
   int ctrl = gliderPtr->controlSize;

   // Get the values of the error variables
   ampl::DataFrame errors = ampl->getVariable("epsilon").getValues();

   //Print errors as string
   //cout << errors.toString() << endl;
   //cout << ampl->getOutput("display epsilon;");
   
   for(int t = 0; t < errors.getNumRows(); t++){
      int row = errors.getRowByIndex(t)[0].dbl();
      int col = errors.getRowByIndex(t)[1].dbl() - 1; //AMPL cols start at 1
      double val = errors.getRowByIndex(t)[2].dbl();
      errorMatrix[row][col] = val;
      if(val > max){
         max = val;
      }
   }

   cout << "Error norm = " << max << endl;

   return max;
}

double vectorMaxNorm(vector<double> vec)
{
   double max = 0;
   for(int i = 0; i < vec.size(); i++){
      if(abs(vec[i]) > max){
         max = abs(vec[i]);
      }
   }
   return max;
}

vector<double> vectorMatrixMultiplication(vector<vector<double>> matrix, vector<double> vec)
{
   int numRows = matrix.size();
   int numCols = matrix[0].size();

   vector<double> result(numRows, 0);

   if(numCols != vec.size()){
      cerr << "The size of arrays do not match for multiplication!" << endl;
      cerr << "FUNCTION vectorMatrixMultiplication (AMPLmodel.cpp)" << endl;
      exit(1);
   }

   for(int i = 0; i < numRows; i++){
      double rowSum = 0;
      for(int j = 0; j < numCols; j++){
         rowSum += matrix[i][j]*vec[j];
      }
      result[i] = rowSum;
   }

   return result;
}

vector<double> subtractVectors(vector<double> vec1, vector<double> vec2)
{
   vector<double> res(vec1.size(), 0);
   
   if(vec1.size() != vec2.size()){
      cerr << "It's not possible to subtract these two vectors" << endl;
      cerr << "FUNCTION subtractVectors (AMPLmodel.cpp)" << endl;
      exit(1);
   }
   
   transform (vec1.begin(), vec1.end(), vec2.begin(), 
                          res.begin(), std::minus<double>());

   return res;
}

double* AMPLmodel::computeRelativeError(Dynamics::dynamics arc)
{
   double* relErrors;
   relErrors = new double [dscrtSize];
   for(int t = 0; t < dscrtSize; t++){
      relErrors[t] = 0.0;
   }

   int state = gliderPtr->stateSize;
   int ctrl = gliderPtr->controlSize;

   //initial and final conditions
   double const* Yo = arc.Yo;
   double const* Uo = arc.Uo;
   
   //set up equilibrium conditions
   vector<double> Yeq {Yo[0],Yo[1],Yo[2],
                       arc.get_Yeq(3),arc.get_Yeq(4),arc.get_Yeq(5)};   
   vector<double> Ueq {arc.get_Ueq(0),arc.get_Ueq(1)};
   
   //set up state and control matrices
   vector<vector<double>> A(state, vector<double>(state, 0));
   vector<vector<double>> B(state, vector<double>(ctrl, 0));

   for(int i = 0; i < state; i++){
      for(int j = 0; j < state; j++){
         A[i][j] = arc.get_A(i,j);
      }
   }

   for(int i = 0; i < state; i++){
      for(int j = 0; j < ctrl; j++){
         B[i][j] = arc.get_B(i,j);
      }
   }

   //compute the relative error (of the first term of Taylor's expansion
   //relative to the magnitude of the error variable)
   for(int t = 0; t < dscrtSize; t++)
   {
      //declaring auxiliary variables for state, control and error values
      vector<double> y_t(state, 0);
      vector<double> u_t(ctrl, 0);
      vector<double> epsilon_t(state, 0);

      for(int s = 0; s < state; s++){
         y_t[s] = solMatrix[t][s];
      }      
      for(int c = 0; c < ctrl; c++){
         u_t[c] = solMatrix[t][c + 6];
      }
      for(int s = 0; s < state; s++){
         epsilon_t[s] = errorMatrix[t][s];
      }

      //computing y_t - y_eq
      vector<double> y_diff = subtractVectors(y_t, Yeq);
      vector<double> u_diff = subtractVectors(u_t, Ueq);

      //computing first term of Taylor's expansion
      vector<double> taylor1st_t(state, 0);
      vector<double> AtimesY = vectorMatrixMultiplication(A, y_diff);
      vector<double> BtimesU = vectorMatrixMultiplication(B, u_diff);

      for(int s = 0; s < state; s++){
        taylor1st_t[s] = AtimesY[s] + BtimesU[s];
      }

      //get the norms of numerator and denominator
      normErrors[t] = vectorMaxNorm(epsilon_t);
      normTaylor1st[t] = vectorMaxNorm(taylor1st_t);

      //compute relative errors
      if(normTaylor1st[t] > 1e-3){
         relErrors[t] = abs(1 - normErrors[t]/normTaylor1st[t]);
      }else{
         relErrors[t] = abs(normErrors[t]);
      }

   }

   return relErrors;
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
    * Set initial guesses
    */
   setInitialGuess(&ampl);
   
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
   cout << "##" << output.getStatus() << endl;
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

      //compute the quality of solutions as the relative error
      //between epsilon and 1st term of Taylor's expansion
      relativeErrors = computeRelativeError(arc);

   }else{
      cout << "\n\n Failed to solve NLP!, line 475 AMPLmodel.cpp \n\n" << endl;
   }

   //reset and close AMPL env
   ampl.eval("reset;");
   ampl.close();

   return status;
}
   
void AMPLmodel::setInitialGuess(ampl::AMPL* ampl)
{
   int stateSize = gliderPtr->stateSize;
   int ctrlSize = gliderPtr->controlSize;

   // Set up dataframes and copy STO solutions
   ampl::DataFrame state = ampl->getVariable("Y").getValues();
   int counter = 0;
   for(int t = 0; t < dscrtSize; t++){
      for(int s = 0; s < stateSize; s++){
         // skip initial configuration
         if(t > 0){
            state.setValue(counter, 2, solMatrix[t][s]);
         }
         counter += 1;
      }
   }
   
   ampl::DataFrame control = ampl->getVariable("U").getValues();
   counter = 0;
   for(int t = 0; t < dscrtSize; t++){
      for(int c = 0; c < ctrlSize; c++){
         // skip initial configuration
         if(t > 0){
            control.setValue(counter, 2, solMatrix[t][6+c]);
         }
         counter += 1;
      }
   }

   ampl->getVariable("Y").setValues(state);
   ampl->getVariable("U").setValues(control);
}
