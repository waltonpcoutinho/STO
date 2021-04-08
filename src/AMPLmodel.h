#ifndef AMPLMODEL_H
#define AMPLMODEL_H

#include <ampl/ampl.h>
#include <string>

#include "Data.h"
#include "Dynamics.h"
#include "Util.h"

using namespace std;

class AMPLmodel{

   public:
      AMPLmodel(Data*,Dynamics*,int); 
      ~AMPLmodel();     

      double** solMatrix;
      double** errorMatrix;
      double* normErrors;
      double* normTaylor1st;
      double* relativeErrors;
      int dscrtSize;

      int solveNLP(Dynamics::dynamics, double&, double&, double&);

      void setLocalSolver(string);
      
   private:
      Data *dataPtr;
      Dynamics *gliderPtr;
      string localSolver;

      double errorNorm(ampl::AMPL* ampl);
      void setInitialGuess(ampl::AMPL* ampl);
      double* computeRelativeError(Dynamics::dynamics);

};

class MyOutputHandler : public ampl::OutputHandler{

   public:

      void output(ampl::output::Kind kind, const char* message)
      {
         if(kind == ampl::output::DISPLAY){
            setStatus(message);
         }
         if(kind == ampl::output::PROMPT){
            setStatus(message);
         }
         if(kind == ampl::output::SOLVE){
            setStatus(message);
         }
         if(kind == ampl::output::OPTION){
            setStatus(message);
         }
      }

      string getStatus()
      {
         return status;
      }

   private:
      string status;

      void setStatus(const char* message)
      {
         status = message;
      }
};

//class MyErrorHandler : public ampl::ErrorHandler{
//
   //public:
//
      //void error(const ampl::AMPLException& exception)
      //{
         //cout << exception.getMessage() << endl;
         //setStatus(exception.getMessage());
      //}
//
      //void warning(const ampl::AMPLException& exception)
      //{
         //cout << exception.getMessage() << endl;
         //setStatus(exception.getMessage());
      //}
      //
      //string getError()
      //{
         //if(!msg.empty()){
            //return msg;
         //}
      //}
//
   //private:
      //string msg;
//
      //void setStatus(const string& message)
      //{
         //if(!message.empty()){
            //msg = message;
         //}
      //}
//
//};

#endif
