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
      int dscrtSize;

      int solveNLP(Dynamics::dynamics, double&, double&, double&);

      void setLocalSolver(string);
      
   private:
      Data *dataPtr;
      Dynamics *gliderPtr;
      string localSolver;

};

class MyOutputHandler : public ampl::OutputHandler{

   public:

      void output(ampl::output::Kind kind, const char* message){
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

      string getStatus(){
         return status;
      }

   private:
      string status;

      void setStatus(const char* message){
         status = message;
      }
};
#endif
