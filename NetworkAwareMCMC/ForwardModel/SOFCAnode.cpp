#include "SOFCAnode.h"
#include "SOFCAnodeIntegrator.h"


const double F=96485.0;
const double pi=22.0/7.0;
using namespace Cantera;
using namespace std;
using namespace Cantera_CXX;

SOFCAnode::SOFCAnode()
{
 
}

SOFCAnode::~SOFCAnode()
{
}

/**
 * Calculates the anode potential difference, Ea-Eelectrolyte
 */
double* SOFCAnode::anodePotentialDifference(int modelProblem,int nReactions,int nSpecies,int nData, double* time,double* rateParameters)
 {

  double* modelOutput;  
    
  SOFCAnodeIntegrator sofcAnodeIntegrator(modelProblem,nReactions,nSpecies,rateParameters);// integration object  
     
  modelOutput=sofcAnodeIntegrator.advance(nData,time);
     
  return modelOutput;
}



