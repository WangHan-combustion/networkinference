#include "SOFCModel.h"

using namespace Cantera;
using namespace std;
using namespace Cantera_CXX;

SOFCModel::SOFCModel()
{

}

SOFCModel::~SOFCModel()
{
}

/**
 *  Method to calculate the cell potential 
 */
double* SOFCModel::calculateCellPotential(int modelProblem,int nReactions, int nSpecies, int nData, double* time,double* rateParameters)
{
  double* modelOutput; 
    
  // call the anode voltage vs. current solver
  modelOutput=sofcAnode.anodePotentialDifference(modelProblem,nReactions,nSpecies,nData,time,rateParameters);
  
  return modelOutput;
}
