#include "fileRead.h"
#include "structDef.h"
#include <string>

/*********************************************************************
 * Controls
 *********************************************************************/

void initializeControls(Controls &primary){

  /* General. */
  cout.precision(16);
  primary.nDim = 0;
  primary.quadMethod = 3;
  primary.dCompTol = 1.0e-14;
  primary.binFile = 0;

  /* PC. */
  primary.thetasLeftSup = NULL;
  primary.thetasRightSup = NULL;
  primary.pOrder = 5;
  primary.nTotalPCETerms = 1;
  primary.nInterFileWrites = -1;
  primary.interFileWriteNReuseNestedFcnEvals = NULL;

  /* Adaptive quadrature controls. */
  primary.reuseNested = 1;
  primary.maxNestedFcnEvalsStorage = 12500;
  primary.integralTol = 1.0e-14;
  primary.maxIndices = 100000;
  primary.maxAbscissasLevels = 23;
  primary.maxReuseNestedFcnEvals = 1000000;
  primary.maxNoReuseNestedFcnEvals = 1000000;
  primary.w = 1.0;
  primary.displayProgress = 0;
  primary.displaySummary = 1;
  primary.saveFinalSets = 0;
  primary.errorNorm = 1;

  /* Test function. */
  primary.nOutputs = 1;
  primary.testFcn = 1;
  primary.case0IntegrandShift = NULL;
  primary.case7Gamma = 0.5;
  primary.case8XSquares = 50;

  /*Fuel cell */
  primary.canteraInputFile="deutsc";
  primary.anodeGasName="Fuelcellgas";
  primary.cathodeGasName="Air";
  primary.fuelRatio="H2:0.2,H2O:0.8";
  primary.dataFile="DATAH2_20.inp";
  primary.dataPoint=1;
  primary.operatingPressure=101325.0;  // One atmospshere pressure
  primary.operatingTemperature=973.0;
  primary.faradaysConstant=96485.0;
  primary.bulkVacancyFraction=0.0;
  primary.nickelLineWidth=0.0;
  primary.nickelLineThickness=0.0;
  primary.YSZLineWidth=0.0;
  primary.YSZLineThickness=0.0;
  primary.lamdaTPB=0.0;
  primary.electrolyteConductivityPreexponential=0.0;
  primary.electrolyteConductivityActivation=0.0;
  primary.anodeSiteDensity=0.0;
  primary.electrolyteSiteDensity=0.0;
  primary.nElectrons=1;
  primary.anodeDiffusivities = new double*[6];
  primary.electrolyteDiffusivities = new double*[5];
  
 
  /* Reaction. */
  primary.pressure = 1.0;
  primary.tFinal = 2.0e0;
  primary.nLnPrefac = 0;
  primary.rxnLnPrefac = NULL;
  primary.nEa = 0;
  primary.rxnEa = NULL;
  primary.errStdDevTauIgnPeakHeight = 0.75;
  primary.enthalpyPeakTol = 1.0e5;
  primary.enthalpySettleTol = 1000.0;
  primary.nSettleEquil = 10;
  primary.writeFirstReactionResults = 0;
  primary.useAir = 0;

  /* Read any Controls input from file. */
  readFileControls(primary);

}

void addMPIInfoControls(Controls &primary, const int nTasks, const int rank){


}

void finalizeControls(Controls &primary){

  /* Free memory. */
  delete [] primary.thetasLeftSup;
  delete [] primary.thetasRightSup;
  delete [] primary.case0IntegrandShift;
  delete [] primary.rxnLnPrefac;
  delete [] primary.rxnEa;
  if (primary.nInterFileWrites>0)
    delete [] primary.interFileWriteNReuseNestedFcnEvals;
  
}
