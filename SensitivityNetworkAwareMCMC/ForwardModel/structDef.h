/*! \file structDef.h 

  \brief Definition of structures.
    
  Controls and Quad1D structures.
*/
#ifndef _STRUCTDEF_H
#define _STRUCTDEF_H

#undef BAND

#include <string>
#include <iostream>
//#include "mpi.h"

using namespace std;

/*! \struct Controls

  \brief Struct for control variables.  

  Note that only shallow copies are used/needed since a single
  reference struct is constructed. No alterations to members should be
  made for the reference or copies after reading the input file.
*/
struct Controls{

  /* General. */
  int nDim;                     //!< Number of dimensions.
  int quadMethod;               //!< 1-GL, 2-GP, 3-CC.
  double dCompTol;              //!< Compare precision for doubles.
  bool binFile;                 //!< Write/read binary file format.

  
  /* PC. */
  double* thetasLeftSup;        //!< Lower supports of the input parameters.
  double* thetasRightSup;       //!< Upper supports of the input parameters.
  int pOrder;                   //!< Global order desired for the PCE.
  int nTotalPCETerms;           //!< Total number of PCE polynomial terms.
  int nInterFileWrites;         //!< Number of intermediate PCE coef file writes (not including final).
  /*! \brief Pointer to the nReuseNestedFcnEvals positions which to
    perform the above coef and index sets file writes. 

    Don't choose these numbers too close together, otherwise some of
    them would be skipped if a single iteration spans multiple of
    these checkpoints. Also, the file names would store these numbers
    instead of the actual nReuseNestedFcnEvals at those points, for an
    easier time at loading these files in MATLAB; but this can be
    changed if desired. Also, nNoReuseNestedFcnEvals or
    nActualFcnEvals can be used to set these checkpoints too if
    desired, but code needs to be (slightly) modified.
  */
  int* interFileWriteNReuseNestedFcnEvals; 

  /* Dimension adaptive. */
  bool reuseNested;             //!< Flag to reuse nested quadrature function values. If true, up to maxNestedFcnEvalsStorage storages can be made.
  /*! \brief Max number of nested function evaluation
    storage. Additional new evaluations are not stored and must be
    computed each time.
    
    This number is estimated based on available memory for each
    CPU. Don't forget to divide by the extra factor of nOutputs when
    computing this. For example, on a 4 GB node with 2 CPUs, give each
    CPU 2 GB. There are nTotalPCETerms=1820 integrals that are
    integrated simultaneously, and from the setup right now, would
    need to store all evaluations for each of them. In the future,
    this needs to be improved by only store the G(\xi) values, instead
    of G(\xi)*\Psi(\xi), reducing the factor of nTotalPCETerms. Give
    or take 2*1.024*10e9/1820(coefficients)/8(bytes per
    double)/10(nOutputs) = 14065 (take 12500 to be safe; this has been
    tested to be working in this scenario). Also note that overhead
    may vary depending on other factors such as the PC order (i.e.,
    number of PCE terms).
  */
  int maxNestedFcnEvalsStorage; //!< Max number of stored function values for reuseNested.
  double integralTol;           //!< Adaptive quadrature error estimate tolerance.
  int maxIndices;               //!< Max number of indices in adaptation.
  int maxAbscissasLevels;       //!< Max number of quadrature levels in any dimension.
  int maxReuseNestedFcnEvals;   //!< Max number of function evaluations as if reuseNested were true.
  int maxNoReuseNestedFcnEvals; //!< Max number of function evaluations as if reuseNested were false.
  double w;                     //!< Weight for computing error indicator.
  bool displayProgress;         //!< Flag to display progress.
  bool displaySummary;          //!< Flag to display final summary.
  bool saveFinalSets;           //!< Flag to save final old and active sets.
  int errorNorm;                //!< Norm selection for error indicator.

  /* Test function. */
  int nOutputs;                 //!< Dimension of the function output.
  int testFcn;                  //!< Test function selection.
  double* case0IntegrandShift;  //!< Shift to integrand of test function 0.
  double case7Gamma;            //!< Parameter gamma in test function 7.
  int case8XSquares;            //!< Number of x^2 terms in test function 8.

  /* Fuel cell */
  string canteraInputFile;
  string anodeGasName;
  string cathodeGasName;
  string fuelRatio;
  string dataFile;
  int dataPoint;
  double operatingPressure;
  double operatingTemperature;
  double faradaysConstant;
  double bulkVacancyFraction;
  double nickelLineWidth;
  double nickelLineThickness;
  double YSZLineWidth;
  double YSZLineThickness;
  double lamdaTPB;
  double electrolyteConductivityPreexponential;
  double electrolyteConductivityActivation;
  double anodeSiteDensity;
  double electrolyteSiteDensity;
  int nElectrons;
  double** anodeDiffusivities;
  double** electrolyteDiffusivities; 

  /* Reaction. */
  double pressure;              //!< Initial pressure.
  double tFinal;                //!< Final time.
  int nLnPrefac;                //!< Number of lnPrefac parameters.
  int* rxnLnPrefac;             //!< List of reaction numbers of lnPrefac.
  int nEa;                      //!< Number of Ea parameters.
  int* rxnEa;                   //!< List of reaction numbers for Ea.
  double errStdDevTauIgnPeakHeight; 
                                //!< % of peak height for computing error standard deviation of ignition delay.
  double enthalpyPeakTol;       //!< Peaking tolerance value of enthalpy release rate.
  double enthalpySettleTol;     //!< Settling tolerance value of enthalpy release rate.
  int nSettleEquil;             //!< Number of settled time steps before declaring equilibrium.
  bool writeFirstReactionResults;
                                //!< Flag to indicate writing the 1st reaction results to file.
  bool useAir;                  //!< Flag to use air (if not, use oxygen).

};

/*! \fn void initializeControls(Controls&);

  \brief Initializes the values for the primary controls, by first
  setting the default values, and then make any corrections by reading
  in from file (see fileRead.h).
 
  \param primary Reference to Controls structure of interest.
*/
void initializeControls(Controls&);

/*! \fn void addMPIInfoControls(Controls&, const int, const int);

  \brief Stores the MPI variables in the Controls structure.
  
  \param primary Reference to controls structure of interest.
  \param nTasks Total number of CPUs.
  \param rank Current CPU rank.
*/
void addMPIInfoControls(Controls&, const int, const int);

/*! \fn void finalizeControls(Controls&);

  \brief Frees memory in the controls structure.
  
  \param primary Reference to Controls structure of interest.
*/
void finalizeControls(Controls&);

/*! \struct Quad1D

  \brief Struct for storing a single 1D quadrature rules.

  This structure is used to store the constructed 1D quadrature
  rule. The program first needs to pass the necessary information
  into this structure, and then use this to request the construction
  of the 1D quadrature rule.
*/
struct Quad1D{

  /* Parameters. */
  int nAbscissas;       //!< Number of abscissas for this 1D rule.
  int quadMethod;       //!< 1-GL, 2-GP, 3-CC.
  double dCompTol;      //!< Comparison precision for doubles.

  /* Pointers for storage (need to allocate memory separately). */
  double* abscissas;    //!< Pointer to storage for abscissas.
  double* weights;      //!< Pointer to storage for weights.

};

#endif
