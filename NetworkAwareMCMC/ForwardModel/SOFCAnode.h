/**
 * @file SOFCAnode.h
 *
 *
 */

/*
 *$ author: Nikhil Galagali
 *
 */ 

#include <cantera/Cantera.h>
#include <cantera/IdealGasMix.h>    // defines class IdealGasMix
#include <cantera/transport.h>      // transport properties
#include <cantera/kernel/DenseMatrix.h>
#include <cantera/numerics.h>
#include <cantera/integrators.h>
#include <cantera/kernel/Array.h>
#include <cantera/kernel/FuncEval.h>
#include <cantera/importPhase.h>
#include <cantera/Interface.h>
#include <cantera/Edge.h>
#include <cantera/Metal.h>
#include <cantera/IncompressibleSolid.h>
#include <string>
#include<math.h>
#include <time.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

#include "structDef.h"
#include "tools.h"

using namespace Cantera;
using namespace Cantera_CXX;
using namespace std;

#ifndef SOFCANODE_H
#define SOFCANODE_H

class SOFCAnode:public DenseMatrix
{
 public:

 // Constructor
 SOFCAnode();

 //Destructor 
 ~SOFCAnode();
 
 /**
  * Calculates the anode potential difference, Ec-Eelectrolyte
  */
 double* anodePotentialDifference(int,int,int,int,double*,double*);

 public:
 
 /* Controls. */
 Controls primary;                     //!< Controls.
 
};

#endif
