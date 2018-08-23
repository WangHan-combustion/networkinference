/*! \file main.h 

  \brief A program to construct the polynomial chaos expansions
  non-intrusively.
  
  This program computes the polynomial chaos expansions of the
  observables of the specified mechanism modelled through
  Cantera. There are a total of 4 parameters, 2 kinetic parameters and
  2 design variables. This is a non-intrusive construction, using
  Galerkin projection that is computed using (dimension-adaptive)
  sparse quadrature. This code may be parallelized.  
*/
#ifndef _MAIN_H
#define _MAIN_H

#undef BAND

#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//#include "adaptiveQuadrature.h"
//#include "SOFC.h"
//#include "PCExpansion.h"
#include "structDef.h"

using namespace std;

#endif
