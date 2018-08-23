/**
 * @file SOFCModel.h
 *
 *
 */

/*
 *$ author: Nikhil Galagali
 *
 */ 

#ifndef SOFCMODEL_H
#define SOFCMODEL_H

#include"SOFCAnode.h"
#include <cantera/Cantera.h>
#include <cantera/IdealGasMix.h>    // defines class IdealGasMix
#include <cantera/transport.h>      // transport properties
#include <cantera/importPhase.h>
#include <cantera/Interface.h>
#include <cantera/Metal.h>
#include <cantera/IncompressibleSolid.h>
#include <string>

#include "structDef.h"
#include "tools.h"

using namespace Cantera;
using namespace std;
using namespace Cantera_CXX;

class SOFCModel
{ 
 public:

 // Constructor
 SOFCModel();

 //Destructor 
  ~SOFCModel();
 
 /**
  *  Initialize the solid-oxide fuel cell model
  */
 void initialize();  
 

/**
 *  Method to calculate the cell potential 
 */
  double* calculateCellPotential(int,int,int,int,double*,double*);
             //!< Controls.
    
  SOFCAnode sofcAnode;

};

#endif
