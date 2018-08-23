/**
 * @file main.cpp
 *
 *
 */

/*
 *$ author: Nikhil Galagali
 *
 */ 

/********************************************************************
This code is used as an odeintegrator for example 1
*********************************************************************/

#include "mex.h"
#include "main.h"
#include "SOFCModel.h"

#include <string>
#include <math.h>
#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;

int main(int modelProblem,int nReactions,int nSpecies,int nData, double* time,double* rateParameters,double* output1)
{ 
          
 double* output;   
    
  SOFCModel sofcModel;
  output=sofcModel.calculateCellPotential(modelProblem,nReactions,nSpecies,nData,time,rateParameters);
  
  for (int j=0;j<nData;j++)
  {
   for (int i=0;i<nReactions+nSpecies;i++)
       output1[j*(nReactions+nSpecies)+i]=output[j*(nReactions+nSpecies)+i];
  }
  
  delete [] output;
  return 0;
  
}

void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int a;


  
  int modelProblem;
  int nReactions;
  int nSpecies;
  int nData;
  double *multiplier;
  double *inMatrix;       /* 1xN input matrix */
  double *outMatrix;      /* output matrix */

  
  /* get the value of the scalar input */
  modelProblem = mxGetScalar(prhs[0]);
  
  /* get the value of the scalar input */
  nReactions = mxGetScalar(prhs[1]);
  
  /* get the value of the second scalar input */
  nSpecies = mxGetScalar(prhs[2]);
  
  /* get the value of the third scalar input */
  nData = mxGetScalar(prhs[3]);
  
  /* create a pointer to the real data in the input matrix  */
  multiplier = mxGetPr(prhs[4]);
    
  /* create a pointer to the real data in the input matrix  */
  inMatrix = mxGetPr(prhs[5]);

  /* initialize ncols */
  mwSize nCols=nData*(nReactions+nSpecies); 
  
  /* create the output matrix */
  plhs[0] = mxCreateDoubleMatrix(1,nCols,mxREAL);

  /* get a pointer to the real data in the output matrix */
  outMatrix = mxGetPr(plhs[0]);
  
  a=main(modelProblem,nReactions,nSpecies,nData,multiplier,inMatrix,outMatrix);
    
}