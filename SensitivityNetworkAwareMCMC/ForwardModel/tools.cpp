#include "tools.h"

/* Functions with templates are defined in header file. */

int nCr(const int n, const int r){

  /* Initializations. */
  int largerDenom = max<int>(r,n-r);
  int smallerDenom = n-largerDenom;
  int result(1);

  /* Cancel the larger of the denominator terms. */
  for (int i=largerDenom; i<n; i++){
    result *= i+1;
    if (result<0)
      return -1;
  }

  /* Final result. */
  result /= factorial(smallerDenom);
  
  return result;

}

int factorial(const int n){
  
  /* Initialization. */
  int prod(1);
  
  for (int i=0; i<n; i++){
    prod *= i+1;
  }
  
  return prod;
}

int binarySearch(const double query, const double* const refList, 
		 const int nLength, const double dCompTol){

  /* Initializations. */
  int first(0), last(nLength-1);

  while (first <= last) {
    int mid = (first+last)/2;  /* Compute mid-point. */
    if (query-refList[mid]>dCompTol)
      first = mid+1;    /* Repeat search in top half. */
    else if (query-refList[mid]<-dCompTol) 
      last = mid-1;     /* Repeat search in bottom half. */
    else
      return mid;       /* Found it, return position. */
  }
  
  return -1;            /* Failed to find key. */
  
}

void recursiveLegendre(const double x, const int maxDegree, 
		       double* const results){
  
  /* Compute term via 3-term recursive formula. */
  for (int i=0; i<maxDegree+1; i++){
    if (i==0){
      results[0] = 1.0;
    }else if (i==1){
      results[1]= x;
    }else{
      results[i] = ((2.0*double(i-1)+1.0)*x*results[i-1]-double(i-1)*results[i-2])/double(i);
    }
  }

}
