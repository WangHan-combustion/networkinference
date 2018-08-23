/*! \file tools.h 

  \brief Various tools.
*/
#ifndef _TOOLS_H
#define _TOOLS_H

#include <iostream>
#include <math.h>
#include <sstream>

using namespace std;

/*! \fn int nCr(const int, const int);

  \brief Computes the combination term, n choose r.
  
  \param n Total elements.
  \param r Number of elements for combination.
  
  \return The computed combination. If this <0, an overflow has
  occurred.
*/
int nCr(const int, const int);

/*! \fn int factorial(const int);

  \brief Computes the factorial of an integer.
  
  \param n Factorial argument.
  
  \return The computed factorial. If this <0, an overflow has
  occurred.
*/
int factorial(const int);

/*! \fn int binarySearch(const double, const double* const, const int,
  const double);
  
  \biref Performs a binary search of a query value in a sorted 1D
  array with distinct entries. The corresponding index value is
  returned.
  
  \param query Search query.
  \param refList Sorted 1D array with distinct entries.
  \param nLength Length of 1D array.
  \param dCompTol Comparison tolerance for double type.
  
  \return Returns the index of the query. If not found, returns -1.
*/
int binarySearch(const double, const double* const, const int, const double);

/*! \fn void recursiveLegendre(const double, const int, double*
  const);
  
  \brief Computes the Legendre polynomial values evaluated for some x
  via the 3-term recursive formula.
  
  \param x Point of evaluation for the Legendre polynomials.
  \param maxDegree Maximum order of the Legendre polynomials.
  \param results Storage for the results.
*/
void recursiveLegendre(const double, const int, double* const);

/*! \fn void quickSort(T* const, U* const, V** const, const int, const
  int);

  \brief Sorts a 1D array of elements in ascending order using a
  recursive implementation of quick sort. It can also sort a secondary
  1D array and a 2D array according to the sorting of the 1st 1D
  array.
  
  \param arrayA A 1D array to be sorted.
  \param arrayB Another 1D array to be sorted according to the sorting
  of arrayA.
  \param arrayC A 2D array to be sorted according to the sorting of
  arrayA.
  \param left Starting position to the section in arrayA to be sorted.
  \param right Ending position to the section in arrayA to be sorted.
*/
template <class T, class U, class V>
void quickSort(T* const, U* const, V** const, const int, const int);

/*! \fn int partition(T* const, U* const, V** const, const int, const
  int, const int);
  
  \brief The bulk of the quick sort. It sorts the pivot, and then
  partitions the array to be sorted into 2 subarrays from the pivot
  (which would be at the correct global position after the previous
  iteration).
  
  \param arrayA A 1D array to be sorted.
  \param arrayB Another 1D array to be sorted according to the sorting
  of arrayA.
  \param arrayC A 2D array to be sorted according to the sorting of
  arrayA.
  \param left Starting position to the section in arrayA to be sorted.
  \param right Ending position to the section in arrayA to be sorted.
  \param pivotIndex The element index selected as the pivot.
  
  \return The final (correct global) element index of the selected
  pivot.
*/
template <class T, class U, class V>
int partition(T* const, U* const, V** const, const int, const int, const int);

/*! \fn void maxOf1DArray(const T* const, const int, int&, T&);
  
  \brief Finds the maximum index and value in a 1D array.
  
  \param V The 1D array of interest.
  \param length Length of the array.
  \param maxIndex Reference to store the index of the maximum element.
  \param maxVal Reference to store the maximum value.
*/
template <class T>
void maxOf1DArray(const T* const, const int, int&, T&);

/*! \fn void maxOfAbs2DArray(const T* const * const, const int, const
  int, int&, int&, T&);

  \brief Finds the index and the maximum absolute value in a 2D array.
  
  \param V The 2D array of interest.
  \param nLen1 Length of 1st dimension.
  \param nLen2 Length of 2nd dimension.
  \param maxIndex1 Reference to store the 1st dimension index of the
  maximum element.
  \param maxIndex2 Reference to store the 2nd dimension index of the
  maximum element.
  \param maxVal Reference to store the maximum value.
*/
template <class T>
void maxOfAbs2DArray(const T* const * const, const int, const int,
		     int&, int&, T&);

/*! \fn void swapPtr(T*&, T*&);
  
  \brief Swaps two pointers (used in sorting multiD arrays). Note that
  the original array must be allocated using a for-loop (ie, not a
  continuous chunk), or else there would be problems when freeing the
  memory after sorting.
  
  \param a Reference to the first pointer.
  \param b Reference to the second pointer.
*/
template <class T>
void swapPtr(T*&, T*&);

/*! \fn string num2String(const T);

  \brief Converts a numerical number to a string (like the inverse of
  atoi or atof).
  
  \param input Numerical number.

  \return String form of the input number.
*/
template <class T>
string num2String(const T);

/*! \fn void minOfVector(const T&, int&, U&);

  \brief Finds the minimum value and its correponding index in a
  vector.
  
  \param V Reference to the vector of interest.
  \param minIndex Reference to the minimum index variable.
  \param minVal Reference to the minimum value variable.
*/
template <class T, class U>
void minOfVector(const T&, int&, U&);



/* Template definitions. */

template <class T, class U, class V>
void quickSort(T* const arrayA, U* const arrayB, V** const arrayC, 
	       const int left, const int right){

  /* Partitioning and sorting. */
  if (left<right){
    int pivotIndex = left;
    int pivotNewIndex = partition<T,U,V>(arrayA,arrayB,arrayC,left,right,pivotIndex);
    quickSort<T,U,V>(arrayA,arrayB,arrayC,left,pivotNewIndex-1);
    quickSort<T,U,V>(arrayA,arrayB,arrayC,pivotNewIndex+1,right);
  }

}

template <class T, class U, class V>
int partition(T* const arrayA, U* const arrayB, V** const arrayC, 
	      const int left, const int right, const int pivotIndex){

  /* Initializations. */
  T pivotValue(arrayA[pivotIndex]);

  swap<T>(arrayA[right],arrayA[pivotIndex]);
  if (arrayB!=NULL)
    swap<U>(arrayB[right],arrayB[pivotIndex]);
  if (arrayC!=NULL)
    swapPtr<V>(arrayC[right],arrayC[pivotIndex]);

  int storeIndex = left;
  for (int i=left; i<right; ++i) {
    if (arrayA[i]<pivotValue) {
      swap<T>(arrayA[i],arrayA[storeIndex]);
      if (arrayB!=NULL)
	swap<U>(arrayB[i],arrayB[storeIndex]);
      if (arrayC!=NULL)
	swapPtr<V>(arrayC[i],arrayC[storeIndex]);
      storeIndex++;

    }
  }
  swap<T>(arrayA[storeIndex],arrayA[right]);
  if (arrayB!=NULL)
    swap<U>(arrayB[storeIndex],arrayB[right]);
  if (arrayC!=NULL)
    swapPtr<V>(arrayC[storeIndex],arrayC[right]);
    
  return storeIndex;
}

template <class T>
void maxOf1DArray(const T* const V, const int length, int &maxIndex, 
		  T &maxVal){
  
  /* Initializations. */
  maxIndex = 0;
  maxVal = V[0];
  
  /* Find the maximum value and its index. */
  for (int i=1; i<length; i++){
    if (V[i]>maxVal){
      maxVal = V[i];
      maxIndex = i;
    }
  }
}

template <class T>
void maxOfAbs2DArray(const T* const * const V, const int nLen1, const int nLen2,
		     int &maxIndex1, int &maxIndex2, T &maxVal){
  
  /* Initializations. */
  maxIndex1 = 0;
  maxIndex2 = 0;
  maxVal = abs(V[0][0]);
  
  /* Find the maximum value and its index. */
  for (int i=0; i<nLen1; i++){
    for (int j=0; j<nLen2; j++){
      if (abs(V[i][j])>maxVal){
	maxVal = abs(V[i][j]);
	maxIndex1 = i;
	maxIndex2 = j;
      }
    }
  }
}

template <class T>
void swapPtr(T* &a, T* &b){

  T* c(a);
  a = b; 
  b = c;
}

template <class T>
string num2String(const T input){

  ostringstream oss;
  oss << input;
  string output = oss.str();
  oss.str("");
  
  return output;

}

template <class T, class U>
void minOfVector(const T &V, int &minIndex, U &minVal){
  
  /* Initializations. */
  minIndex = 0;
  minVal = V[0];
  
  /* Find the minimum value and its index. */
  for (unsigned int i=1; i<V.size(); i++){
    if (V[i]<minVal){
      minVal = V[i];
      minIndex = i;
    }
  }
}

#endif
