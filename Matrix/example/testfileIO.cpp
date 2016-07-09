/*
 *  testfileIO.cpp
 *
 *    Purpose: sample code to show how to use the fileIO and the matrix lib
 *
 *    Command-line inputs:
 *        - name of the data file (the data file shouldn't contain any
 *          textural information. Head row should be manurally deleted if any)
 *        - number of features (nf)
 *
 *    Author: Hairong Qi
 *
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "Matrix.h"             // include the Matrix class definition
#include "Pr.h"

using namespace std;

#define Usage "Usage: ./testfileIO filename nr_of_feature \n Note 1: We assume the data file is in the following format:\n \t\tfeature1 feature2 ... featuren type \n Note 2: All feature values are floating point, all type values are positive integers.\n"

int main(int argc, char **argv)
{
  int nf;                // number of features
  Matrix A, B, A1, A2;

  // check to see if the number of argument is correct
  if (argc < 3) {
    cout << Usage;
    exit(1);
  }

  // read in the number of features
  nf = atoi(argv[2]);

  // read in data from the data file and save in matrix A
  // A should have nf+1 columns with the last column indicating the 
  // category of the samples
  A = readData(argv[1], nf+1);

  // output 3 elements of the matrix for testing
  // A(0,0) refers to the item in the first row and the first column of A
  // pay attention that the index starts with 0, not 1 as used in MATLAB
  // ... and A(0,1) is the item in the first row and 2nd column, etc.
  cout << "\nOutput three elements of the matrix for testing ... \n";
  cout << "A(0,0)=" << A(0,0) << ' '
       << "A(0,1)=" << A(0,1) << ' ' 
       << "A(0,2)=" << A(0,2) << endl;

  // you can also use "<<" to output the whole matrix
  cout << "\nOutput the whole matrix for testing ... \n";
  cout << A; 

  // or output the first two rows and first two columns of the matrix
  cout << "\nOutput the upper-left 2x2 matrix for testing ...\n";
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++)
      cout << A(i,j) << ' ';
    cout << endl;             // return to another line
  }

  // or you can get the submatrix and then print it out
  cout << "\nOutput the same submatrix but use the 'subMatrix' function ...\n";
  B = subMatrix(A, 0, 0, 1, 1);
  cout << B;

  // get the samples that belong to type 1 and save it to matrix A1
  // note that you cannot use this function if the last column doesn't 
  // indicate the type of the sample
  // also note that A1 would have nf columns
  A1 = getType(A, 1);

  // get the samples who belong to class type 2 and form another matrix A2
  A2 = getType(A, 2);

  // output matrix A1 and A2
  cout << "\nType 1 samples are:  \n";
  cout << A1;

  cout << "\nType 2 samples are:  \n";
  cout << A2;

  // write matrix A1 and A2 to a file. 
  // you can see the content of files to see if the data are correctly stored
  writeData(A1, "type1.dat");
  writeData(A2, "type2.dat");

  return 0;
}
      

  
