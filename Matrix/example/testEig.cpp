/*
 *  testMatrix.cpp - test the correctness of routines in 
 *        Matrix.cpp, matrixProcessing.cpp, and eigen.cpp
 *
 *    Author: Hairong Qi
 *
 *    Date: 01/11/08
 *
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "Matrix.h"             // include the Matrix class definition
#include "Pr.h"

using namespace std;

#define Usage "Usage: ./testEig filename feature-number\n"

int main(int argc, char **argv)
{

  // check to see if the number of argument is correct
  if (argc > 3) {
    cout << Usage;
    exit(1);
  }

  int n = atoi(argv[2]);

  Matrix A = readData(argv[1], n, n);
  cout << A.getCol() << ' ' << A.getRow() << endl;
  cout << A << endl;
  Matrix eigenvalue(1,n), eigenvector(n,n);

  jacobi(A, eigenvalue, eigenvector);
  eigsrt(eigenvalue, eigenvector);

  cout << eigenvalue << endl;
  cout << eigenvector << endl;

  return 0;
}
      

  
