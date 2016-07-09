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

#define Usage "Usage: ./testMatrix\n"

int main(int argc, char **argv)
{

  // check to see if the number of argument is correct
  if (argc > 1) {
    cout << Usage;
    exit(1);
  }

  // to create a matrix, you can either use
  //       Matrix A(3,4);
  // or 
  //       Matrix B;
  //       B.createMatrix(2,3);
  // where the declaration and creation (or memory allocation) procedures
  // are separated
  Matrix A(5,6);      // declare and create a 3x4 matrix
  Matrix B;           // declare a matrix of zero row and zero column
  B.createMatrix(2,3);

  // assign A and B
  for (int i=0; i<A.getRow(); i++)
    for (int j=0; j<A.getCol(); j++) 
      A(i,j) = i+j;

  cout << "A = " << endl;
  cout << A << endl;

  for (int i=0; i<B.getRow(); i++)
    for (int j=0; j<B.getCol(); j++) 
      B(i,j) = i*j;

  cout << "B = " << endl;
  cout << B << endl;

  // test submatrix
  Matrix iZ = subMatrix(A,0,0,1,2)->*transpose(B);
  cout << "iZ = " << endl << iZ;

  // assign A to C
  Matrix C = A;

  cout << "C = A " << endl;
  cout << C << endl;

  // test operators +=, -=, /=, *=, ->*=
  A += C;
  cout << "A += C" << endl << A << endl;

  A -= C;
  cout << "A -= C" << endl << A << endl;

  A *= C;
  cout << "A *= C" << endl << A << endl;

  A /= C;
  cout << "A /= C" << endl << A << endl;

  // crop matrix A
  Matrix D = subMatrix(A, 0, 0, 2, 3);
  
  cout << "D = subMatrix(A, 0, 0, 2, 3)" << endl;
  cout << D << endl;
  
  // test transpose()
  Matrix tD = transpose(D);

  cout << "Transpose of D (tD):" << endl;
  cout << tD << endl;

  // test inverse()
  Matrix tmp = D ->* tD;    // matrix multiplication, tmp is a square matrix
  Matrix iD = inverse(tmp);
  cout << "inverse of D*tD:" << endl;
  cout << iD << endl;

  // test getType()
  Matrix E = getType(A, 5);
  cout << "samples of type 5: " << endl;
  cout << E << endl;

  // test operators
  tmp = A + C;
  cout << "A+C" << endl;
  cout << tmp << endl;

  tmp = A - C;
  cout << "A-C" << endl;
  cout << tmp << endl;

  tmp = A * C;
  cout << "A*C" << endl;
  cout << tmp << endl;

  tmp = A / C;
  cout << "A/C" << endl;
  cout << tmp << endl;

  tmp = A + 2;
  cout << "A+2" << endl;
  cout << tmp << endl;

  tmp = A - 2;
  cout << "A-2" << endl;
  cout << tmp << endl;

  tmp = A * 2;
  cout << "A*2" << endl;
  cout << tmp << endl;

  tmp = A / 2;
  cout << "A/2" << endl;
  cout << tmp << endl;

  // test mean()
  tmp = mean(A, 6);
  cout << "mean(A)" << endl;
  cout << tmp << endl;

  // test cov()
  tmp = cov(A, 6);
  cout << "cov(A)" << endl;
  cout << tmp << endl;

  // test det()
  Matrix F(4,4);
  F(0,0) = 2; F(0,1) = 3; F(0,2) = 4; F(0,3) = 5;
  F(1,0) = 1; F(1,1) = 2; F(1,2) = 1; F(1,3) = 4;
  F(2,0) = -1; F(2,1) = 2; F(2,2) = 3; F(2,3) = 3;
  F(3,0) = 2; F(3,1) = 3; F(3,2) = 1; F(3,3) = 0;
  cout << "F" << endl;
  cout << F << endl;
  cout << "det(F)=" << endl;
  cout << det(F) << endl;

  // test detLU()
  // better performance. For large matrices, detLU() is faster
  cout << "detLU(F)=" << endl
       << detLU(F) << endl << endl;

  // test insertsort()
  Matrix G = subMatrix(F, 1, 0, 1, 3);
  Matrix sG(1,4), pG(1,4);
  insertsort(G, sG, pG);
  cout << "G=" << endl;
  cout << G << endl;
  cout << "sorted and index of G" << endl;
  cout << sG << endl;
  cout << pG << endl;

  // test eig-system analysis
  Matrix H(2,2);
  H(0,0) = 1;
  H(0,1) = 3;
  H(1,0) = 3;
  H(1,1) = 1;
  cout << "H" << endl;
  cout << H << endl;

  Matrix V(2,2), d(1,2);
  jacobi(H, d, V);
  eigsrt(d, V);
  cout << "eigen" << endl;
  cout << d << endl;
  cout << V << endl;

  Matrix Z(5,5);
  for (int i=0; i<5; i++)
    for (int j=0; j<5; j++)
      Z(i,j) = i+j;

  cout << "det(Z) = " << det(Z) << endl;

  return 0;
}
      

  
