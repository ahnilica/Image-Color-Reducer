/********************************************************************
 * Matrix.h - header file of the Matrix library which defines 
 *           a new class "Matrix" and the associated member functions
 *
 * Note: 
 *   This is a simple C++ library for matrix processing. 
 *   The purpose is not for high performance, but to show how 
 *   the algorithm works through programming.
 *
 * Copyright (C) hqi@utk.edu, Jan. 2002
 *
 * Other contributors:
 *   - Rui Guo: rewrote jacobi() such that results are consistent with matlab
 *   - Steven Clukey: make it more convenient for matrix manipulation
 *   - Ryan Kerekes: determinant()
 *   - Xiaoling Wang: jacobi() based on "Numerical Recipe"
 *   - Xiaoling Wang: eigsrt() based on "Numerical Recipe"
 *   - Xiaoling Wang: feature extraction related functions
 *
 * Modifications:
 *   - 10/25/13: jacobi() is rewritten by Rui Guo by decomposing a symmetric
 *               positive-definite matrix. See the following link for detail
 *               http://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm
 *   - 09/24/13: the following changes are made by Steven Clukey
 *               1. change the friend functions (matrix + - / * scalar) 
 *               into member functions and add left side operations 
 *               2. add "const" to functions to allow chaining operations
 *               3. add mtod() to convert 1x1 matrix to double
 *   - 01/12/08: separate jacobi() and eigsrt(), it does have a reason!
 *   - 01/11/08: compress the Matrix class member function,
 *               make most C-style functions
 *   - 05/03/05: combined jacobi() and eigsrt() 
 *               added distance() and insertsort()
 *   - 04/07/05: fix jacobi() that does not alter matrix (Tom Karnowski)
 *   - 04/07/05: fix readPPM and writePPM using C style file I/O
 *   - 01/13/05: delete the Data class, only keep Matrix
 *   - 01/13/05: change all the overloading operators to pass-by-reference
 *   - 01/13/05: add overloading for "<<" 
 *
 ********************************************************************/

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
using namespace std;

/** 
 * Matrix is an m x n matrix, where
 *     m - # of rows (or # of samples)
 *     n - # of columns (or # or features)
 *         (might have one more column that indicates the sample label/class)
 *     channel - # of channels
 **/   
class Matrix {
  friend ostream & operator<<(ostream &, const Matrix &);
 
 public:
  // constructors and destructor
  Matrix();                          // default constructor 
  Matrix(int,                        // constructor with row
         int);                       // column
  Matrix(const Matrix &);            // copy constructor 
  ~Matrix();                         // destructor 

  // create a matrix
  void createMatrix(int,             // row 
                    int c=1);        // column (default 1, a column vector)
  void initMatrix(float init=0.0);   // initiate the matrix component
                                     // the default is 0.0

  // get and set functions
  int getRow() const;                // get row number 
  int getCol() const;                // get column number
  void setRow(int);                  // set row number 
  void setCol(int);                  // set column number 


  // operator overloading functions
  double & operator()(int,                  // row index
                      int) const;           // column index
  const Matrix operator=(const Matrix &);   // = operator overloading
  Matrix  operator+(const Matrix &) const;  // overloading + operator
  Matrix  operator+(const double) const;    // matrix added by a scalar
  Matrix &operator+=(const Matrix &);       // overloading += operator
  Matrix  operator-(const Matrix &) const;  // overloading - operator
  Matrix  operator-(const double) const;    // matrix subtracted by a scalar
  Matrix &operator-=(const Matrix &);       // overloading -= operator
  Matrix  operator/(const Matrix &) const;  // overloading / operator
  Matrix  operator/(const double) const;    // matrix divided by a scalar
                                            // (element-wised division)
  Matrix &operator/=(const Matrix &);       // overloading /= operator
  Matrix  operator*(const Matrix &) const;  // overloading * operator 
  Matrix  operator*(const double) const;    // matrix multiplied by a scalar
                                            // (element-wised multiplication)
  Matrix &operator*=(const Matrix &);       // overloading *= operator
  Matrix  operator->*(const Matrix &) const;// overloading ->* operator 
                                            // (matrix multiplication) 
  //  Matrix & operator->*=(const Matrix &);// overloading ->*= operator

 private:
  int row;                           // number of rows / samples 
  int col;                           // number of columns / features
  double *matrix;                    // matrix buffer
};

// Left side operator overloads
Matrix operator+(const double, const Matrix &);    // matrix added by a scalar
Matrix operator-(const double, const Matrix &);    // matrix subtracted by a scalar
Matrix operator/(const double, const Matrix &);    // scalar divided by a matrix
Matrix operator*(const double, const Matrix &);    // matrix multiplied by a scalar

////////////////////////////////////
// matrix manipulation
Matrix transpose(const Matrix &);        // matrix transpose
Matrix inverse(const Matrix &);          // matrix inverse
Matrix mean(const Matrix &, const int nf);     // return a column vector with
                                   // mean value of each column
                                   // nf is the number of features
Matrix cov(const Matrix &, const int nf);      // return the covariance of a matrix

double det(const Matrix &);              // return the determinant of a matrix
double detLU(const Matrix &);            // determinant of a matrix using LU decomposition
double deteig(Matrix &);           // det of a matrix using eigenvalue decomposition
Matrix ludcmp(const Matrix &, double &); // LU decomposition

Matrix subMatrix(const Matrix &,         // crop a matrix
                 const int,                // starting row index
                 const int,                // starting column index
                 const int,                // ending row index
                 const int);               // ending column index
Matrix getType(const Matrix &,           // get samples of a certain class
               const int);               // the data type
                                   // can only be called when the last col 
                                   // is the data type

double mtod(const Matrix &);       // convert 1x1 matrix to a double precision number

// matrix sorting (only for 1-D row vectors)
void insertsort(Matrix &,          // the vector
                const Matrix &,    // the sorted vector (row vector), 
                const Matrix &);   // the index of the sorted matrix

// calculate the eigensystem
void jacobi(Matrix &,              // the input matrix
            const Matrix &,        // generated eigenvalues in vector format
            Matrix &);             // generated eigenvectors
// sort the eigenvalues
void eigsrt(Matrix &,              // the eigenvalue vector
            Matrix &);             // the eigenvector matrix
// find the maximum entry value of the off-diagnal matrix
void findOffDiagMax(Matrix &inmtx, 
		    double &maxvalue, 
		    int &xind, 
		    int &yind);

#endif











