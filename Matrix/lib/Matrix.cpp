/**********************************************************************
 * Matrix.cpp - the matrix library which implements
 *              the member functions defined in Matrix.h
 *
 * Author: Hairong Qi, ECE, University of Tennessee
 *
 * Created: January 2002
 *
 * Copyright (C) hqi@utk.edu
 *********************************************************************/ 
#include "Matrix.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iomanip>
using namespace std;


/**
 * Default constructor.
 */ 
Matrix::Matrix()
{
  createMatrix(0, 0);
}


/**
 * Constructor when the number of rows and columns are known.
 * @param nr Numbers of rows.
 * @param nc Number of columns.
 * @return The created matrix.
 */
Matrix::Matrix(int nr, int nc)
{
  createMatrix(nr, nc);
}


/**
 * Copy constructor. 
 * @param m Copy matrix.
 * @return The created matrix.
 */
Matrix::Matrix(const Matrix &m)
{
  int i, j;

  createMatrix(m.getRow(), m.getCol());             // allocate memory
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      matrix[i*col+j] = m(i,j);
}


/**
 * Allocate memory for the matrix and initialize it to be a zero matrix.
 * @param nr The number of rows.
 * @param nc The number of columns.
 */
void Matrix::createMatrix(int nr, int nc)
{
  int i;

  row = nr;
  col = nc;
  matrix = (double *) new double [row * col];
  if (!matrix) {
    cout << "CREATEMATRIX: Out of memory.\n";
    exit(1);
  }

  initMatrix();
}


/**
 * Initialize the matrix.
 * @para init The value the matrix elements are initialized to. Default is 0.0.
 */
void Matrix::initMatrix(float init)
{
  int i;

  for (i=0; i<row*col; i++)
    matrix[i] = init;
}


/**
 * Destructor.  Frees up memory.
 */
Matrix::~Matrix()
{
  if (matrix)
  {
    delete [] matrix;       // free the matrix buffer
  }
}


/**
 * Returns the total number of rows in the matrix.
 * @return Total number of rows.
 * \ingroup getset
 */
int Matrix::getRow() const
{
  return row;
}


/**
 * Returns the total number of columns in the matrix.
 * @return Total number of columns.
 * \ingroup getset
 */
int Matrix::getCol() const
{
  return col;
}


/**
 * Sets the total number of rows in a matrix.
 * @param r Total number of rows.
 * \ingroup getset
 */
void Matrix::setRow(int nr)
{
  row = nr;
}


/**
 * Sets the total number of columns in a matrix.
 * @param r Total number of columns.
 * \ingroup getset
 */
void Matrix::setCol(int nc)
{
  col = nc;
}


/**
 * Overloading () operator
 * \ingroup overload
 * @param i Row
 * @param j Column
 */
double & Matrix::operator()(int i, int j) const
{
  return matrix[i*col+j];
}


/**
 * Overloading = operator
 * \ingroup overload
 * @param m The right operand
 * @return The assignment
 */
const Matrix Matrix::operator=(const Matrix &m)
{
  int i, j;

  if (this == &m)
    return *this;

  if (matrix)
    delete [] matrix;

  createMatrix(m.getRow(), m.getCol());             // allocate memory
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      matrix[i*col+j] = m(i,j);

  return *this;
}


/**
 * Overloading + operator
 * \ingroup overload
 * @param m The right operand of +
 * @return The addition
 */
Matrix Matrix::operator+(const Matrix &m) const
{
  int i, j, nr, nc;
  Matrix temp;

  nr = m.getRow();
  nc = m.getCol();

  if (nr != row || nc != col) {
    cerr << "OPERATOR+: "
	 << "Matrices are not of the same size, cannot do addition\n";
    exit(3);
  }

  temp.createMatrix(row, col);             
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      temp(i,j) = matrix[i*col+j] + m(i,j);

  return temp;
}


/**
 * Overloading += operator
 * \ingroup overload
 * @param m The right operand
 * @return The addition
 */
Matrix & Matrix::operator+=(const Matrix &m)
{
  int i, j, nr, nc;

  nr = m.getRow();
  nc = m.getCol();

  if (nr != row || nc != col) {
    cerr << "OPERATOR+: "
	 << "Matrices are not of the same size, cannot do addition\n";
    exit(3);
  }
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      matrix[i*col+j] += m(i,j);

  return *this;
}


/**
 * Overloading - operator
 * \ingroup overload
 * @param m The right operand
 * @return The subtraction
 */
Matrix Matrix::operator-(const Matrix &m) const
{
  int i, j, nr, nc;
  Matrix temp;

  nr = m.getRow();
  nc = m.getCol();

  if (nr != row || nc != col) {
    cerr << "OPERATOR-: "
	 << "Matrices are not of the same size, cannot do subtraction\n";
    exit(3);
  }

  temp.createMatrix(row, col);             
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      temp(i,j) = matrix[i*col+j] - m(i,j);

  return temp;
}


/**
 * Overloading -= operator
 * \ingroup overload
 * @param m The right operand
 * @return The subtraction
 */
Matrix & Matrix::operator-=(const Matrix &m)
{
  int i, j, nr, nc;

  nr = m.getRow();
  nc = m.getCol();

  if (nr != row || nc != col) {
    cerr << "OPERATOR-: "
	 << "Matrices are not of the same size, cannot do subtraction\n";
    exit(3);
  }
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      matrix[i*col+j] -= m(i,j);

  return *this;
}


/**
 * Overloading * operator (piecewise multiplication)
 * \ingroup overload
 * @param m The right operand
 * @return The component-wise multiplication
 */
Matrix Matrix::operator*(const Matrix &m) const
{
  int i, j, nr, nc;
  Matrix temp;

  nr = m.getRow();
  nc = m.getCol();

  if (nr != row || nc != col) {
    cerr << "OPERATOR*: "
	 << "Matrices are not of the same size, cannot do multiplication\n";
    exit(3);
  }

  temp.createMatrix(row, col);             
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      temp(i,j) = matrix[i*col+j] * m(i,j);

  return temp;
}


/**
 * Overloading *= operator (piecewise multiplication)
 * \ingroup overload
 * @param m The right operand
 * @return The component-wise multiplication
 */
Matrix & Matrix::operator*=(const Matrix &m)
{
  int i, j, nr, nc;

  nr = m.getRow();
  nc = m.getCol();

  if (nr != row || nc != col) {
    cerr << "OPERATOR*: "
	 << "Matrices are not of the same size, cannot do multiplication\n";
    exit(3);
  }
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      matrix[i*col+j] *= m(i,j);

  return *this;
}


/**
 * Overloading / operator (piecewise division)
 * \ingroup overload
 * @param m The right operand
 * @return The division
 */
Matrix Matrix::operator/(const Matrix &m) const
{
  int i, j, nr, nc;
  Matrix temp;

  nr = m.getRow();
  nc = m.getCol();

  if (nr != row || nc != col) {
    cerr << "OPERATOR/: "
	 << "Matrices are not of the same size, cannot do multiplication\n";
    exit(3);
  }

  temp.createMatrix(row, col);             
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++) {
      if (m(i,j)!=0)
	temp(i,j) = matrix[i*col+j] / m(i,j);
      else
	temp(i,j) = matrix[i*col+j] / 0.000001;
    }

  return temp;
}


/**
 * Overloading /= operator (piecewise division)
 * \ingroup overload
 * @param m The right operand
 * @return The division
 */
Matrix & Matrix::operator/=(const Matrix &m)
{
  int i, j, nr, nc;

  nr = m.getRow();
  nc = m.getCol();

  if (nr != row || nc != col) {
    cerr << "OPERATOR/: "
	 << "Matrices are not of the same size, cannot do multiplication\n";
    exit(3);
  }
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++) {
      if (m(i,j)!=0)
	matrix[i*col+j] /= m(i,j);
      else
	matrix[i*col+j] /= 0.000001;
    }

  return *this;
}


/**
 * Overloading ->* operator (matrix multiplication)
 * \ingroup overload
 * @param m The right operand
 * @return The matrix multiplication
 */
Matrix Matrix::operator->*(const Matrix &m) const 
{
  int i, j, nr, nc, p;
  Matrix temp;
  double tmp;

  nr = m.getRow();
  nc = m.getCol();

  if (col != nr) {
    cerr << "OPERATOR->*: "
	 << "Matrix size is not consistent, cannot do multiplication\n";
    exit(3);
  }

  temp.createMatrix(row, nc);             
  
  for (i=0; i<row; i++)
    for (j=0; j<nc; j++) {
      tmp = 0.0;
      for (p=0; p<col; p++) {
	tmp += matrix[i*col+p] * m(p,j);
      }
      temp(i,j) = tmp;
    }

  return temp;
}


/**
 * Overloading ->*= operator (matrix multiplication)
 * \ingroup overload
 * @param m The right operand
 * @return The matrix multiplication
 */
/*
Matrix & Matrix::operator->*=(const Matrix &m) 
{
  int i, j, nr, nc, p;
  double tmp;
  Matrix temp;

  nr = m.getRow();
  nc = m.getCol();

  if (col != nr) {
    cerr << "OPERATOR->*: "
	 << "Matrix size is not consistent, cannot do multiplication\n";
    exit(3);
  }
  
  temp.createMatrix(nr, nc);

  for (i=0; i<row; i++)
    for (j=0; j<nc; j++) {
      tmp = 0.0;
      for (p=0; p<col; p++) {
	tmp += matrix[i*col+p] * m(p,j);
      }
      temp(i,j) = tmp;
    }

  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      matrix[i*col+j] = temp(i,j);

  return *this;
}
*/

/**
 * Overloading << operator.  Output the matrix to the specified destination.
 * \ingroup overload
 * @param out The specified output stream (or output destination).
 * @param m Matrix to be output.
 * @result Output matrix to the specified file destination.
 */
ostream & operator<<(ostream &out, const Matrix &m)
{
  int i, j;
  
  for (i=0; i<m.getRow(); i++) {
    for (j=0; j<m.getCol(); j++)
      out << setw(4) << m(i,j) << ' ';
    out << endl;
  }

  return out; 
}


/**
 * Overloading / operator.  The left operand is the matrix and the right
 * is the dividend (a double point number). Each element in the matrix is 
 * divided by the double point number.
 * \ingroup overload
 * @param left Matrix as the left operand.
 * @param dev A double point number as the right operand.
 * @result Matrix divided by a double point number.
 */
Matrix Matrix::operator/(const double dev) const
{
  int i, j;
  Matrix temp;

  temp.createMatrix(row, col);

  for (i=0; i<row; i++) 
    for (j=0; j<col; j++) 
      temp(i,j) = matrix[i*col+j] / dev;

  return temp;
}


/**
 * Overloading * operator.  The left operand is the matrix and the right
 * is a double point number. Each element in the matrix is 
 * multiplied by the double point number.
 * \ingroup overload
 * @param left Matrix as the left operand.
 * @param dev A double point number as the right operand.
 * @result Matrix multiplied by a double point number.
 */

// overloading * operator
Matrix Matrix::operator*(double s) const
{
  int i, j;
  Matrix temp;

  temp.createMatrix(row, col);

  for (i=0; i<row; i++) 
    for (j=0; j<col; j++) 
      temp(i,j) = matrix[i*col+j] * s;

  return temp;
}


/**
 * Overloading + operator.  The left operand is the matrix and the right
 * is a double point number. Each element in the matrix is 
 * added by the double point number.
 * \ingroup overload
 * @param left Matrix as the left operand.
 * @param dev A double point number as the right operand.
 * @result Matrix added by a double point number.
 */
Matrix Matrix::operator+(const double s) const
{
  int i, j;
  Matrix temp;

  temp.createMatrix(row, col);

  for (i=0; i<row; i++) 
    for (j=0; j<col; j++) 
      temp(i,j) = matrix[i*col+j] + s;

  return temp;
}


/**
 * Overloading - operator.  The left operand is the matrix and the right
 * is a double point number. Each element in the matrix is 
 * subtracted by the double point number.
 * \ingroup overload
 * @param left Matrix as the left operand.
 * @param dev A double point number as the right operand.
 * @result Matrix subtracted by a double point number.
 */
Matrix Matrix::operator-(const double s) const
{
  int i, j;
  Matrix temp;

  temp.createMatrix(row, col);

  for (i=0; i<row; i++) 
    for (j=0; j<col; j++) 
      temp(i,j) = matrix[i*col+j] - s;

  return temp;
}


/**
 * Overloading - operator.  The left operand is a double point number
 * and the right is the matrix. Each element in the matrix is 
 * added with the double point number.
 * \ingroup overload
 * @param dev A double point number as the right operand.
 * @param right Matrix as the left operand.
 * @result Matrix subtracted by a double point number.
 */
Matrix operator+(const double s, const Matrix &right)
{
  return right + s;
}


/**
 * Overloading - operator.  The left operand is the matrix and the right
 * is a double point number. Each element in the matrix is 
 * subtracted by the double point number.
 * \ingroup overload
 * @param dev A double point number as the right operand.
 * @param right Matrix as the left operand.
 * @result Matrix subtracted by a double point number.
 */
Matrix operator-(const double s, const Matrix &right)
{
  int i, j;
  int nc, nr;
  Matrix temp;

  nr = right.getRow();
  nc = right.getCol();
  temp.createMatrix(nr, nc);

  for (i=0; i<nr; i++) 
    for (j=0; j<nc; j++) 
      temp(i,j) = s - right(i,j);

  return temp;
}


/**
 * Overloading / operator.  The left operand is a double point number
 * and the right is the matrix. The double point number is divided by
 * each element in the matrix.
 * \ingroup overload
 * @param dev A double point number as the right operand.
 * @param right Matrix as the left operand.
 * @result Matrix divided by a double point number.
 */
Matrix operator/(const double dev, const Matrix &right)
{
  int i, j;
  int nc, nr;
  Matrix temp;

  nr = right.getRow();
  nc = right.getCol();
  temp.createMatrix(nr, nc);

  for (i=0; i<nr; i++) 
    for (j=0; j<nc; j++) 
      temp(i,j) = dev / right(i,j);

  return temp;
}


/**
 * Overloading * operator.  The left operand is a double point number
 * and the right is the matrix. Each element in the matrix is 
 * multiplied by the double point number.
 * \ingroup overload
 * @param dev A double point number as the right operand.
 * @param right Matrix as the left operand.
 * @result Matrix multiplied by a double point number.
 */
Matrix operator*(const double s, const Matrix &right)
{
  return right * s;
}





