/**********************************************************
 * matrixProcessing.cpp - matrix manipulation 
 *
 *   - transpose: matrix transpose
 *   - inverse: matrix inverse
 *   - subMatrix: crop a matrix
 *   - insertsort: sort an array
 *   - cov: covariance of a matrix
 *   - det: determinant of a matrix based on recursive method
 *   - detLU: determinant of a matrix based on LU decomposition
 *   - deteig: determinant of a matrix based on eigenvalue decomposition
 *   - mean: mean of each column of a matrix
 *   - getType: get the subset of the matrix of a certain type
 * 
 * Author: Hairong Qi (C) hqi@utk.edu
 * 
 * Other contributors:
 *   - Ryan Kerekes: det()
 *   - Sangwoo Moon: detLU() and deteig()
 *   - Steven Clukey: mtod()
 *
 * Created: 01/11/08
 * 
 * Modified: 
 *   - 09/24/13: add "const" to functions to allow chaining operation
 *   - 09/24/13: add mtod() to convert 1x1 matrix to double
 *   - 01/22/09: add Sangwoo's implementation on det() calculation
 *
 **********************************************************/
#include "Matrix.h"
#include "Pr.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;


/**
 * Returns the samples of a certain type and form another matrix
 * this function can only be called if the last column of the original matrix
 * indicates the class label
 * @return Samples of a certain class.
 */
Matrix getType(const Matrix &inmtx, const int t)
{
  int i, j, p, q, nr, nc, mr, mc;

  mr = inmtx.getRow();
  mc = inmtx.getCol();

  nr = 0;                     // initialize the row number of the new matrix 
  nc = inmtx.getCol() - 1;    // the column number of the new matrix

  for (i=0; i<mr; i++)
    if (inmtx(i,mc-1) == t)
      nr++;

  Matrix temp(nr, nc);

  i = 0;
  j = 0;
  for (p=0; p<mr; p++) {
    if (inmtx(p,mc-1) == t) {
      for (q=0; q<nc; q++) {
	temp(i,j) = inmtx(p,q);
	j++;
      }
      j = 0;
      i++;
    }
  }

  return temp;
}


// functions used by matrix inverse calculation
int findPivot(Matrix &, int);
void switchRow(Matrix &, int, int);
void dividePivot(Matrix &, int);
void eliminate(Matrix &, int);

/**
 * Calculate transpose of a matrix
 * @param inmtx The input matrix.
 * @return The transpose of the input matrix.
 */
Matrix transpose(const Matrix &inmtx) {
  int i, j;
  int nr, nc;  
  Matrix temp;

  nr = inmtx.getRow();
  nc = inmtx.getCol();    
  temp.createMatrix(nc, nr);

  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++)
      temp(j,i) = inmtx(i,j);

  return temp;
}


/**
 * Calculate the inverse of a matrix using Gauss-Jordan elimination
 * @param inmtx The input matrix.
 * @return The inverse of the input matrix.
 */
Matrix inverse(const Matrix &inmtx) {
  int i, j, m, k, pivot;
  int nr, nc, nchan;
  Matrix temp, mtx;

  nr = inmtx.getRow();
  nc = inmtx.getCol();
  if (nr != nc) {
    cout << "inverse: Cannot compute the inverse of a non-square matrix\n";
    exit(3);
  }
    
  temp.createMatrix(nr, 2*nc);

  // construct the A | I matrix
  k = 0;
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++)
      temp(i,j) = inmtx(i,j);
    temp(i,nc+k++) = 1;
  }

  pivot = 0;  
  while (pivot < nr) {
    m = findPivot(temp, pivot);    // partial pivoting

    if (m != pivot) {
      switchRow(temp, pivot, m);
    }
    else if (temp(pivot,pivot) == 0) {
      cout << "inverse: This is a singular matrix, noninvertible\n";
      exit(3);
    }

    if (temp(pivot, pivot) != 1.0) {
      dividePivot(temp, pivot);
    }
      
    eliminate(temp, pivot);

    pivot++;
  }

  mtx = subMatrix(temp, 0, nc, nr-1, 2*nc-1);

  return mtx;
}

// find the maximum absolute value in the pivotal column
int findPivot(Matrix &mtx, int pivot) {
  float maxi;
  int i, index;
  int nr;

  nr = mtx.getRow();
  maxi = fabs(mtx(pivot,pivot));
  index = pivot;

  for (i=pivot+1; i<nr; i++)
    if (fabs(mtx(i,pivot)) > maxi) {     // Chris reported problem with abs
      index = i;
      maxi = fabs(mtx(i,pivot));
    }

  return index;
}

// switch two rows pivot <--> m
void switchRow(Matrix &mtx, int pivot, int m) {
  int i, j;
  float tmp;
  int nc;

  nc = mtx.getCol();
  for (j=pivot; j<nc; j++) {
    tmp = mtx(pivot,j);
    mtx(pivot,j) = mtx(m,j);
    mtx(m,j) = tmp;
  }
}

// divide the pivotal row by pivot
void dividePivot(Matrix &mtx, int pivot) {
  float scale;
  int j;
  int nc;

  nc = mtx.getCol();
  scale = mtx(pivot,pivot);

  for (j=pivot; j<nc; j++) 
    mtx(pivot,j) /= scale;
}

// eliminate elements on the pivotal column
void eliminate(Matrix &mtx, int pivot) {
  int i, j;
  int nr, nc;
  float scale;

  nr = mtx.getRow();
  nc = mtx.getCol();
  
  for (i=pivot+1; i<nr; i++) {
    scale = mtx(i,pivot);
    for (j=pivot; j<nc; j++) 
      mtx(i,j) -= mtx(pivot,j) * scale;
  }

  for (i=pivot-1; i>=0; i--) {
    scale = mtx(i,pivot);
    for (j=pivot; j<nc; j++) 
      mtx(i,j) -= mtx(pivot,j) * scale;
  }
}
 

/**
 * Crop a matrix.
 * @param inmtx The input matrix.
 * @param srow The row index of the top-left corner, starts from 0.
 * @param scol The column index of the top-left corner.
 * @param erow The row index of the lower-right corner.
 * @param ecol The column index of the lower-right corner.
 * @return The cropped matrix.
 */
Matrix subMatrix(const Matrix &inmtx, const int srow, const int scol, const int erow, const int ecol) {
  int i, j;
  int nr, nc;
  Matrix temp;

  nr = inmtx.getRow();
  nc = inmtx.getCol();
  if (srow > erow || scol > ecol || erow > nr -1 || ecol > nc - 1) {
    cout << "subMatrix: Check the cropping index.\n";
    exit(3);
  }
  temp.createMatrix(erow-srow+1, ecol-scol+1);

  for (i=srow; i<=erow; i++)
    for (j=scol; j<=ecol; j++)
      temp(i-srow,j-scol) = inmtx(i,j);

  return temp;
}


/**
 * Calculate the mean vector of a matrix.
 * @param inmtx The input matrix.
 * @param nf The number of features.
 * @return A column vector of mean of each column.
 */
Matrix mean(const Matrix &inmtx, const int nf)
{
  Matrix temp;
  int nr, nc;
  int i, j;

  // each item in the mean vector is the mean value of each column
  nr = inmtx.getRow();
  nc = inmtx.getCol();
  if (nf > nc) {
    cout << "MEAN: the input feature number is too large.\n";
    exit(3);
  }
  
  temp.createMatrix(nf, 1);

  for (i=0; i<nr; i++) {
    for (j=0; j<nf; j++)
      temp(j,0) += inmtx(i,j);
  }

  for (j=0; j<nf; j++)
    temp(j,0) /= (double)nr;

  return temp;
}


/**
 * Calculate the covariance of a matrix.
 * @param inmtx The input matrix.
 * @param nf The number of features.
 * @return The covariance matrix.
 */
Matrix cov(const Matrix &inmtx, const int nf)
{    
  Matrix mu, temp, x, y;
  int nr, nc;

  nr = inmtx.getRow();
  nc = inmtx.getCol();
  if (nf > nc) {
    cout << "COV: the input feature number is too large.\n";
    exit(3);
  }

  // get the mean vector
  mu = mean(inmtx, nf);          

  // the size of the covariance matrix should be col x col 
  temp.createMatrix(nf, nf);

  /*
  for (i=0; i<nr; i++) {
    x = subMatrix(inmtx, i, 0, i, nc-1);
    y = transpose(x) - mu;
    temp = temp + y ->* transpose(y);
  }
  */

  // calculate the upper-right corner of the covariance matrix
  for (int i=0; i<nr; i++) 
    for (int m=0; m<nf; m++) 
      for (int n=m; n<nf; n++) 
	temp(m,n) += (inmtx(i,m) - mu(m,0)) * (inmtx(i,n) - mu(n,0));
    
  // copy the lower-left corner
  for (int m=0; m<nf; m++)
    for (int n=0; n<m; n++)
      temp(m,n) = temp(n,m);

  // normalize by nr-1
  if (nr > 1) {
    for (int i=0; i<nf; i++)
      for (int j=0; j<nf; j++)
	temp(i, j) /= (double)(nr-1);
  }

  return temp;
}


/**
 * Calculate the determinant of a matrix.
 * @param inmtx The input matrix.
 * @return The determinant.
 */
double det(const Matrix &inmtx)
{
  double d, cof;
  int size, i, j1, j2, c, nr, nc;
  Matrix *temp;

  nr = inmtx.getRow();
  nc = inmtx.getCol();

  if (nr != nc) {
    cout << "DET: Tried to compute determinant of non-square matrix\n";
    exit(3);
  }
  if (nr == 0) return (0);
  if (nr == 1) return (inmtx(0,0));
  if (nr == 2) return (inmtx(0,0)*inmtx(1,1)-inmtx(0,1)*inmtx(1,0));
  else {
    d = 0;
    size = nr;
    for(c=0; c<size; c++) {
      temp = new Matrix(size-1,size-1);
      for(i=1;i<size;i++) {
	j2 = 0;
	for(j1=0;j1<size;j1++) {
	  if(j1!=c) {
	    (*temp)(i-1,j2) = inmtx(i,j1);
	    j2++;
	  }
	}
      }
      cof = det(*temp);
      cof *= inmtx(0,c);
      cof *= pow((double)-1,(double)c);
      d += cof;
      delete temp;
    }
    return (d);
  }
}


/**
 * Calculate the determinant of a matrix based on LU decomposition.
 * @param inmtx The input matrix.
 * @return The determinant.
 */
double detLU(const Matrix& A)
{
  double val;
  Matrix Ap;

  Ap = ludcmp(A,val);
  for (int i=0; i<Ap.getRow(); i++) 
    val *= Ap(i,i);

  return val;
}


/**
 * Calculate the determinant of a matrix based on eigenvalue decomposition.
 * Note that it cannot handle complex eigenvalues.
 * @param inmtx The input matrix.
 * @return The determinant.
 */
double deteig(Matrix& A)
{
  int n = A.getCol();
  Matrix d(1,n),   	// eigenvalue (a row vector) 
         V(n,n);      // eigenvector with each col an eigenvector

  jacobi(A, d, V);
  double val=1;
  for (int i=0; i<n; i++) 
    val=val*d(1,i);

  return val;
}


/**
 * LU decomposition.
 * Reference: Numerical Recipes in C, 2nd Edition, Cambridge University Press
 */
Matrix ludcmp(const Matrix &src, double& d)
{
  int n = src.getRow();
  if (n != src.getCol()) {
    cout << "ludcmp: No squre matrix input.\n";
    exit(1);
  }

  double tiny=1.0e-20;
  int i,imax,j,k;
  double big,dum,sum,temp;
  Matrix vv(1,n);

  Matrix a=src;
  d=1;

  for (i=0; i<n; i++) {
    big=0;
    for (j=0; j<n; j++)
      if ((temp=fabs(a(i,j))) > big) big=temp;
    if (big==0) {
      cout << "ludcmp: singular matrix.\n";
      exit(2);
    }
    vv(0,i)=1.0/big;
  }

  for (j=0; j<n; j++) {
    for (i=0; i<j; i++)	{
      sum=a(i,j);
      for (k=0; k<i; k++) sum-=a(i,k)*a(k,j);
      a(i,j)=sum;
    }
    big=0;
    for (i=j; i<n; i++) {
      sum=a(i,j);
      for (k=0; k<j; k++) sum-=a(i,k)*a(k,j);
      a(i,j)=sum;
      if ((dum=vv(0,i)*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j!=imax) {
      for (k=0; k<n; k++) {
	dum=a(imax,k);
	a(imax,k)=a(j,k);
	a(j,k)=dum;
      }
      d=-d;
      vv(0,imax)=vv(0,j);
    }
    if (a(j,j)==0) a(j,j)=tiny;
    if (j!=n) {
      dum=1.0/a(j,j);
      for (i=j+1; i<n; i++) a(i,j)*=dum;
    }
  }
  
  return a;
}

/**
 * Convert a 1x1 matrix to a double precision number
 * @param inmtx The input 1x1 matrix.
 * @return The double precision value of the matrix
 */
double mtod(const Matrix &inmtx)
{
  if (inmtx.getRow() != 1 || inmtx.getCol() != 1) {
    cout << "MTOD: The matrix is not 1x1.\n";
    exit(3);
  }
  
  return inmtx(0,0);
  
}

/**
 * Insert sort.
 * @param inmtx The input 1-row array.
 * @param sorted The sorted array in the ascending order.
 * @param pos The index of the sorted array.
 */
void insertsort(Matrix &inmtx, const Matrix &sorted, const Matrix &pos)
{
  double inserter;
  int i, l, index;
  int nr, nc;

  nr = inmtx.getRow();
  nc = inmtx.getCol();

  if (nr>1 && nc>1) {
    cout << "INSERTSORT: The matrix to be sorted needs to be a vector\n";
    exit(3);
  }
  if (nc==1) 
    inmtx = transpose(inmtx);

  for (i=0; i<nc; i++) {
    sorted(0,i) = inmtx(0,i);
    pos(0,i) = i;
  }
  
  for (i=1; i<nc; i++) {
    inserter = sorted(0,i);
    l = (int)pos(0,i);
    index = i-1;
    while (index>=0 && inserter<sorted(0,index)) {
      sorted(0,index+1) = sorted(0,index);
      pos(0,index+1) = pos(0,index);
      index--;
    }
    sorted(0,index+1) = inserter;
    pos(0,index+1) = l;
  }
}

