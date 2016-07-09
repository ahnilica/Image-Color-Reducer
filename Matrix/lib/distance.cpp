/**********************************************************
 * distance.cpp - distance calculation 
 *
 *   - euc: Euclidean distance between two vectors
 *   - mah: Mahalanobis distance between a sample and a cluster
 * 
 * Author: Hairong Qi (C) hqi@utk.edu
 *
 * Created: 02/17/08
 *
 * Modification:
 *    - 10/25/13: add "const" to matrix arguments
 *
 **********************************************************/
#include "Matrix.h"
#include "Pr.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

/**
 * Calculate the Euclidean distance between two vectors.
 * @param x Vector x (row vector or col vector).
 * @param y Vector y (row vector or col vector).
 * @return Euclidean distance between x and y.
 */
double euc(const Matrix &x, const Matrix &y)
{
  int nrx, nry, ncx, ncy;
  int i, j;
  double sum = 0.0;

  nrx = x.getRow();
  nry = y.getRow();
  ncx = x.getCol();
  ncy = y.getCol();

  if (!(nrx*nry==1 || ncx*ncy==1)) {
    cout << "The euclidean() routine can only calculate distance of vectors\n";
    exit(3);
  }

  if (nrx!=nry || ncx!=ncy) {
    cout << "Two vectors are not of the same dimension\n";
    exit(3);
  }

  for (i=0; i<nrx; i++)
    for (j=0; j<ncx; j++)
      sum += (x(i,j)-y(i,j)) * (x(i,j)-y(i,j));

  return sqrt(sum);
}


/**
 * Calculate the Mahalanobis distance between a sample and a cluster.
 * @param x The sample (a column vector).
 * @param C The covariance matrix.
 * @param mu The mean (a column vector).
 * @return Mahalanobis distance between the sample and the cluster 
 *         characterized the mean and covariance.
 */
double mah(const Matrix &x, const Matrix &C, const Matrix &mu)
{
  int nrx, nrmu, nrC, ncx, ncmu, ncC;
  //  Matrix invC, diff, diffT, tmp, mdist;

  nrx = x.getRow();
  nrmu = mu.getRow();
  nrC = C.getRow();
  ncx = x.getCol();
  ncmu = mu.getCol();
  ncC = C.getCol();

  if (ncx!=1 || ncmu!=1) {
    cout << "Mahalanobis: "
	 << "the input sample and mean need to be column vectors\n";    
    exit(3);
  }

  if (nrC!=nrmu || nrC!=nrx) {
    cout << "Mahalanobis: "
	 << "the dimension of the input parameters are not consistent\n";
    exit(3);
  }
  /*
  invC = inverse(C);
  diff = x - mu;
  diffT = transpose(diff);

  tmp = diffT->*invC;
  mdist = tmp->*diff;
  */
  return sqrt(mtod(transpose(x-mu)->*inverse(C)->*(x-mu)));
}
