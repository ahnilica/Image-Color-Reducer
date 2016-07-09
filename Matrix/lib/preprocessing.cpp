/**********************************************************
 * preprocessing.cpp  
 *
 *   - normalize: normalize training and test set
 *   - pca: principal component analysis
 * 
 * Author: Hairong Qi (C) hqi@utk.edu
 *
 **********************************************************/
#include "Matrix.h"
#include "Pr.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;


/**
 * Matrix normalization.
 * @param tr The training set.
 * @param te The test set. 
 * @param nf The number of features.
 * @param flag If flag is on, it's supervised learning; otherwise, it's
 *             unsupervised learning and the second argument can be empty.
 */
void normalize(Matrix &tr, Matrix &te, const int nf, const int flag)
{
  Matrix mu, Sigma, sigma;

  // get the statistics from the training set
  mu = mean(tr, nf);

  Sigma = cov(tr, nf);
  sigma.createMatrix(nf,1);
  for (int j=0; j<nf; j++)
    sigma(j,0) = sqrt(Sigma(j,j));

  // normalize the training set
  for (int i=0; i<tr.getRow(); i++) {
    for (int j=0; j<nf; j++)
      tr(i,j) = (tr(i,j)-mu(j,0)) / sigma(j,0);
  }

  // normalize the test set
  if (flag) {
    for (int i=0; i<te.getRow(); i++) {
      for (int j=0; j<nf; j++)
	te(i,j) = (te(i,j)-mu(j,0)) / sigma(j,0);
    }
  }
}


/**
 * Principal component analysis.
 * @param tr The training set.
 * @param te The test set. 
 * @param nf The number of features.
 * @param err The error rate needs to be satisfied.
 * @param flag If flag is on, it's supervised learning; otherwise, it's
 *             unsupervised learning and the second argument can be empty.
 * @return The number of features after PCA based on "err"
 */
int pca(Matrix &tr, Matrix &te, const int nf, const float err, const int flag)
{
  Matrix Sigma, temp;
  Matrix d(1,nf),   // eigenvalue (a row vector) 
    V(nf,nf),       // eigenvector with each col an eigenvector
    pV;             // eigenvectors selected based on "err"
  int p, pnf;
  float psum, sum;

  Sigma = cov(tr, nf);
  jacobi(Sigma, d, V);
  eigsrt(d, V);   // sort the eigenvalue in the ascending order

  // determine the number of principal components to keep based on "err" given
  sum = 0.0;
  for (int i=0; i<nf; i++)
    sum += d(0,i);

  p = 0;
  psum = 0.0;
  while (psum/sum < err && p < nf) {
    psum += d(0,p);
    p++;
  }
  pnf = nf - p;

  pV = subMatrix(V,0,p,nf-1,nf-1);
  
  // perform the transformation 
  for (int i=0; i<tr.getRow(); i++) {          // for training set
    temp = subMatrix(tr,i,0,i,nf-1);
    temp = temp->*pV;
    for (int j=0; j<pnf; j++)
      tr(i,j) = temp(0,j);
  }
    
  if (flag) {
    for (int i=0; i<te.getRow(); i++) {          // for test set
      temp = subMatrix(te,i,0,i,nf-1);
      temp = temp->*pV;
      for (int j=0; j<pnf; j++)
	te(i,j) = temp(0,j);
    }
  }

  return pnf;
}
