/**********************************************************
 * eigen.cpp - eigen system analysis 
 *   - jacobi: computes eigenvalues and eigenvectors of a
 *             real symmetric matrix
 *   - eigsrt: sort the eigenvalues and eigenvectors
 * 
 * Author: Xiaoling Wang, based on "Numerical Recipe"
 * 
 * Created: Jan. 2002
 *
 * Modified:
 *   - Rui Guo: rewrote jacobi() to fix inconsistent results
 *         between matlab and the original function
 *         added 
 *   - to make it C-stype (H.Q.)
 *
 **********************************************************/
#include "Matrix.h"
#include "Pr.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

/**
 * Computes all eigenvalues and eigenvectors of a real symmetric matrix
 * @param inmtx The input matrix.
 * @param d The eigenvalue in a row vector
 * @param V The eigenvector (each column is an eigenvector)
 */
// Modification: Rui Guo rewrote the function
// Modification: Tom Karnowski found that after jacobi,
//     the matrix is changed. The code actually outputs the matrix with
//     its upper diagonal components set to zero
//     Made an easy but dumb fix. 
void jacobi(Matrix &inmtx, const Matrix &D, Matrix &V)
{
  int i, j, xind, yind, ip, iq;
  double maxvalue;
  int nr, nc;
  double epsilon = 0.0000001;
	
  nr = inmtx.getRow();
  nc = inmtx.getCol();

  if (nr != nc) {
    cout << "jacobi: the input matrix is not square.\n";
    exit(3);
  }
  
  for (i=0;i<nr;i++){
    for (j=0;j<nc;j++){
      if (fabs(inmtx(i,j)-inmtx(j,i))>epsilon){
	cout<<"jacobi: the input matrix should be symmetric.\n";
	exit(3);		
      }
    }
  }

  for (ip=0; ip<nr; ip++) {
    for (iq=0; iq<nr; iq++)  
      V(ip,iq) = 0.0;
    V(ip,ip) = 1.0;
  }
  maxvalue = 0.1;
  xind = 0;
  yind = 0;
  findOffDiagMax(inmtx,maxvalue,xind,yind);

  Matrix G(nr,nc);
  while (fabs(maxvalue)>epsilon) { 	
    maxvalue = 0.0;
    for (ip=0; ip<nr; ip++) {
      for (iq=0; iq<nr; iq++)  
	G(ip,iq) = 0.0;
      G(ip,ip) = 1.0;
    }
    
    double phi = atan2(2*inmtx(xind,yind),inmtx(xind,xind)-inmtx(yind,yind))/2;
    G(xind,xind) = cos(phi);
    G(xind,yind) = -sin(phi);
    G(yind,xind) = sin(phi);
    G(yind,yind) = cos(phi);
	
    inmtx = transpose(G)->*inmtx->*G;
    V = V->*G;
    for (i=0;i<nr;i++){
      double temp = inmtx(i,i);
      D(0,i) = temp;
    }
    
    findOffDiagMax(inmtx,maxvalue,xind,yind); 
  } 	
}


/**
 * Find the maximum entry value of the off-diagnal matrix
 */
void findOffDiagMax(Matrix &inmtx, double &maxvalue, int &xind, int &yind)
{
  int nr,nc;	

  maxvalue = 0.0;
  xind = 0;
  yind = 0;
  nr = inmtx.getRow();
  nc = inmtx.getCol();

  for (int i=0;i<nr-1;i++){
    for(int j=i+1;j<nc;j++){
      if(fabs(maxvalue)<fabs(inmtx(i,j))){
	maxvalue=inmtx(i,j);
	xind = i;
	yind = j;
      }
    }
  }
}

/*
#define ROTATE(a,i,j,k,l) g=a(i,j);h=a(k,l);a(i,j)=g-s*(h+g*tau);\
   a(k,l)=h+s*(g-h*tau);
*/
/**
 * Computes all eigenvalues and eigenvectors of a real symmetric matrix
 * @param inmtx The input matrix.
 * @param d The eigenvalue in a row vector
 * @param V The eigenvector (each column is an eigenvector)
 */
// Modification: Tom Karnowski found that after jacobi,
//     the matrix is changed. The code actually outputs the matrix with
//     its upper diagonal components set to zero
//     Made an easy but dumb fix. 
/*
void jacobi(Matrix &inmtx, const Matrix &d, const Matrix &V)
{
  int j, iq, ip, i;
  double tresh, theta, tau, t, sm, s, h, g, c;
  Matrix b, z;
  Matrix tmp;
  int nr, nc;

  nr = inmtx.getRow();
  nc = inmtx.getCol();

  if (nr != nc) {
    cout << "jacobi: the input matrix is not square.\n";
    exit(3);
  }

  // to keep a copy of original matrix
  tmp = inmtx;

  b.createMatrix(1, nr);
  z.createMatrix(1, nr);
  
  // initialize V to an identity matrix
  for (ip=0; ip<nr; ip++) {
    for (iq=0; iq<nr; iq++)  
      V(ip,iq) = 0.0;
    V(ip,ip) = 1.0;
  }

  // initialize b and d to the diagonal of matrix
  for (ip=0; ip<nr; ip++) {
    b(0,ip) = inmtx(ip,ip);;
    d(0,ip) = inmtx(ip,ip);
    z(0,ip) = 0.0;
  }

  for (i=0; i<50; i++) {
    sm = 0.0;
    // sum of off-diagonal elements
    for (ip=0; ip<nr-1; ip++) {
      for (iq=ip+1; iq<nr; iq++)
        sm += (double)fabs(inmtx(ip,iq));
    }
    if (sm==0.0) {
      inmtx = tmp;
      return;
    }
    if (i<3)
      tresh = 0.2*sm/(double)(nr*nr);
    else
      tresh = 0.0;
    for (ip=0; ip<nr-1; ip++) {
      for (iq=ip+1; iq<nr; iq++) {
        g = 100.0*(double)fabs(inmtx(ip,iq));
        if (i>3 && ((double)fabs(d(0,ip))+g)==(double)fabs(d(0,ip))
           && ((double)fabs(d(0,iq))+g)==(double)fabs(d(0,iq)))
          inmtx(ip,iq) = 0.0;
        else if ((double)fabs(inmtx(ip,iq))>tresh) {
          h = d(0,iq) - d(0,ip);
          if ((double)fabs(h)+g==(double)fabs(h))
            t = inmtx(ip,iq)/h;
          else {
            theta = 0.5*h/inmtx(ip,iq);
            t = 1.0/((double)fabs(theta)+sqrt(1.0+theta*theta));
            if (theta<0.0) t = -t;
          }
          c = 1.0/sqrt(1+t*t);
          s = t*c;
          tau = s/(1.0+c);
          h = t*inmtx(ip,iq);
          z(0,ip) -= h;
          z(0,iq) += h;
          d(0,ip) -= h;
          d(0,iq) += h;
          inmtx(ip,iq) = 0.0;
          for (j=0; j<ip-1; j++) {
            ROTATE(inmtx,j,ip,j,iq);
          }
          for (j=ip+1; j<iq-1; j++) {
            ROTATE(inmtx,ip,j,j,iq);
          }
          for (j=iq+1; j<nr; j++) {
            ROTATE(inmtx,ip,j,iq,j);
          }
          for (j=0; j<nr; j++) {
            ROTATE(V,j,ip,j,iq);
          }
        }
      }
    }
    for (ip=0; ip<nr; ip++) {
      b(0,ip) += z(0,ip);
      d(0,ip) = (double)b(0,ip);
      z(0,ip) = 0.0;
    }
  }

  // set the upper diagonal of matrix back to the original
  inmtx = tmp;
}
*/

/**
 * Sort the eigenvalues and rearrange eigenvectors
 * @param
 */
void eigsrt(Matrix &d, Matrix &V)
{
  Matrix sd, pd, sV;
  int nc, i, j;

  nc = d.getCol();
  sd.createMatrix(1, nc);
  pd.createMatrix(1, nc);
  sV.createMatrix(nc, nc);

  insertsort(d, sd, pd);

  d = sd;

  for (i=0; i<nc; i++) 
    for (j=0; j<nc; j++) 
      sV(j,i) = V(j,(int)pd(0,i));    

  V = sV;
}

