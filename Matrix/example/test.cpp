/*
  Testing matrix library

*/

#include <iostream>
#include <cmath>
#include "../include/Matrix.h"
#include "../include/Pr.h"
using namespace std;

int main() {

  Matrix a, b;
  
  a.createMatrix(2, 2);
  b.createMatrix(2, 2);
  
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(1, 0) = 3;
  a(1, 1) = 4;
  
  b(0, 0) = 1;
  b(0, 1) = 1;
  b(1, 0) = 1;
  b(1, 1) = 1;
  
  Matrix inva = inverse(a);
  Matrix identity = a ->* inva;
  
  for(int i = 0; i < identity.getRow(); i++) {
    for(int j=0; j < identity.getCol(); j++)   {
      if(identity(i, j) <= .00001) {
        identity(i, j) = 0;
      }
    }
  }
  
  cout << identity << endl;

  a += b;
  cout << a << endl;
  
  a -= b;
  cout << a << endl;
  
  b(1,1) = 2;
  a *= b;
  cout << a << endl;
  
  a /= b;
  cout << a << endl;
  
  a = transpose(a);
  cout << a << endl;
  
  Matrix aInv = inverse(a);
  cout << a << endl;
  
  a = a ->* aInv;
  cout << a << endl;
}