/*
 *  testMpp.cpp - test routine to use MPP to process the synthetic dataset
 *
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "Matrix.h"             
#include "Pr.h"

using namespace std;

#define Usage "Usage: ./testMpp training_set test_set classes features cases \n\t training_set: the file name for training set\n\t test_set: the file name for test set\n\t classes: number of classes\n\t features: number of features (dimension)\n\t cases: used for MPP, cases can only be 1, 2, or 3\n"

int main(int argc, char **argv)
{
  int nrTr, nrTe,        // number of rows in the training and test set
    nc;                  // number of columns in the data set; both the 
                         // training and test set should have the same
                         // column number                         
  Matrix Tr, Te;

  // check to see if the number of argument is correct
  if (argc < 6) {
    cout << Usage;
    exit(1);
  }

  int classes = atoi(argv[3]);   // number of classes
  int nf = atoi(argv[4]);        // number of features (dimension)
  int cases = atoi(argv[5]);     // used by MPP
 
  // read in data from the data file
  nc = nf+1; // the data dimension; plus the one label column
  Tr = readData(argv[1], nc);
  nrTr = Tr.getRow();   // get the number of rows
  Te = readData(argv[2], nc);
  nrTe = Te.getRow();   // get the number of columns

  // prepare the labels and error count
  Matrix labelMPP(nrTe, 1);  // a col vector to hold result for MPP
  int errCountMPP = 0;       // calcualte error rate for MPP

  // assign prior probability
  Matrix Pw(classes, 1);        
  for (int i=0; i<classes; i++)
    Pw(i,0) = (float)1/classes;   // if assuming equal prior probability

  // perform classification
  for (int i=0; i<nrTe; i++) {
    // classify one test sample at a time, get one sample from the test data
    Matrix sample = transpose(subMatrix(Te, i, 0, i, nf-1)); 

    // call MPP to perform classification
    labelMPP(i,0) = mpp(Tr, sample, classes, cases, Pw);
    
    // check if the classification result is correct or not
    if (labelMPP(i,0) != Te(i,nf))
      errCountMPP++;
  }

  // calculate accuracy
  cout << "The error rate using MPP is = " << (float)errCountMPP/nrTe << endl;
 
  return 0;
}
      

  
