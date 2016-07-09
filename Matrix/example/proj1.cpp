/*================================

  Andrew Hnilica
  ECE - 471: Project 1
  proj1.cpp

=================================*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "Matrix.h"
#include "Pr.h"

#define _USE_MATH_DEFINES

using namespace std;

const int nf = 2;
const int nc = 2;

double likelihood2D(const Matrix &pos, const Matrix &mean, const Matrix &cov) {

  double W1, W2, W3;

  if(pos.getRow() != nf || pos.getCol() != 1) {
    cerr << "likelihood2D::invalid pos dimensions\n";
    exit(1);
  }
  if(mean.getRow() != nf || mean.getCol() != 1) {
    cerr << "likelihood2D::invalid mean dimensions\n";
    exit(1);
  }
  if(cov.getRow() != nf || cov.getCol() != nf) {
    cerr << "likelihood2D::invalid cov dimensions\n";
    exit(1);
  }
  
  //multivariate normal distribution
  W1 = 1.0 / (2 * M_PI * sqrt(det(cov)));
  W2 = mtod((transpose(pos - mean)) ->* (inverse(cov)) ->* (pos - mean));
  W3 = pow(M_E, (-.5 * W2));
   
  return(W1 * W3);  
}

int likelihoodRatio(const Matrix &pos, const Matrix &M0, const Matrix &cov0, const Matrix &M1, const Matrix &cov1) {

  double L0, L1;

  if(pos.getRow() != nf || pos.getCol() != 1) {
    cerr << "likelihood2D::invalid pos dimensions\n";
    exit(1);
  }
  
  L0 = likelihood2D(pos, M0, cov0);
  L1 = likelihood2D(pos, M1, cov1);
  
  //returns type guess based on highest likelihood ratio
  if((L0 / L1) > (L1 / L0)) {
    return 0;
  }
  else if((L0 / L1) == (L1 / L0)) {
    return -1;
  }
  else {
    return 1;
  }
}


int main(int argc, char *argv[]) {

  Matrix tr0, tr1;
  Matrix trM0, trM1;
  Matrix trCov0, trCov1;
  Matrix trData, teData;
  
  Matrix pos(2, 1);
  Matrix prior(2, 1);
  
  int i, j, k;
  int trRows, trCols;
  int teRows, teCols;
  int classVal, guess;
  
  double percentErr;
  double correct[4] = {0.0, 0.0, 0.0, 0.0};
  
  ofstream resultFile;
   
  //read data sets
  trData = readData("synth.tr", nf + 1);
  teData = readData("synth.te", nf + 1);
 
  //get size of data sets
  trRows = trData.getRow();
  trCols = trData.getCol();
  teRows = teData.getRow();
  teCols = teData.getCol();
     
  //separate by types
  tr0 = getType(trData, 0);
  tr1 = getType(trData, 1);
  
  //get mean values
  trM0 = mean(tr0, nf);
  trM1 = mean(tr1, nf);
  
  //get cov matrix
  trCov0 = cov(tr0, nf);
  trCov1 = cov(tr1, nf);
  
  //show mean values
  cout << "\n\nParameters from Training Data Set";
  cout << endl << endl;
  cout << setprecision(3);
  cout << "Mean values for Class 0:\n";
  cout << "X: " << trM0(0, 0) << endl;
  cout << "Y: " << trM0(1, 0) << endl;
  cout << endl;
  cout << "Mean values for Class 1:\n";
  cout << "X: " << trM1(0, 0) << endl;
  cout << "Y: " << trM1(1, 0) << endl;
  
  //show cov matrices
  cout << endl;
  cout << "Covariance Matrix for Class 0:\n";
  cout << left << setw(6) << trCov0(0, 0) << " " << trCov0(0, 1) << endl;
  cout << left << setw(6) << trCov0(1, 0) << " " << trCov0(1, 1) << endl;
  
  cout << "\nCovariance Matrix for Class 1:\n";
  cout << left << setw(6) << trCov1(0, 0) << " " << trCov1(0, 1) << endl;
  cout << left << setw(6) << trCov1(1, 0) << " " << trCov1(1, 1) << endl;
  
  //create table of percent error for each classifier with changing prior probabilities
  cout << "\n\nPercent Error on Testing Data Set for each Classifier\n";
  cout << "\n Prior Prob      Percent Error\n";
  cout << "Pr0:  Pr1:   LR:    Case1: Case2: Case3:  " << endl;
  
  for(i = 0; i < 11; i++) {
    //change prior probabilities for each iteration
    prior(0, 0) = i / 10.0;
    prior(1, 0) = 1 - prior(0, 0);
    cout << left << setw(4) << prior(0, 0) << "  " << setw(4) << prior(0, 1) << "   ";
    for(j = 0; j < teRows; j++) {
      pos(0, 0) = teData(j, 0);
      pos(1, 0) = teData(j, 1);
      
      classVal = teData(j, 2);
      
      //running sum of correct guesses for each classifier
      guess = likelihoodRatio(pos, trM0, trCov0, trM1, trCov1);
      if(guess == classVal) {correct[0]++;}
      
      guess = mpp(trData, pos, nc, 1, prior);
      if(guess == classVal) {correct[1]++;}
      
      guess = mpp(trData, pos, nc, 2, prior);
      if(guess == classVal) {correct[2]++;}
      
      guess = mpp(trData, pos, nc, 3, prior);
      if(guess == classVal) {correct[3]++;}
    }
    
    for(k = 0; k < 4; k++) {
      //calculate and print percent error for each classifier
      percentErr = 1.0 - (correct[k] / (double) teRows); 
      cout << setw(5) << percentErr << "  ";
      correct[k] = 0;
    }
    cout << endl;
  }
  
  //create output file for decision rules
  resultFile.open("RESULT_FILE.txt");
  
  //reset to equal prior probability
  prior(0, 0) = 0.5;
  prior(1, 0) = 0.5;
  
  //create file for each classifier
  //'X_Val, Y_Val, Guesses' is data in file
  //data can be used to plot decision rule for each classifier
  for(i = 0; i < 301; i++) {
    for(j = 0; j < 141; j++) {
      pos(0, 0) = -1.5 + (i / 100.0);
      pos(0, 1) = -0.2 + (j / 100.0);
      
      //plot over y = [-.2, 1.2] x = [-1.5, 1.5]
      //dimensions same as testing data in plotsynth.m
      
      guess = likelihoodRatio(pos, trM0, trCov0, trM1, trCov1);
      resultFile << pos(0, 0) << " "  << pos(0, 1) << " " << guess << " ";
      
      guess = mpp(trData, pos, nc, 1, prior);
      resultFile << guess << " ";
      
      guess = mpp(trData, pos, nc, 2, prior);
      resultFile << guess << " ";
      
      guess = mpp(trData, pos, nc, 3, prior);
      resultFile << guess << endl;  
    }
  }
  
  resultFile.close();
}

