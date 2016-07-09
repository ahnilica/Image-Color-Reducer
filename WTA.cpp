/*
  An implementation of the winner-takes-all clustering algorithm to approximate the optimal subset of
  colors found in the image which can be used to reduce its color set. Running with 32
  clusters will result in a new image that looks similar to the old image, but only has
  32 colors present in it.
*/

#include <ctime>
#include <iostream>
#include <cstdlib>
#include "Matrix.h"
#include "Pr.h"

using namespace std;

int main(int argc, char *argv[]) {

  //parameter values
  int clusterNumber;
  double learningPar;
  string filenamePrefix;

  bool end = false;
  int changeCount = 0;

  if(argc != 4) {
    cerr << "Invalid Command Line Args" << endl << "Usage: clusterNumber learningRate filename";
    exit(1);
  }

  cout << "Starting..." << endl << endl;

  //parse command line
  clusterNumber = atoi(argv[1]);
  learningPar = atof(argv[2]);
  filenamePrefix = argv[3];

  //read image
  int numRows, numCols;
  Matrix image = readImage("flowers.ppm", &numRows, &numCols);

  /*  image is returned as a 2d matrix where the number of rows is the total number of
      pixels in the image, and the number of columns is 3 (for r, g, b values)        */

  int numPixels = image.getRow();

  //seed random number generator
  srand(time(NULL));

  Matrix clusters(clusterNumber, 3);
  Matrix clusterAssignments(numPixels, 1);

  //initialize cluster values to random pixel value from image
  for(int i=0; i < clusterNumber; i++) {
      int randCoordinate = rand() % (120 * 120);
      for(int j=0; j < 3; j++) {
          clusters(i, j) = image(randCoordinate, j);
      }
  }

  cout << "Number of pixel-cluster assignment changes per round: " << endl;

  while(!end) {
      int winner;

      //assign each sample to a cluster
      for(int i=0; i < numPixels; i++) {

          int minDistCoordinate;

          //get sample value
          Matrix X = subMatrix(image, i, 0, i, 2);

          //for each pixel, find cluster with closest rgb value (measured in euclidean distance)
          for(int j=0; j < clusterNumber; j++) {
              double minDist;

              Matrix Y = subMatrix(clusters, j, 0, j, 2);
              int dist = euc(X, Y);


              if(j == 0) {
                  minDist = dist;
                  minDistCoordinate = j;
              }
              else if(dist < minDist) {
                  minDist = dist;
                  minDistCoordinate = j;
              }
          }

          if(clusterAssignments(i, 0) != minDistCoordinate) {
              clusterAssignments(i, 0) = minDistCoordinate;
              changeCount++;
          }
      }

      //update each cluster value for each sample
      for(int i=0; i < numPixels; i++) {
          Matrix X = subMatrix(image, i, 0, i, 2);

          winner = clusterAssignments(i, 0);
          Matrix Y = subMatrix(clusters, winner, 0, winner, 2);

          Y = Y + learningPar * (X - Y);
          clusters(winner, 0) = Y(0, 0);
          clusters(winner, 1) = Y(0, 1);
          clusters(winner, 2) = Y(0, 2);
      }

      //see if converged
      if(changeCount < 50) {
          end = true;
      }
      cout << changeCount << endl;
      changeCount = 0;
  }

  //change the value of each pixel to its cluster's value (reduces color space of image)
  for(int i=0; i < numPixels; i++) {
      int clusterIndex = clusterAssignments(i, 0);

      image(i, 0) = clusters(clusterIndex, 0);
      image(i, 1) = clusters(clusterIndex, 1);
      image(i, 2) = clusters(clusterIndex, 2);
  }

  cout << "Finished, writting image" << endl;

  string filename = filenamePrefix + ".ppm";
  writeImage(filename.c_str(), image, numRows, numCols);
}
