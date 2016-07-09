/*
 *  testimageIO.cpp - test routines in fileIO.cpp
 *
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "Matrix.h"             
#include "Pr.h"

using namespace std;

#define Usage "Usage: ./testimageIO in_filename out_filename\n"

int main(int argc, char **argv)
{
  int nr,                // number of rows in the image
    nc;                  // number of columns in the image
  Matrix data;

  // check to see if the number of argument is correct
  if (argc < 3) {
    cout << Usage;
    exit(1);
  }

  // read in data from the image file and save in a matrix
  data = readImage(argv[1], &nr, &nc);
  cout << "The size of the input image is " << nr << " x " << nc << endl;
  cout << "The size of the matrix data is " << data.getRow() << " x " 
       << data.getCol() << endl;

  // write data to an image file
  writeImage(argv[2], data, nr, nc);

  return 0;
}
      

  
