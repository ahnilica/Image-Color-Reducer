/**********************************************************
 * fileIO.cpp - read/write data from/to a file
 *
 *   - readData: read data from a file to a matrix
 *   - writeData: write data to a file  
 *   - readImage: read image from a file to a matrix 
 *   - writeImage: write image to a file
 * 
 * Author: Hairong Qi (C) hqi@utk.edu
 * 
 * Created: 01/11/08
 *
 * Modified:
 *   - 09/24/13: add const to filename parameters to remove warning msg
 *               in readData() (Steven Clukey)
 **********************************************************/
#include "Matrix.h"
#include "Pr.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;


/**
 * Read data from a file.
 * @param fname The file name.
 * @param nc The number of columns of the matrix.
 * @return The data in a matrix format.
 */
Matrix readData(const char *fname, int nc)
{
  Matrix data;                // matrix to store data
  char line[1024];            // temp array to store matrix line 
  char s[100];                // temp array to store one number
  FILE *fd;                   // file pointer
  char *p;
  int i, j, k;

  //open file for reading
  if ((fd = fopen(fname, "r"))==NULL) {
    perror("readData");
    exit(1);
  }

  //get the number of rows in the file
  i = 0;
  while (fgets(line,1024,fd))
    i++;
  data.createMatrix(i, nc);

  //reset file descriptor
  rewind(fd);
  i = 0;
  while (fgets(line,1024,fd)) {
    j = 0;
    p = line;
    while (j < nc) {
      k = 0;
      // bypass white space (32) or comma (44) or tab (9)
      while ((*p==32) || (*p==44) || (*p==9))    
	p++;
      while (*p==46 || (*p<58 && *p>47) || *p==43 || *p==45) { 
	            // only take numerical values and +,-,.
	s[k++] = *p++;
      }
      s[k] = '\0';
      data(i,j) = atof(s);
      j++;
    }
    i++;
  }

  return data;
  
  /*
  Matrix temp;
  ifstream infile;
  int i, j;
  int nd,           // number of data elements in total
      nr;
  double tmp;

  infile.open(fname, ios::in);
  if (!infile) {
    cerr << "READDATA: Can't read data file: " << fname << endl;
    exit(1);
  }

  // determine the number of rows
  nd = 0;
  while (infile >> tmp) 
    nd++;
  infile.clear();
  infile.seekg(0);
  nr = nd / nc;

  // allocate memory
  temp.createMatrix(nr, nc);

  // read in each data element
  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++)
      infile >> temp(i,j);

  return temp;
  */
}


/**
 * Read data from a file - overloading readData
 * @param fname The file name.
 * @param nc The number of columns of the matrix.
 * @param nr The number of rows of the matrix.
 * @return The data in a matrix format.
 */
Matrix readData(const char *fname, int nc, int nr)
{
  Matrix data;                // matrix to store data
  char line[1024];            // temp array to store matrix line 
  char s[100];                // temp array to store one number
  FILE *fd;                   // file pointer
  char *p;
  int i, j, k;

  //open file for reading
  if ((fd = fopen(fname, "r"))==NULL) {
    perror("readData");
    exit(1);
  }

  data.createMatrix(nr, nc);

  i = 0;
  while (fgets(line,1024,fd)) {
    j = 0;
    p = line;
    while (j < nc) {
      k = 0;
      // bypass white space (32) or comma (44) or tab (9)
      while ((*p==32) || (*p==44) || (*p==9))    
	p++;
      while (*p==46 || (*p<58 && *p>47) || *p==43 || *p==45) { 
	            // only take numerical values and +,-,.
	s[k++] = *p++;
      }
      s[k] = '\0';
      data(i,j) = atof(s);
      j++;
    }
    i++;
  }

  return data;
  
  /*
  Matrix temp;
  ifstream infile;
  int i, j;
  double tmp;

  infile.open(fname, ios::in);
  if (!infile) {
    cerr << "READDATA: Can't read data file: " << fname << endl;
    exit(1);
  }

  // allocate memory
  temp.createMatrix(nr, nc);

  // read in each element
  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++)
      infile >> temp(i,j);

  return temp;
  */
}


/**
 * Read data from a file to a 1-row matrix, again overloading readData
 * @param fname The file name.
 * @return The data in a 1-row matrix format.
 */
Matrix readData(const char *fname)
{
  Matrix temp;
  ifstream infile;
  int i, j;
  int nd;             // number of data elements in the file
  double tmp;

  infile.open(fname, ios::in);
  if (!infile) {
    cerr << "READDATA: Can't read data file: " << fname << endl;
    exit(1);
  }

  // find out the number of data elements in the file
  nd = 0;
  while (infile >> tmp)
    nd++;
  infile.clear();
  infile.seekg(0);

  temp.createMatrix(1, nd);

  for (i=0; i<1; i++)
    for (j=0; j<nd; j++)
      infile >> temp(i,j);

  return temp;
}


/**
 * Write data to a file
 * @param fname The file name.
 */
void writeData(Matrix &inmtx, const char *fname)
{
  ofstream outfile;
  int i, j;
  int nr, nc; 

  outfile.open(fname, ios::out | ios::binary);
  if (!outfile) {
    cerr << "WRITEDATA: Can't write data file: " << fname << endl;
    exit(1);
  }

  nr = inmtx.getRow();
  nc = inmtx.getCol();

  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++)
      outfile << inmtx(i,j) << ' ';
    outfile << endl;
  }
}


/**
 * Read image from a file. Can only read PPM image format for now.
 * @param fname The image file name.
 * @param nrImg Record the number of rows in the image.
 * @param ncImg Record the number of columns in the image.
 * @return The image read into a matrix with 3 columns (R, G, B color)
 *         and the number of row equals to the number of pixels in the image.
 */
Matrix readImage(const char *fname, int *nrImg, int *ncImg)
{
  Matrix temp;
  FILE *fp;
  int i, j;
  int ncData, nrData;
  char dummy[80];

  fp = fopen(fname, "r");
  if (!fp) {
    cerr << "READIMAGE: Can't read image file: " << fname << endl;
  }

  // identify image format
  fgets(dummy, 80, fp);

  if (!(dummy[0] == 'P' && dummy[1] == '6')) {
    cerr << "READIMAGE: Can't identify image format." << endl;
  }

  // skip the comments
  fgets(dummy, 80, fp);

  while (dummy[0] == '#') {
    fgets(dummy, 80, fp);
  }

  // read the row number and column number
  sscanf(dummy, "%d %d", &*ncImg, &*nrImg);
  nrData = *nrImg * *ncImg;      // convert an image to an nrxnc by 3 matrix
  ncData = 3;

  // read the maximum pixel value
  fgets(dummy, 80, fp);

  // read in the data elements
  temp.createMatrix(nrData, ncData);

  for (i=0; i<nrData; i++)
    for (j=0; j<ncData; j++)
      temp(i,j) = (float)fgetc(fp);

  fclose(fp);

  return temp;
}


/**
 * Write image to a file. Can only write PPM image format for now.
 * @param fname The image file name.
 * @param nrImg The number of rows in the image.
 * @param ncImg The number of columns in the image.
 */
void writeImage(const char *fname, Matrix &img, int nrImg, int ncImg)
{
  FILE *fp;
  char dummy[80];
  int i, j;

  fp = fopen(fname, "w");

  if (!fp) {
    cerr << "WRITEIMAGE: Can't write image: " << fname << endl;
  }

  fprintf(fp, "P6\n%d %d\n255\n", ncImg, nrImg);
  
  // write the data
  for (i=0; i<img.getRow(); i++)
    for (j=0; j<img.getCol(); j++)
      fputc((char)img(i,j), fp);
      
  fclose(fp);
}
