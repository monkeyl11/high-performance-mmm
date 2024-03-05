
#include "project.h"
#include <math.h>


int roundUp(int lcm, int n) {
  if (n % lcm != 0) {
    return n + lcm - (n % lcm);
  }
  return n;
}

// int padMatrix(double* M, double** newMat, int r, int c, int rAlign, int cAlign, int colStride) {

//   int dimR = roundUp(rAlign, r);
//   int dimC = roundUp(cAlign, c);
//   if (dimR == r && dimC == c)
//     return 0;
//   *newMat = calloc(dimR * dimC, sizeof(double));
//   for (int i = 0; i < c; i++) {
//     for (int j = 0; j < r; j++) {
//       (*newMat)[j + i * dimR] = M[j + i * colStride];
//     }
//   }
//   return 1;
// }

// void copyBack(double* orig, double* padded, int r, int c, int rAlign, int cAlign, int colStride) {
//   int dimR = roundUp(rAlign, r);
//   int dimC = roundUp(cAlign, c);
//   for (int i = 0; i < c; i++) {
//     for (int j = 0; j < r; j++) {
//       orig[j + i * colStride] = padded[j + i * dimR];
//     }
//   }
// }

void transposeMat(double** mat, int r, int c) {
  //to do for row-ordered matries
  double temp = 0;
  double* newMat = malloc(r * c * sizeof(double));
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < c; j++) {
      newMat[i * c + j] = (*mat)[j * r + i];
    }
  }
  free(*mat);
  *mat = newMat;
}



void MyGemm( int m, int n, int k, double *A, int rsA, int csA,
	     double *B, int rsB, int csB, double *C, int rsC, int csC )
{
  //Assume row-stored for all matrices if C is row-stored
  if (rsC != 1 && csC == 1) {
    // transposeMat(A, m, k);
    // transposeMat(B, k, n);
    // transposeMat(C, m, n);
    double* CTemp = calloc(m * n, sizeof(double));
    fiveloops( n, m, k, B, csB, rsB, A, csA, rsA, CTemp, csC, n);
    for (int i = 0; i < m * n; i++) {
      C[(i / m) * rsC + i % m] += CTemp[i];
    }
    free(CTemp);
    // printf("\n\n\n\n");
    // printMat(C, 4, 4);
    //transposeMat(&C, m, n, rsC);
  }
  else
    fiveloops( m, n, k, A, rsA, csA, B, rsB, csB, C, rsC, csC);
}