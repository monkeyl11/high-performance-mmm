
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

void transposeMat(double* mat, int r, int c) {
  //to do for row-ordered matries
  int temp = 0;
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < c; j++) {
      temp = mat[i * c + j];
      mat[i * c + j] = mat[j * r + i];
      mat[j * r + i] = temp;
    }
  }
}



void MyGemm( int m, int n, int k, double *A, int rsA, int csA,
	     double *B, int rsB, int csB,  double *C, int rsC, int csC )
{
  //Assume row-stored for all matrices if C is row-stored
  if (rsC != 1 && csC == 1) {
    transposeMat(A, m, k);
    transposeMat(B, k, n);
    transposeMat(C, m, n);
    fiveloops( n, m, k, B, csB, rsB, A, csA, rsA, C, csC, rsC);
    transposeMat(C, n, m);
  }
  else
    fiveloops( m, n, k, A, rsA, csA, B, rsB, csB, C, rsC, csC);
}

// void MyGemm2( int m, int n, int k, double *A, int rsA, int csA,
// 	     double *B, int rsB, int csB,  double *C, int rsC, int csC )
// {
//     printf("\n\n\n\n\n");
//     printMat(C, 48, 48);
//     printf("\n\n\n\n\n");
//   fiveloops( m, n, k, A, rsA, csA, B, rsB, csB, C, rsC, csC);
//     printf("\n\n\n\n\n");
//     printMat(C, 48, 48);
//     printf("\n\n\n\n\n");
// }




  
