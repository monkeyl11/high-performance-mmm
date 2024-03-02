
#include "project.h"
#include <math.h>



int lcm_() {
  int i = 1;
  while (i % MR != 0 || i % NR != 0 || i % KC != 0) {
    i++;
  }
  return i;
}

int roundUp(int lcm, int n) {
  if (n % lcm != 0) {
    return n + lcm - (n % lcm);
  }
  return n;
}

int padMatrix(double* M, double** newMat, int r, int c, int rAlign, int cAlign) {

  int dimR = roundUp(rAlign, r);
  int dimC = roundUp(cAlign, c);
  if (dimR == r && dimC == c)
    return 0;
  *newMat = calloc(dimR * dimC, sizeof(double));
  for (int i = 0; i < c; i++) {
    for (int j = 0; j < r; j++) {
      (*newMat)[j + i * dimR] = M[j + i * r];
    }
  }
  return 1;
  //printf("\nPRINTING NEWMAT\n");
  //printMat(*newMat, dimR, dimC);
}

void copyBack(double* orig, double* padded, int r, int c, int rAlign, int cAlign) {
  int dimR = roundUp(rAlign, r);
  int dimC = roundUp(cAlign, c);
  for (int i = 0; i < c; i++) {
    for (int j = 0; j < r; j++) {
      orig[j + i * r] = padded[j + i * dimR];
    }
  }
}



void MyGemm( int m, int n, int k, double *A, int rsA, int csA,
	     double *B, int rsB, int csB,  double *C, int rsC, int csC )
{

  int lcm = lcm_();
  //printf("LCM: %d\n", lcm);
  double* APad = NULL; double* BPad = NULL; double* CPad = NULL;
  char padA = 0; char padB = 0; char padC = 0;


 // printMat(A, m, k);
  //printf("\n\n\n\n\n");
  if (padMatrix(A, &APad, m, k, MC, KC)) {
    padA = 1;
  }
  else {
    APad = A;
  }
  if (padMatrix(B, &BPad, k, n, KC, NC)) {
    padB = 1;
  }
  else {
    BPad = B;
  }
  if (padMatrix(C, &CPad, m, n, MC, NC)) {
    padC = 1;
  }
  else {
    CPad = C;
  }

  // if (padA && m < 96) {
  //   printf("\n\n\n\n\n");
  //   printMat(CPad, 96, 96);
  //   printf("\n\n\n\n\n");
  // }

  if (!(padA || padB || padC))
    fiveloops( m, n, k, A, rsA, csA, B, rsB, csB, C, rsC, csC);
  else {
    fiveloops( roundUp(lcm, m), roundUp(lcm, n), roundUp(lcm, k), 
                  APad, rsA, roundUp(lcm, csA), BPad, rsB, roundUp(lcm, csB), CPad, rsC, roundUp(lcm, csC));
  }
  if (padA) {
    //copyBack(A, APad, m, k, lcm);
    free(APad);
  }
  if (padB) {
    //copyBack(B, BPad, k, n, lcm);
    free(BPad);
  }
  if (padC) {
    copyBack(C, CPad, m, n, MC, NC);
    free(CPad);
  }

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




  
