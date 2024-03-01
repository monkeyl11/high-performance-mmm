
#include "project.h"



int lcm_() {
  int i = 1;
  while (i % MR != 0 || i % NR != 0 || i % KC != 0) {
    i++;
  }
  return i;
}

void padMatrix(double* M, double** newMat, int r, int c, int lcm) {
  int dimR = ((r / lcm) + 1) * lcm;
  int dimC = ((c / lcm) + 1) * lcm;
  *newMat = calloc(dimR * dimC, sizeof(double));
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < c; j++) {
      (*newMat)[i * dimC + j] = M[i * c + j];
    }
  }
}

void copyBack(double* orig, double* padded, int r, int c, int lcm) {
  int dimR = ((r / lcm) + 1) * lcm;
  int dimC = ((c / lcm) + 1) * lcm;
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < c; j++) {
      orig[i * c + j] = padded[i * dimC + j];
    }
  }
}

int roundUp(int lcm, int n) {
  return ((n / lcm) + 1) * lcm;
}

void MyGemm( int m, int n, int k, double *A, int rsA, int csA,
	     double *B, int rsB, int csB,  double *C, int rsC, int csC )
{

  int lcm = lcm_();
  double* APad = NULL; double* BPad = NULL; double* CPad = NULL;
  char padA = 0; char padB = 0; char padC = 0;




  if (m % lcm != 0 || k % lcm != 0) {
    padMatrix(A, &APad, m, k, lcm);
    padA = 1;
  }
  else {
    APad = A;
  }
  if (n % lcm != 0 || k % lcm != 0) {
    padMatrix(B, &BPad, k, n, lcm);
    padB = 1;
  }
  else {
    BPad = B;
  }
  if (m % lcm != 0 || n % lcm != 0) {
    padMatrix(C, &CPad, m, n, lcm);
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
    copyBack(C, CPad, m, n, lcm);
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




  
