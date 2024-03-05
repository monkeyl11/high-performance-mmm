//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

#include "project.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define DOUBLES_STORED_IN_REG 4

#define MCMIN ((m - c > MC) ? MC : (m - c))
#define NCMIN ((n - a > NC) ? NC : (n - a))
#define KCMIN ((b + KC > k) ? (k - b) : KC)

void ukernel(int k, double *A, int rsA, int csA, 
	     double *B, int rsB, int csB,  double *C, int rsC, int csC )
{
  __m256d gamma_0123_0, gamma_0123_1, gamma_0123_2, gamma_0123_3, gamma_0123_4, gamma_0123_5;
  __m256d gamma_0123_01, gamma_0123_11, gamma_0123_21, gamma_0123_31, gamma_0123_41, gamma_0123_51;
  __m256d alpha_0123_p, alpha_0123_p1, beta_p_j;


  gamma_0123_0 = _mm256_loadu_pd( &gamma(0, 0) ) ;
  gamma_0123_1 = _mm256_loadu_pd( &gamma(0, 1) ) ;
  gamma_0123_2 = _mm256_loadu_pd( &gamma(0, 2) ) ;
  gamma_0123_3 = _mm256_loadu_pd( &gamma(0, 3) ) ;
  gamma_0123_4 = _mm256_loadu_pd( &gamma(0, 4) ) ;
  gamma_0123_5 = _mm256_loadu_pd( &gamma(0, 5) ) ;

  gamma_0123_01 = _mm256_loadu_pd( &gamma(DOUBLES_STORED_IN_REG, 0) ) ;
  gamma_0123_11 = _mm256_loadu_pd( &gamma(DOUBLES_STORED_IN_REG, 1) ) ;
  gamma_0123_21 = _mm256_loadu_pd( &gamma(DOUBLES_STORED_IN_REG, 2) ) ;
  gamma_0123_31 = _mm256_loadu_pd( &gamma(DOUBLES_STORED_IN_REG, 3) ) ;
  gamma_0123_41 = _mm256_loadu_pd( &gamma(DOUBLES_STORED_IN_REG, 4) ) ;
  gamma_0123_51 = _mm256_loadu_pd( &gamma(DOUBLES_STORED_IN_REG, 5) ) ;


  for ( int p=0; p < k; p+=2){ //PAD KC FOR MULTIPLES OF 2!!!!!!
    alpha_0123_p = _mm256_loadu_pd( &alpha(p * MR, 0) ) ;
    alpha_0123_p1 = _mm256_loadu_pd( &alpha(p * MR + DOUBLES_STORED_IN_REG, 0) ) ;

    beta_p_j     = _mm256_broadcast_sd( &beta( NR * p, 0) );
    gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );
    gamma_0123_01 = _mm256_fmadd_pd( alpha_0123_p1, beta_p_j, gamma_0123_01 );

    beta_p_j     = _mm256_broadcast_sd( &beta( NR * p + 1, 0) );
    gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );
    gamma_0123_11 = _mm256_fmadd_pd( alpha_0123_p1, beta_p_j, gamma_0123_11 );

    beta_p_j     = _mm256_broadcast_sd( &beta( NR * p + 2, 0) );
    gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );
    gamma_0123_21 = _mm256_fmadd_pd( alpha_0123_p1, beta_p_j, gamma_0123_21 );

    beta_p_j     = _mm256_broadcast_sd( &beta( NR * p + 3, 0) );
    gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 );
    gamma_0123_31 = _mm256_fmadd_pd( alpha_0123_p1, beta_p_j, gamma_0123_31 );

    beta_p_j     = _mm256_broadcast_sd( &beta( NR * p + 4, 0) );
    gamma_0123_4 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_4 );
    gamma_0123_41 = _mm256_fmadd_pd( alpha_0123_p1, beta_p_j, gamma_0123_41 );

    beta_p_j     = _mm256_broadcast_sd( &beta( NR * p + 5, 0) );
    gamma_0123_5 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_5 );
    gamma_0123_51 = _mm256_fmadd_pd( alpha_0123_p1, beta_p_j, gamma_0123_51 );



    alpha_0123_p = _mm256_loadu_pd( &alpha((p+1) * MR, 0) ) ;
    alpha_0123_p1 = _mm256_loadu_pd( &alpha((p+1) * MR + DOUBLES_STORED_IN_REG, 0) ) ;

    beta_p_j     = _mm256_broadcast_sd( &beta( NR * (p+1), 0) );
    gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );
    gamma_0123_01 = _mm256_fmadd_pd( alpha_0123_p1, beta_p_j, gamma_0123_01 );

    beta_p_j     = _mm256_broadcast_sd( &beta( NR * (p+1) + 1, 0) );
    gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );
    gamma_0123_11 = _mm256_fmadd_pd( alpha_0123_p1, beta_p_j, gamma_0123_11 );

    beta_p_j     = _mm256_broadcast_sd( &beta( NR * (p+1) + 2, 0) );
    gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );
    gamma_0123_21 = _mm256_fmadd_pd( alpha_0123_p1, beta_p_j, gamma_0123_21 );

    beta_p_j     = _mm256_broadcast_sd( &beta( NR * (p+1) + 3, 0) );
    gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 );
    gamma_0123_31 = _mm256_fmadd_pd( alpha_0123_p1, beta_p_j, gamma_0123_31 );

    beta_p_j     = _mm256_broadcast_sd( &beta( NR * (p+1) + 4, 0) );
    gamma_0123_4 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_4 );
    gamma_0123_41 = _mm256_fmadd_pd( alpha_0123_p1, beta_p_j, gamma_0123_41 );

    beta_p_j     = _mm256_broadcast_sd( &beta( NR * (p+1) + 5, 0) );
    gamma_0123_5 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_5 );
    gamma_0123_51 = _mm256_fmadd_pd( alpha_0123_p1, beta_p_j, gamma_0123_51 );

  
  }
  _mm256_storeu_pd( &gamma(0,0), gamma_0123_0 );
  _mm256_storeu_pd( &gamma(0,1), gamma_0123_1 );
  _mm256_storeu_pd( &gamma(0,2), gamma_0123_2 );
  _mm256_storeu_pd( &gamma(0,3), gamma_0123_3 );
  _mm256_storeu_pd( &gamma(0,4), gamma_0123_4 );
  _mm256_storeu_pd( &gamma(0,5), gamma_0123_5 );
  _mm256_storeu_pd( &gamma(DOUBLES_STORED_IN_REG,0), gamma_0123_01 );
  _mm256_storeu_pd( &gamma(DOUBLES_STORED_IN_REG,1), gamma_0123_11 );
  _mm256_storeu_pd( &gamma(DOUBLES_STORED_IN_REG,2), gamma_0123_21 );
  _mm256_storeu_pd( &gamma(DOUBLES_STORED_IN_REG,3), gamma_0123_31 );
  _mm256_storeu_pd( &gamma(DOUBLES_STORED_IN_REG,4), gamma_0123_41 );
  _mm256_storeu_pd( &gamma(DOUBLES_STORED_IN_REG,5), gamma_0123_51 );


}

//packing MCxKC matrix
void packMatrixA(double *A, int rsA, int csA, double* newMatrix, int m, int k, int kcmin) {
  int i = 0;
  memset(newMatrix, 0, KC * MC * sizeof(double));
  int mRounded = roundUp(MR, m);
    while (i < KC * mRounded) {
      if ((i / MR) % KC < kcmin && i % MR +  MR * (i / (MR * KC)) < m)
        newMatrix[i] = A[(i % MR) * rsA + csA * (((i / MR)) % KC) + rsA * MR * (i / (MR * KC))];
      i ++;
    }
}

//NCxKC
void packMatrixB(double *B, int rsB, int csB, double* newMatrix, int n, int k, int kcmin) {
  memset(newMatrix, 0, KC * NC * sizeof(double));
  int i = 0;
  while (i < KC * n) {
      if ((i / NR) % KC < kcmin && i % NR + NR * (i / (NR * KC)) < n)
        newMatrix[i] = B[(i % NR) * csB + ((i / NR) * rsB) % KC + csB * NR * (i / (NR * KC))];
      i ++;
    }

}

int padMat(double* C, double* cTemp, int r, int c, int csC) {

  for (int i = 0; i < c; i++) {
    for (int j = 0; j < r; j++) {
      cTemp[j + i * MR] = C[j + i * csC];
    }
  }
  return 1;
}

int copyB(double* C, double* cTemp, int r, int c, int csC) {
    for (int i = 0; i < c; i++) {
      for (int j = 0; j < r; j++) {
        C[j + i * csC] = cTemp[j + i * MR];
      }
    }
    return 1;
}

void innerloop( int m, int n, int k, double *A, int rsA, int csA, 
	     double *B, int rsB, int csB, double *C, int rsC, int csC) 
{
  for ( int j=0; j<n; j += NR )
    for ( int i=0; i<m; i += MR )
    {
      if ((j + NR > n) || (i + MR > m)) {
          double* cTemp = calloc(MR * NR, sizeof(double));
          padMat(C + rsC * i + csC * j, cTemp, MIN(m - i, MR), MIN(n - j, NR), csC);
          ukernel(k, A + i * KC, rsA, csA, B + j * KC, rsB, csB, cTemp, rsC, MR); //cTemp???
          copyB(C + rsC * i + csC * j, cTemp, MIN(m - i, MR), MIN(n - j, NR), csC);
          free(cTemp);
      }
      else
        ukernel(k, A + i * KC, rsA, csA, B + j * KC, rsB, csB, C + rsC * i + csC * j, rsC, csC );
    }
}

void fiveloops( int m, int n, int k, double *A, int rsA, int csA, 
	     double *B, int rsB, int csB,  double *C, int rsC, int csC )
{

  double* packedB = NULL;
  double* packedA = NULL;
  packedB = malloc(KC * NC * sizeof(double));
  packedA = malloc(KC * MC * sizeof(double));
  for (int a = 0; a < n; a += NC) {
    for (int b = 0; b < k; b += KC) {
      packMatrixB(&beta(b, a), rsB, csB, packedB, NCMIN, k, KCMIN);
      for (int c = 0; c < m; c += MC) {
        packMatrixA(&alpha(c, b), rsA, csA, packedA, MCMIN, k, KCMIN);
        innerloop(MCMIN, NCMIN, KC, packedA, rsA, csA, packedB, rsB, csB, &gamma(c, a), rsC, csC);
      }
    }
  }
  free(packedA);
  free(packedB);

}

//helper method, prints a matrix
void printMat(double* x, int r, int c) {
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < c; j++) {
      printf("%f ", (x[j * r + i]));
    }
    printf("\n");
  }
}

//helper method for debugging padding
double sumMat(double* mat, int size) {
  double total = 0;
  for (int i = 0; i < size; i++) {
    total += mat[i];
  }
  return total;
}

int inBounds(int num, int max) {
  if (num >= max && num < 0) {
    printf("OUT OF BOUNDS FOR NUM %d UNDER MAX %d\n", num, max);
    return 0;
  }
  return 1;
}






  
