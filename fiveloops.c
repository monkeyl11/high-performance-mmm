//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

#include "project.h"

#define DOUBLES_STORED_IN_REG 4

void ukernel( int m, int n, int k, double *A, int rsA, int csA, 
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


  for ( int p=0; p < k; p++){
    alpha_0123_p = _mm256_loadu_pd( &alpha(p * MR, 0) ) ;
    alpha_0123_p1 = _mm256_loadu_pd( &alpha(p * MR + DOUBLES_STORED_IN_REG, 0) ) ;
    //printf("ROW: %f %f %f %f\n", alpha_0123_p[0], alpha_0123_p[1], alpha_0123_p[2], alpha_0123_p[3]);
    //printf("ROW2: %f %f %f %f\n", alpha_0123_p1[0], alpha_0123_p1[1], alpha_0123_p1[2], alpha_0123_p1[3]);
    //printf("P: %d\n", p);

    beta_p_j     = _mm256_broadcast_sd( &beta( p, 0) );
    gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );
    gamma_0123_01 = _mm256_fmadd_pd( alpha_0123_p1, beta_p_j, gamma_0123_01 );

    beta_p_j     = _mm256_broadcast_sd( &beta( p, 1) );
    gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );
    gamma_0123_11 = _mm256_fmadd_pd( alpha_0123_p1, beta_p_j, gamma_0123_11 );

    beta_p_j     = _mm256_broadcast_sd( &beta( p, 2) );
    gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );
    gamma_0123_21 = _mm256_fmadd_pd( alpha_0123_p1, beta_p_j, gamma_0123_21 );

    beta_p_j     = _mm256_broadcast_sd( &beta( p, 3) );
    gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 );
    gamma_0123_31 = _mm256_fmadd_pd( alpha_0123_p1, beta_p_j, gamma_0123_31 );

    beta_p_j     = _mm256_broadcast_sd( &beta( p, 4) );
    gamma_0123_4 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_4 );
    gamma_0123_41 = _mm256_fmadd_pd( alpha_0123_p1, beta_p_j, gamma_0123_41 );

    beta_p_j     = _mm256_broadcast_sd( &beta( p, 5) );
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
void packMatrixA(double *A, int rsA, int csA, double* newMatrix) {
    for (int i = 0; i < KC * MC; i++) {
      newMatrix[i] = A[(i % MR) * rsA + csA * ((int)(i / MR) % KC) + (MR * (int)(i / (MR * KC))) * rsA];
    }
}

void packMatrixB(double *B, int rsB, int csB, double* newMatrix) {
  
}

void innerloop( int m, int n, int k, double *A, int rsA, int csA, 
	     double *B, int rsB, int csB, double *C, int rsC, int csC) 
{
  for ( int j=0; j<n; j += NR )
    for ( int i=0; i<m; i += MR )
    {
      //gamma( i,j ) += alpha( i,p ) * beta( p,j );
      //printf("NUM: %d\n", i * KC);
      ukernel(MR, NR, k, A + i * KC, rsA, csA, B + j * csB, rsB, csB, C + rsC * i + csC * j, rsC, csC );
    }
}

void fiveloops( int m, int n, int k, double *A, int rsA, int csA, 
	     double *B, int rsB, int csB,  double *C, int rsC, int csC )
{
  //REMEMBER TO FREE MATRICES
  for (int a = 0; a < n; a += NC) {
    for (int b = 0; b < k; b += KC) {
      for (int c = 0; c < m; c += MC) {
        double* packedA = malloc(sizeof(double) * KC * MC);
        packMatrixA(&alpha(c, b), rsA, csA, packedA);
        innerloop(MC, NC, KC, packedA, rsA, csA, &beta(b, a), rsB, csB, &gamma(c, a), rsC, csC );
        free(packedA);
      }
    }
  }
}


//helper method, prints a matrix
void printMat(double* x, int r, int c) {
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < c; j++) {
      printf("%d ", (int)x[j * r + i]);
    }
    printf("\n");
  }
}






  
