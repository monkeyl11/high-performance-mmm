#include "project.h"

void ukernel( int m, int n, int k, double *A, int rsA, int csA, 
	     double *B, int rsB, int csB,  double *C, int rsC, int csC )
{
  __m256d gamma_0123_0, gamma_0123_1, gamma_0123_2, gamma_0123_3;
  __m256d alpha_0123_p, beta_p_j;


  gamma_0123_0 = _mm256_loadu_pd( &gamma(0, 0) ) ;
  gamma_0123_1 = _mm256_loadu_pd( &gamma(0, 1) ) ;
  gamma_0123_2 = _mm256_loadu_pd( &gamma(0, 2) ) ;
  gamma_0123_3 = _mm256_loadu_pd( &gamma(0, 3) ) ;

  for ( int p=0; p < k; p++){
    alpha_0123_p = _mm256_loadu_pd( &alpha(0, p) ) ;

    beta_p_j     = _mm256_broadcast_sd( &beta( p, 0) );
    gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );

    beta_p_j     = _mm256_broadcast_sd( &beta( p, 1) );
    gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );

    beta_p_j     = _mm256_broadcast_sd( &beta( p, 2) );
    gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );

    beta_p_j     = _mm256_broadcast_sd( &beta( p, 3) );
    gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 );

  }
  _mm256_storeu_pd( &gamma(0,0), gamma_0123_0 );
  _mm256_storeu_pd( &gamma(0,1), gamma_0123_1 );
  _mm256_storeu_pd( &gamma(0,2), gamma_0123_2 );
  _mm256_storeu_pd( &gamma(0,3), gamma_0123_3 );

}

// void bad_ukernel( int m, int n, int k, double *A, int rsA, int csA, 
// 	     double *B, int rsB, int csB,  double *C, int rsC, int csC ) {

//       for ( int i=0; i<MR; i++ )
//         for ( int j=0; j<NR; j++)
//           for ( int p=0; p<NR; p++) {
//             C[i + csC * j] += A[i + csA * p] * B[p + csB * j];
//           }
//        }

void fiveloops( int m, int n, int k, double *A, int rsA, int csA, 
	     double *B, int rsB, int csB,  double *C, int rsC, int csC )
{
  for ( int i=0; i<m; i += MR )
    for ( int j=0; j<n; j += NR )
    {
      //gamma( i,j ) += alpha( i,p ) * beta( p,j );
      ukernel(MR, NR, k, &alpha(i, 0), rsA, csA, &beta(0, j), rsB, csB, &gamma(i, j), rsC, csC );
    }

}



  
