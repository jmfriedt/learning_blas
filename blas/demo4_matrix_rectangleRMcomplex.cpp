// gcc -o demo_matrix demo_matrix.cpp -L/usr/local/lib -lopenblas
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cblas.h>
#include <complex>

void affiche_matrice(std::complex<double>*mat,int x,int y)
{int l,m;
  for (m=0;m<x;m++) 
   {for (l=0;l<y;l++)
       // printf("%.2lf+j%.2f ",mat[l*x+m].real(),mat[l*x+m].imag());    // Column Major
       printf("%.2lf+j%.2lf ",mat[l+y*m].real(),mat[l+y*m].imag());      // Row Major
    printf("\n");
   }
  printf("\n");
}

int main(int argc, char *argv[])
{ int nobs=3;
  int nlag=2;
  int l,m;
  std::complex<double> *host_val,*host_res;
  std::complex<double> alpha,beta; 
  host_val=(std::complex<double>*)malloc(sizeof(std::complex<double>) * nobs * nlag);
  host_res=(std::complex<double>*)malloc(sizeof(std::complex<double>) * nobs * nobs); // max of nlag & nobs
  alpha=1.;beta=0.;
  for (m=0;m<nlag;m++)
    for (l=0;l<nobs;l++)
      {host_val[m*nobs+l].real((double)(3*m-l));
       host_val[m*nobs+l].imag((double)(2*m-l));
      }
  for (m=0;m<nobs*nlag;m++) printf("%.2lf+j%.2lf ",host_val[m].real(),host_val[m].imag()); // 0.00 1.00 -1.00 0.00 -2.00 -1.00
  printf("\n");
// CblasRowMajor ou CblasColMajor ? row-major (C) or column-major (Fortran) data ordering
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nlag , nlag, nobs, &alpha, host_val, nobs, host_val, nobs, &beta, host_res, nlag ); // change trans/notrans
  affiche_matrice(host_res,nlag,nlag);
  cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nobs , nobs, nlag, &alpha, host_val, nobs, host_val, nobs, &beta, host_res, nobs ); // change trans/notrans
  affiche_matrice(host_res,nobs,nobs);
  // C := alpha*op( A )*op( B ) + beta*C,
  // alpha and beta are scalars, and A, B and C are matrices, with op( A )
  // an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  // (T/N,T/N,m,n,k,alpha, A,m si N ou k si T, B, k si N ou n si T, beta, C, m)
}

// compare with a=[0 -1 -2; 2 1 0] ; a*a' ; a'*a
