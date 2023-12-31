// gcc -o demo_matrix demo_matrix.cpp -L/usr/local/lib -lopenblas
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cblas.h>
#include <string.h> // memset

void affiche_matrice(float*mat,int x,int y)
{int l,m;
  for (m=0;m<x;m++) 
   {for (l=0;l<y;l++)
       // printf("%.2f ",mat[l*x+m]);    // Column Major
       printf("%.2f ",mat[l+y*m]);       // Row Major
    printf("\n");
   }
  printf("\n");
}

int main(int argc, char *argv[])
{ int nobs=3;
  int nlag=2;
  int l,m;
  float *host_val,*host_res;
  float alpha,beta; 
  host_val=(float*)malloc(sizeof(float) * nobs * nlag);
  host_res=(float*)malloc(sizeof(float) * nobs * nobs); // max of nlag & nobs
  alpha=1.;beta=0.;
  for (m=0;m<nlag;m++)
    for (l=0;l<nobs;l++)
      host_val[m*nobs+l]=(float)(2*m-l);  // change assignement order
  for (m=0;m<nobs*nlag;m++) printf("%.2f ",host_val[m]); // 0.00 1.00 -1.00 0.00 -2.00 -1.00
  printf("\n");
// CblasRowMajor ou CblasColMajor ? row-major (C) or column-major (Fortran) data ordering
  memset(host_res,0,sizeof(float)*nobs * nlag);
  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nlag , nlag, nobs, alpha, host_val, nobs, host_val, nobs, beta, host_res, nlag ); // change trans/notrans
  affiche_matrice(host_res,nlag,nlag);
  memset(host_res,0,sizeof(float)*nobs * nlag);
  cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nobs , nobs, nlag, alpha, host_val, nobs, host_val, nobs, beta, host_res, nobs ); // change trans/notrans
  affiche_matrice(host_res,nobs,nobs);
  // C := alpha*op( A )*op( B ) + beta*C,
  // alpha and beta are scalars, and A, B and C are matrices, with op( A )
  // an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  // (T/N,T/N,m,n,k,alpha, A,m si N ou k si T, B, k si N ou n si T, beta, C, m)
}

// compare with a=[0 -1 -2; 2 1 0] ; a*a' ; a'*a
