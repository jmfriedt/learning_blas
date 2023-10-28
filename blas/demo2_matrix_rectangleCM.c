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
       printf("%.2f ",mat[l*x+m]);
    printf("\n");
   }
  printf("\n");
}

int main(int argc, char *argv[])
{ int nobs=3;
  int nlag=2;
  int l,m;
  float tab2d[nobs][nlag];  // [nbre lignes][nbre colonnes]
  float *t;
  t=(float*)tab2d;
  for (m=0;m<nlag;m++) 
    for (l=0;l<nobs;l++)  
      tab2d[l][m]=(float)(2*m-l);
  for (m=0;m<nobs*nlag;m++) printf("%.2f ",t[m]); // 0.00 1.00 -1.00 0.00 -2.00 -1.00
  printf("\n");                                 // donc tab2d[l][m+1] est apres tab2d[l][m]

  float *host_val,*host_res;
  float alpha,beta; 
  host_val=(float*)malloc(sizeof(float) * nobs * nlag);
  host_res=(float*)malloc(sizeof(float) * nobs * nobs); // max of nlag & nobs
  alpha=1.;beta=0.;
  for (m=0;m<nlag;m++)
    for (l=0;l<nobs;l++)
      host_val[m+nlag*l]=(double)(2*m-l);
  for (m=0;m<nobs*nlag;m++) printf("%.2f ",host_val[m]); // 0.00 1.00 -1.00 0.00 -2.00 -1.00
  printf("\n");
// CblasRowMajor ou CblasColMajor ? row-major (C) or column-major (Fortran) data ordering
  memset(host_res,0,sizeof(float)*nobs * nlag);
  cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, nlag , nlag, nobs, alpha, host_val, nobs, host_val, nobs, beta, host_res, nlag ); 
  affiche_matrice(host_res,nlag,nlag);
  memset(host_res,0,sizeof(float)*nobs * nlag);
  cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, nobs , nobs, nlag, alpha, host_val, nobs, host_val, nobs, beta, host_res, nobs ); 
  affiche_matrice(host_res,nobs,nobs);
  // C := alpha*op( A )*op( B ) + beta*C,
  // alpha and beta are scalars, and A, B and C are matrices, with op( A )
  // an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  // (T/N,T/N,m,n,k,alpha, A,m si N ou k si T, B, k si N ou n si T, beta, C, m)
}

// compare with a=[0 1; 2 -2; -1 0] ; a'*a ; a*a'
