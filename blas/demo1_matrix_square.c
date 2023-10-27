// g++ -o demo_matrix demo_matrix.cpp -L/usr/local/lib -lopenblas
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cblas.h>

void affiche_matrice(float*mat,int x,int y)
{int l,m;
  for (m=0;m<x;m++) 
   {for (l=0;l<y;l++)          // time shifted copies of the code
       printf("%.2f ",mat[m*y+l]);
    printf("\n");
   }
  printf("\n");
}

int main(int argc, char *argv[])
{ int nobs=2;
  int nlag=2;
  int l,m;
  float tab2d[nobs][nlag];  // [nbre lignes][nbre colonnes]
  float *t;
  t=&tab2d[0][0];
  for (m=0;m<nlag;m++) 
    for (l=0;l<nobs;l++)  
      tab2d[l][m]=(float)(2*m-l);
  for (m=0;m<nobs*nlag;m++) printf("%.2f ",t[m]); // 0.00 1.00 -1.00 0.00 -2.00 -1.00
  printf("\n\n");                                 // donc tab2d[l][m+1] est apres tab2d[l][m]

  float *host_val,*host_res;
  float alpha,beta; 
  host_val=(float*)malloc(sizeof(float) * nobs * nlag);
  host_res=(float*)malloc(sizeof(float) * nlag * nlag);
  alpha=1.;beta=0.;
  for (m=0;m<nobs;m++) 
    for (l=0;l<nlag;l++)          // time shifted copies of the code
      host_val[m*nlag+l]=(double)(2*m-l);
  printf("\n\n");
// CblasRowMajor ou CblasColMajor ? row-major (C) or column-major (Fortran) data ordering
  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nlag , nlag, nobs, alpha, host_val, nlag, host_val, nobs, beta, host_res, nlag ); 
  affiche_matrice(host_res,nlag,nlag);
  cblas_sgemm(CblasRowMajor, CblasTrans, CblasTrans, nlag , nlag, nobs, alpha, host_val, nlag, host_val, nobs, beta, host_res, nlag ); 
  affiche_matrice(host_res,nlag,nlag);
  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,   nlag , nlag, nobs, alpha, host_val, nlag, host_val, nobs, beta, host_res, nlag ); 
  affiche_matrice(host_res,nlag,nlag);
  cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans,   nlag , nlag, nobs, alpha, host_val, nlag, host_val, nobs, beta, host_res, nlag ); 
  affiche_matrice(host_res,nlag,nlag);

  cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nlag , nlag, nobs, alpha, host_val, nlag, host_val, nobs, beta, host_res, nlag ); 
  affiche_matrice(host_res,nlag,nlag);
  cblas_sgemm(CblasColMajor, CblasTrans, CblasTrans, nlag , nlag, nobs, alpha, host_val, nlag, host_val, nobs, beta, host_res, nlag ); 
  affiche_matrice(host_res,nlag,nlag);
  cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans,   nlag , nlag, nobs, alpha, host_val, nlag, host_val, nobs, beta, host_res, nlag ); 
  affiche_matrice(host_res,nlag,nlag);
  cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,   nlag , nlag, nobs, alpha, host_val, nlag, host_val, nobs, beta, host_res, nlag ); 
  affiche_matrice(host_res,nlag,nlag);
     // C := alpha*op( A )*op( B ) + beta*C,
     // alpha and beta are scalars, and A, B and C are matrices, with op( A )
     // an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
     // (T/N,T/N,m,n,k,alpha, A,m/k selon N ou T, B, k/n selon N ou T, beta, C, m)
}
