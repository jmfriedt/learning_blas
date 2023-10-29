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
    printf(" ;\n");
   }
}

int main(int argc, char *argv[])
{ int nobs=3;
  int nlag=2;
  int l,m;

  float *host_mem,*host_res,*host_val;
  float alpha,beta; 
  host_mem=(float*)malloc(sizeof(float) * nobs * nlag);
  host_res=(float*)malloc(sizeof(float) * nlag);
  host_val=(float*)malloc(sizeof(float) * nobs);
  alpha=1.;beta=0.;
  for (m=0;m<nobs;m++)
    for (l=0;l<nlag;l++)
      host_mem[m+nobs*l]=(float)(2*m-l);
      // host_mem[m+nobs*l]=(float)(m);  // all columns identical
  for (l=0;l<nobs;l++) host_val[l]=(float)(3*m-l);
  for (m=0;m<nobs*nlag;m++) printf("%.0f ",host_mem[m]);
  printf("\n");
  printf("a=[\n");
  affiche_matrice(host_val,1,nobs);
  printf("];\n");
  printf("b=[\n");
  affiche_matrice(host_mem,nobs,nlag);
  printf("];\n");
// CblasRowMajor ou CblasColMajor ? row-major (C) or column-major (Fortran) data ordering
  printf("a*b\n");
  memset(host_res,0,sizeof(float)*nlag);
  cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, 1, nlag, nobs, alpha, host_val, nobs, host_mem, nobs, beta, host_res, 1 ); 
  affiche_matrice(host_res,1,nlag);
  memset(host_res,0,sizeof(float)*nlag);
  cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, nlag ,1, nobs, alpha, host_mem, nobs, host_val, nobs, beta, host_res, nlag ); 
  affiche_matrice(host_res,1,nlag);
  // C := alpha*op( A )*op( B ) + beta*C,
  // alpha and beta are scalars, and A, B and C are matrices, with op( A )
  // an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  // (T/N,T/N,m,n,k,alpha, A,m si N ou k si T, B, k si N ou n si T, beta, C, m)
  memset(host_res,0,sizeof(float)*nlag);
  cblas_sgemv(CblasColMajor, CblasTrans, nobs ,nlag, alpha, host_mem, nobs, host_val, 1, beta, host_res, 1 ); 
  affiche_matrice(host_res,1,nlag);
}
