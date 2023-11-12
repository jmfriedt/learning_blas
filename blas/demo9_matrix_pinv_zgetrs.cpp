// g++ -o demo_matrix demo_matrix.cpp -L/usr/local/lib -lopenblas
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#include <cblas.h>
#include <lapacke.h>
#include <complex>

#undef debug

void affiche_matrice(std::complex<double>*mat,int x,int y)
{int l,m;
  for (m=0;m<x;m++) 
   {for (l=0;l<y;l++)
       printf("(%.5lf)+j*(%.5lf) ",real(mat[l*x+m]),imag(mat[l*x+m]));
    printf(" ;\n");
   }
  printf("\n");
}

int main(int argc, char *argv[])
{ int nobs=2100;
  int nlag=30;
  const int N=2*nlag+1;
  int l,m;
  int *IPIV;
  struct timeval tv1,tv2;
  std::complex<double> *host_mem,*host_code,*host_val,*host_res,*host_out,*host_final,*host_Id;
  std::complex<double> alpha,beta,pwr; 
  host_val=(std::complex<double>*)malloc(sizeof(std::complex<double>) * nobs);
  alpha=1.;beta=0.;
//https://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
  host_code=(std::complex<double>*)malloc(sizeof(std::complex<double>)*nobs);                  // known code
  host_res=(std::complex<double>*)malloc(sizeof(std::complex<double>)*(nlag*2+1)*(nlag*2+1));  // intermediate matrix A^h*A
  host_Id=(std::complex<double>*)malloc(sizeof(std::complex<double>)*(nlag*2+1)*(nlag*2+1));  // intermediate matrix A^h*A
  host_mem=(std::complex<double>*)malloc(sizeof(std::complex<double>)*(2*nlag+1)*nobs);        // time delayed copies of the code
  host_out=(std::complex<double>*)malloc(sizeof(std::complex<double>)*(2*nlag+1)*nobs);        // A*(A^h*A)^-1
  host_final=(std::complex<double>*)malloc(sizeof(std::complex<double>)*(2*nlag+1));           // final solution x*(A*(A^h*A)^-1)
  IPIV=(int*)malloc(sizeof(int) * N);
  for (m=0;m<nobs;m++) 
      {host_val[m].real((double)(random()/pow(2,31))-0.5);
       host_val[m].imag((double)(random()/pow(2,31))-0.5);
       host_code[m].real((double)(random()/pow(2,31))-0.5);
       host_code[m].imag((double)(random()/pow(2,31))-0.5);
      }
  printf("\n");
  for (m=0;m<nobs-5;m++)          // time shifted copies of the code
      {host_val[m+2]+=0.5*host_code[m]; // m+4 = ...
       host_val[m+5]+=1.2*host_code[m];
       host_val[m]+=0.3*host_code[m+3];
      }
#ifdef debug
  printf("a=[\n");
  for (m=0;m<nobs;m++) printf("(%.2lf)+j*(%.2lf) ",real(host_val[m]),imag(host_val[m]));
  printf("];\n");
#endif
  for (m=0;m<(2*nlag+1)*nobs;m++) host_mem[m]=0.;
  for (l=-nlag;l<=nlag;l++)
    for (m=0;m<nobs-(l+nlag);m++)
       if (l<0) host_mem[(m)+nobs*(l+nlag)]=host_code[m-l];
          else  host_mem[(m+l)+nobs*(l+nlag)]=host_code[m];
  // for (m=0;m<(2*nlag+1)*(2*nlag+1);m++) {host_Id[m].real(0.);host_Id[m].imag(0.);}
  for (m=0;m<(2*nlag+1)*(2*nlag+1);m++) host_Id[m]=0.;
  for (m=0;m<(2*nlag+1);m++) host_Id[m*(2*nlag+1)+m]=1.;
#ifdef debug
//  printf("b=[\n"); for (m=0;m<nlag*nobs;m++) printf("(%.2lf)+j*(%.2lf) ",real(host_mem[(m)]),imag(host_mem[(m)]));
//  printf("];\n");
  printf("b=[\n");
  affiche_matrice(host_mem,nobs,2*nlag+1);
  printf("];\n");
#endif
  gettimeofday(&tv1,NULL);
  cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, 2*nlag+1, 2*nlag+1, nobs, &alpha, host_mem, nobs, host_mem, nobs, &beta, host_res, 2*nlag+1 ); 
  // C := alpha*op( A )*op( B ) + beta*C,
  // alpha and beta are scalars, and A, B and C are matrices, with op( A )
  // an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  // (T/N,T/N,m,n,k,alpha, A,m si N ou k si T, B, k si N ou n si T, beta, C, m)
  LAPACKE_zgetrf(LAPACK_COL_MAJOR,N,N,reinterpret_cast <__complex__ double*>(host_res),N,IPIV); // LU decomposition: modifies input to output
  LAPACKE_zgetrs(LAPACK_COL_MAJOR,'N',N,N,reinterpret_cast <__complex__ double*>(host_res),N,IPIV,reinterpret_cast <__complex__ double*>(host_Id),N); // general inverse by solving A*X=B, here B=Id
#ifdef debug
  affiche_matrice(host_res,(nlag*2+1), 2*nlag+1);
#endif
  cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nobs, 2*nlag+1, 2*nlag+1, &alpha, host_mem, nobs, host_Id, 2*nlag+1, &beta, host_out, nobs ); 
#ifdef debug
  affiche_matrice(host_out,nobs, 2*nlag+1);
#endif
  cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, 1, 2*nlag+1, nobs, &alpha, host_val, nobs, host_out, nobs, &beta, host_final, 1 ); 
  gettimeofday(&tv2,NULL);
  printf("\ntime %ld\n",tv2.tv_usec-tv1.tv_usec);
  affiche_matrice(host_final,1, 2*nlag+1);
  printf("\n");
}
