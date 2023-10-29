// g++ -o demo_matrix demo_matrix.cpp -L/usr/local/lib -lopenblas
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#include <cblas.h>
#include <lapacke.h>
#include <complex>

#define debug

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
{ int nobs=5;
  int nlag=5;
  const int N=nobs;
  int l,m,info;
  int *IPIV;
  int LWORK = N*N;
  std::complex<double> *WORK;

  std::complex<double> *host_mem;
  std::complex<double> alpha,beta,pwr; 
  host_mem=(std::complex<double>*)malloc(sizeof(std::complex<double>) * nobs * nlag);
  WORK=(std::complex<double>*)malloc(sizeof(std::complex<double>) * LWORK);
  IPIV=(int*)malloc(sizeof(int) * N);

  for (l=0;l<nlag;l++)
    for (m=0;m<nobs;m++)
      {host_mem[m+nobs*l].real((double)(random()/pow(2,31))-0.5);
       host_mem[m+nobs*l].imag((double)(random()/pow(2,31))-0.5);
      }
#ifdef debug
  printf("b=[\n");
  affiche_matrice(host_mem,nobs,nlag);
  printf("];\n");
#endif
  // LAPACKE_zgetrf(LAPACK_COL_MAJOR,N,N,reinterpret_cast <__complex__ double*>(host_mem),N,IPIV);
  zgetrf_(&N,&N,reinterpret_cast <__complex__ double*>(host_mem),&N,IPIV,&info); // LU decomposition: modifies input to output
  if (info!=0) printf("error: info = %d\n",info);
  // LAPACKE_zgetri(LAPACK_COL_MAJOR,N,reinterpret_cast <__complex__ double*>(host_mem),N,IPIV);
  zgetri_(&N,reinterpret_cast <__complex__ double*>(host_mem),&N,IPIV,reinterpret_cast <__complex__ double*>(WORK),&LWORK,&info); // inverse
  if (info!=0) printf("error: info = %d\n",info);
  affiche_matrice(host_mem,nobs, nlag);
}

/*
> inv(b)
ans =
   0.880897 + 1.127919i   0.056371 - 0.677401i   0.228673 - 0.187015i  -0.799350 - 1.048781i   0.446246 - 0.215708i
  -1.108056 + 1.096574i  -0.116187 - 0.444804i   1.059831 - 0.955830i   0.922420 - 0.811728i   0.222927 + 0.039027i
  -0.668254 + 0.643562i  -0.484864 - 0.641085i   0.214633 - 0.301413i  -0.465733 - 1.050222i   1.246947 - 0.496074i
  -0.045621 - 0.534575i  -0.235588 + 0.674173i   0.093298 - 0.210545i  -1.461550 - 1.171325i   1.189031 - 1.377971i
   0.152451 - 0.088931i  -0.392460 - 0.610878i  -0.714504 + 0.109109i   0.054098 + 0.133626i  -0.731174 - 0.550664i
*/
