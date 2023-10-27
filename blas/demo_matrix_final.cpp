// g++ -o demo_matrix demo_matrix.cpp -L/usr/local/lib -lopenblas
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#include <cblas.h>
#include <complex>
#include <cstring>

int main(int argc, char *argv[])
{ int nobs=100;
  int nlag=20;
  int l,m;
  struct timeval tv1,tv2;
  std::complex<double> *host_mem,*host_code,*host_val,*host_res;
  std::complex<double> alpha,beta,pwr; 
  host_val=(std::complex<double>*)malloc(sizeof(std::complex<double>) * nobs);
  alpha=1.;beta=0.;
  for (m=0;m<nobs;m++)
    {host_val[m].real((double((m)%8)));
     host_val[m].imag((double((m)%9)));
    }
  gettimeofday(&tv1,NULL);
  pwr=cblas_zdotc(nobs, host_val, 1, host_val, 1);
  gettimeofday(&tv2,NULL); printf("\ntime %d\n",tv2.tv_usec-tv1.tv_usec);
  printf("power %lf\n",sqrt(pwr.real()*pwr.real()+pwr.imag()+pwr.imag()));


//https://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
  host_code=(std::complex<double>*)malloc(sizeof(std::complex<double>)*nobs);
  host_res=(std::complex<double>*)malloc(sizeof(std::complex<double>)*nlag);
  host_mem=(std::complex<double>*)malloc(sizeof(std::complex<double>)*nlag*nobs);
  for (m=0;m<nobs;m++) 
      {host_val[m].real((double)(random()/pow(2,31))-0.5);
       host_val[m].imag((double)(random()/pow(2,31))-0.5);
       host_code[m].real((double)(random()/pow(2,31))-0.5);
       host_code[m].imag((double)(random()/pow(2,31))-0.5);
      }
  for (m=0;m<nlag;m++) printf("%lf ",real(host_val[m]));
  printf("\n");
  for (m=0;m<nobs-12;m++)          // time shifted copies of the code
      {host_val[m+12]+=host_code[m];
       host_val[m+3]+=host_code[m];
       host_val[m]+=host_code[m+12];
       host_val[m]+=host_code[m+3];
      }

  memset(host_mem , 0x0, sizeof(std::complex<double>) * nobs * nlag);
  for (l=-nlag;l<nlag;l++)
    for (m=0;m<nobs-l;m++)
      host_mem[(l+nlag)*nobs+m+l+nlag]=host_code[m];
//correlation et position du max
  cblas_zgemm(CblasColMajor, CblasTrans, CblasNoTrans, nlag , nobs, nobs, &alpha, host_val, nobs, host_mem, nobs, &beta, host_res, nlag ); // CblasRowMajor ou CblasColMajor ? row-major (C) or column-major (Fortran) data ordering
// TODO : remplacer zgemm par zgemv
//  cudaMemcpy(host_res, dev_res, sizeof(cuDoubleComplex) * nlag       , cudaMemcpyDeviceToHost);
  // pk_idx -= 1;
  for (m=0;m<nlag;m++) printf("%lf ",abs(host_res[m]));
  printf("\n");
}

/*
        A = gsl_matrix_view_array(ci[i].dev_wav, ci[i].nobs, ci[i].nlag * 2 + 1);
        B = gsl_matrix_view_array(dev_obs, ci[i].nobs, (ci[i].bps - 1) * 2);
        C = gsl_matrix_view_array(ci[i].dev_res, ci[i].nlag * 2 + 1,(ci[i].bps - 1) * 2 );
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, alpha, &A.matrix, &B.matrix, beta, &C.matrix);
*/
//printf("after dgemm: %lf\n",ci[i].dev_res[0]);
