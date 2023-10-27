//nvcc demo_matrix.cu -o demo_matrix -lcublas -lm
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <cublas_v2.h>

#include <gsl/gsl_fit.h>

#include <sys/time.h>

#define PI 3.141592653589793

int main()
{
  //cublasHandle_t handle;
  int nobs=100;
  int nlag=20;
  int l,m;
  struct timeval tv1,tv2;
  cuDoubleComplex *dev_mem, *dev_res, *dev_val;     // .x, .y
  std::complex<double> *host_mem,*host_res,*host_val,*host_code;  // real(), imag()
  cublasHandle_t handle;
  cuDoubleComplex pwr;
  cuDoubleComplex alpha,beta;
  alpha.x=1.;alpha.y=0.;
  beta.x=0;beta.y=0;

  cudaSetDevice(0);
  cublasCreate(&handle);
  cudaMalloc((void **)&dev_mem, sizeof(cuDoubleComplex) * nobs * nlag*2);
  host_mem=(std::complex<double>*)malloc(sizeof(std::complex<double>)*nobs*nlag*2);
  for (m=0;m<nobs;m++)
    {host_mem[m].real((double((m)%8)));
     host_mem[m].imag((double((m)%9)));
    }
  cudaMemcpy(dev_mem, host_mem, sizeof(cuDoubleComplex) * nobs , cudaMemcpyHostToDevice);
  gettimeofday(&tv1,NULL);
  cublasZdotc(handle, nobs, dev_mem, 1, dev_mem, 1, &pwr);
  gettimeofday(&tv2,NULL); printf("\ntime %d\n",tv2.tv_usec-tv1.tv_usec);
  printf("power %lf\n",sqrt(pwr.x*pwr.x+pwr.y+pwr.y));

/*
octave:9> a=mod([0:99],8)+j*mod([0:99],9);
octave:10> (a*a')^2
3938
aussi sum(abs(a).^2)^2
power 3938.000000
*/

//https://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
  host_val=(std::complex<double>*)malloc(sizeof(std::complex<double>)*nobs);
  host_code=(std::complex<double>*)malloc(sizeof(std::complex<double>)*nobs);
  host_res=(std::complex<double>*)malloc(sizeof(std::complex<double>)*nlag);
  cudaMalloc((void **)&dev_res, sizeof(cuDoubleComplex) * nlag);
  cudaMalloc((void **)&dev_val, sizeof(cuDoubleComplex) * nobs);
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
  cublasSetMatrix (nobs, nlag*2, sizeof(*host_mem), host_mem, nobs, dev_mem, nobs);
  cublasSetMatrix (1, nobs, sizeof(*host_val), host_val, 1   , dev_val, 1   );
//  cudaMemcpy(dev_mem, host_mem, sizeof(cuDoubleComplex) * nobs * nlag, cudaMemcpyHostToDevice);
//  cudaMemcpy(dev_val, host_val, sizeof(cuDoubleComplex) * nobs       , cudaMemcpyHostToDevice);
//correlation et position du max
  //cublasZgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 1, nlag, nobs, &alpha, dev_val, nobs, dev_mem, nobs, &beta, dev_res, 1);
  //                  transpose    no transpose m   n      k     1.     mxk       k    kxn     k      0.     res    m
  cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 1, nlag, nobs, &alpha, dev_val,   1 , dev_mem, nobs, &beta, dev_res, 1);
  //                  transpose    no transpose m   n      k     1.     mxk       m    kxn     k      0.     res    m
  //cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, 1, nlag, nobs, &alpha, dev_val,   1 , dev_mem, nlag, &beta, dev_res, 1);
  //                  transpose    no transpose m   n      k     1.     mxk       m    kxn     n      0.     res    m
  //cublasZgemm(handle, CUBLAS_OP_T, CUBLAS_OP_T, 1, nlag, nobs, &alpha, dev_val, nobs, dev_mem, nlag, &beta, dev_res, 1);
  //                  transpose    no transpose m   n      k     1.     mxk       k    kxn     n      0.     res    m
     // C := alpha*op( A )*op( B ) + beta*C,
     // alpha and beta are scalars, and A, B and C are matrices, with op( A )
     // an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
     // (T/N,T/N,m,n,k,alpha, A,m/k selon N ou T, B, k/n selon N ou T, beta, C, m)
// ** On entry to ZGEMM  parameter number 8 had an illegal value
//  cudaMemcpy(host_res, dev_res, sizeof(cuDoubleComplex) * nlag       , cudaMemcpyDeviceToHost);
  cublasGetMatrix (1, nlag, sizeof(*host_res), dev_res, 1, host_res, 1);
  // cublasIdamax(handle, ci[i].nlag * 2 + 1, ci[i].dev_cor + p * (ci[i].nlag * 2 + 1), 1, &pk_idx);
  // pk_idx -= 1;
  for (m=0;m<nlag;m++) printf("%lf ",abs(host_res[m]));
  printf("\n");
}


