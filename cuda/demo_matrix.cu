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
  int nobs=2100;
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
  cudaMalloc((void **)&dev_mem, sizeof(cuDoubleComplex) * nobs * (nlag*2+1));
  host_mem=(std::complex<double>*)malloc(sizeof(std::complex<double>)*nobs*(nlag*2+1));
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
  host_res=(std::complex<double>*)malloc(sizeof(std::complex<double>)*(2*nlag+1));
  cudaMalloc((void **)&dev_res, sizeof(cuDoubleComplex) * (2*nlag+1));
  cudaMalloc((void **)&dev_val, sizeof(cuDoubleComplex) * nobs);
  for (m=0;m<nobs;m++) 
      {host_val[m].real((double)(random()/pow(2,31))-0.5);
       host_val[m].imag((double)(random()/pow(2,31))-0.5);
       host_code[m].real((double)(random()/pow(2,31))-0.5);
       host_code[m].imag((double)(random()/pow(2,31))-0.5);
      }
//  for (m=0;m<nlag;m++) printf("%lf ",real(host_val[m]));
//  printf("\n");
  for (m=0;m<nobs-12;m++)          // time shifted copies of the code
      {host_val[m+10]+=host_code[m];
       host_val[m+5]+=host_code[m];
       host_val[m]+=host_code[m+12];
       host_val[m]+=host_code[m+3];
      }

  memset(host_mem , 0x0, sizeof(std::complex<double>) * nobs * nlag);
  for (l=-nlag;l<=nlag;l++)
    for (m=0;m<nobs-(l+nlag);m++)
       if (l<0) host_mem[(m)+nobs*(l+nlag)]=host_code[m-l];
          else  host_mem[(m+l)+nobs*(l+nlag)]=host_code[m];
      // host_mem[(l+nlag)*nobs+m+l+nlag]=host_code[m];
  cublasSetMatrix (nobs, nlag*2+1, sizeof(*host_mem), host_mem, nobs, dev_mem, nobs);
  cublasSetMatrix (1, nobs, sizeof(*host_val), host_val, 1   , dev_val, 1   );
//  cudaMemcpy(dev_mem, host_mem, sizeof(cuDoubleComplex) * nobs * nlag, cudaMemcpyHostToDevice);
//  cudaMemcpy(dev_val, host_val, sizeof(cuDoubleComplex) * nobs       , cudaMemcpyHostToDevice);
//correlation et position du max
  cublasZgemm(handle, CUBLAS_OP_C, CUBLAS_OP_N, 1, 2*nlag+1, nobs, &alpha, dev_val,  nobs, dev_mem, nobs, &beta, dev_res, 1);
//   ** On entry to ZGEMM  parameter number 8 had an illegal value si 1 au lieu de nobs apres dev_val
//  cudaMemcpy(host_res, dev_res, sizeof(cuDoubleComplex) * nlag       , cudaMemcpyDeviceToHost);
  cublasGetMatrix (1, 2*nlag+1, sizeof(*host_res), dev_res, 1, host_res, 1);
// cublasIdamax(handle, ci[i].nlag * 2 + 1, ci[i].dev_cor + p * (ci[i].nlag * 2 + 1), 1, &pk_idx);
// pk_idx -= 1;
  for (m=0;m<2*nlag+1;m++) printf("%.2lf ",abs(host_res[m]));
  printf("\n");
}
