#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

#include <sys/time.h>

#define PI 3.141592653589793

#undef debug

// https://stackoverflow.com/questions/15997888/creating-identity-matrix-with-cuda
__global__ void initIdentityGPU(cuDoubleComplex *devMatrix, int N) {
int x = blockDim.x*blockIdx.x + threadIdx.x;
int y = blockDim.y*blockIdx.y + threadIdx.y;
if(y < N && x < N) {
  if(x == y) {devMatrix[y*N+x].x = 1.; devMatrix[y*N+x].y = 0.;}
    else {devMatrix[y*N+x].x = 0.; devMatrix[y*N+x].x = 0.;}
  }
}

// https://stackoverflow.com/questions/22887167/cublas-incorrect-inversion-for-matrix-with-zero-pivot
int main()
{
  int nobs=2100;
  int nlag=15; // MUST BE < 16 !
  int l,m;
  const int N=2*nlag+1;
  cuDoubleComplex *dev_mem, *dev_res, *dev_val, *dev_inv, *dev_Id, *dev_in;     // .x, .y
  std::complex<double> *host_mem,*host_res,*host_val,*host_code;  // real(), imag()
  cublasHandle_t handle;
  cuDoubleComplex alpha,beta;
  alpha.x=1.;alpha.y=0.;
  beta.x=0;beta.y=0;
  
  cudaSetDevice(0);
  cublasCreate(&handle);

  cudaMalloc((void **)&dev_mem, sizeof(cuDoubleComplex) * nobs * (nlag*2+1));
  host_mem=(std::complex<double>*)malloc(sizeof(std::complex<double>)*nobs*(nlag*2+1));
  host_val=(std::complex<double>*)malloc(sizeof(std::complex<double>)*nobs);
  host_code=(std::complex<double>*)malloc(sizeof(std::complex<double>)*nobs);
  host_res=(std::complex<double>*)malloc(sizeof(std::complex<double>)*(2*nlag+1)*(2*nlag+1));
  cudaMalloc((void **)&dev_res, sizeof(cuDoubleComplex) * (2*nlag+1));
  cudaMalloc((void **)&dev_val, sizeof(cuDoubleComplex) * nobs);
  cudaMalloc((void **)&dev_inv, sizeof(cuDoubleComplex) * (2*nlag+1)*(2*nlag+1));
  cudaMalloc((void **)&dev_in, sizeof(cuDoubleComplex) * (2*nlag+1)*(2*nlag+1));
  cudaMalloc((void **)&dev_Id, sizeof(cuDoubleComplex) * (2*nlag+1)*(2*nlag+1));
  for (m=0;m<nobs;m++) 
      {host_val[m].real((double)(random()/pow(2,31))-0.5);
       host_val[m].imag((double)(random()/pow(2,31))-0.5);
       host_code[m].real((double)(random()/pow(2,31))-0.5);
       host_code[m].imag((double)(random()/pow(2,31))-0.5);
      }
//  for (m=0;m<nlag;m++) printf("%lf ",real(host_val[m]));
//  printf("\n");
  for (m=0;m<nobs-12;m++)          // time shifted copies of the code
      {host_val[m+10]+=host_code[m]*.3;
       host_val[m+5]+=host_code[m]*.5;
       host_val[m]+=host_code[m+12]*.7;
       host_val[m]+=host_code[m+3];
      }

  memset(host_mem , 0x0, sizeof(std::complex<double>) * nobs * (2*nlag+1));
  for (l=-nlag;l<=nlag;l++)
    for (m=0;m<nobs-(l+nlag);m++)
       if (l<0) host_mem[(m)+nobs*(l+nlag)]=host_code[m-l];
          else  host_mem[(m+l)+nobs*(l+nlag)]=host_code[m];
      // host_mem[(l+nlag)*nobs+m+l+nlag]=host_code[m];
  cublasSetMatrix (nobs, nlag*2+1, sizeof(*host_mem), host_mem, nobs, dev_mem, nobs);
  cublasSetMatrix (1, nobs, sizeof(*host_val), host_val, 1   , dev_val, 1   );
  if (cublasZgemm(handle, CUBLAS_OP_C, CUBLAS_OP_N, N, N, nobs, &alpha, dev_mem,  nobs, dev_mem, nobs, &beta, dev_in, N) != CUBLAS_STATUS_SUCCESS)
//	 printf("error 0\n");
  cudaDeviceSynchronize();
  int *P, *INFO;
  cudaMalloc((void **)&P, sizeof(int) * (2*nlag+1));
  cudaMalloc((void **)&INFO, sizeof(int));
//  (cudaMalloc<int>(&P,N * sizeof(int)));
//  (cudaMalloc<int>(&INFO,sizeof(int)));
/*
  cublasZgetrfBatched(handle, N, A_d, N, P, INFO,1)
  cublasZgetriBatched(handle, N, A_d, N, P, C_d, N , INFO, 1)
*/

  cusolverDnHandle_t handlegetrs = NULL;
  int bufferSize = 0;
  cuDoubleComplex *buffer = NULL;
//  initIdentityGPU<<<128, 128>>>(dev_Id,N);
  for (m=0;m<(2*nlag+1);m++)
    for (l=0;l<(2*nlag+1);l++) 
      if (m!=l) {host_res[m+l*(2*nlag+1)].real(0.);host_res[m+l*(2*nlag+1)].imag(0.);}
      else {host_res[m+l*(2*nlag+1)].real(1.);host_res[m+l*(2*nlag+1)].imag(0.);}
  cudaMemcpy(dev_Id,host_res,sizeof(cuDoubleComplex) * (2*nlag+1)*(2*nlag+1),cudaMemcpyHostToDevice);
  memset(host_res , 0x0, sizeof(std::complex<double>) * (2*nlag+1) * (2*nlag+1));
#ifdef debug
  cudaMemcpy(host_res,dev_in,sizeof(cuDoubleComplex) * (2*nlag+1)*(2*nlag+1),cudaMemcpyDeviceToHost);
  printf("Id\n");
  for (m=0;m<(2*nlag+1);m++)
    {for (l=0;l<(2*nlag+1);l++) printf("%.0lf ",real(host_res[m+l*(2*nlag+1)]));
     printf("; \n");
    }
#endif
  cusolverDnCreate(&handlegetrs);
  cusolverDnZgetrf_bufferSize(handlegetrs, N, N, dev_in, N, &bufferSize);
  cudaMalloc(&buffer, sizeof(cuDoubleComplex) * bufferSize );
// https://docs.nvidia.com/cuda/cusolver/index.html
  if (cusolverDnZgetrf(handlegetrs, N, N, dev_in, N, buffer, P, INFO) != CUSOLVER_STATUS_SUCCESS)
     printf("error 1\n");
  if (cusolverDnZgetrs(handlegetrs, CUBLAS_OP_N, N, N, dev_in, N, P, dev_Id, N, INFO) != CUSOLVER_STATUS_SUCCESS)
     printf("error 2\n");
  cudaMemcpy(host_res,dev_Id,sizeof(cuDoubleComplex) * (2*nlag+1)*(2*nlag+1),cudaMemcpyDeviceToHost);
#ifdef debug
  printf("res\n");
  for (m=0;m<(2*nlag+1);m++)
    {for (l=0;l<(2*nlag+1);l++) printf("%.4lf ",real(host_res[m+l*(2*nlag+1)]));
     printf("; \n");
    }
#endif
  if (cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nobs, 2*nlag+1, 2*nlag+1, &alpha, dev_mem,  nobs, dev_Id, 2*nlag+1, &beta, dev_mem, nobs) != CUBLAS_STATUS_SUCCESS)
     printf("error 3\n");
  cublasZgemm(handle, CUBLAS_OP_C, CUBLAS_OP_N, 1, 2*nlag+1, nobs, &alpha, dev_val,  nobs, dev_mem, nobs, &beta, dev_res, 1);
  cudaMemcpy(host_res,dev_res,sizeof(cuDoubleComplex) * (2*nlag+1),cudaMemcpyDeviceToHost);
  for (m=0;m<2*nlag+1;m++) printf("%.2lf ",abs(host_res[m]));
  printf("\n");
  cudaFree(P), cudaFree(INFO), cublasDestroy(handle);
}

/*
#include <lapacke.h>
  cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, 2*nlag+1, 2*nlag+1, nobs, &alpha, host_mem, nobs, host_mem, nobs, &beta, host_res, 2*nlag+1 ); 
  zgetrf_(&N,&N,reinterpret_cast <__complex__ double*>(host_res),&N,IPIV,&info); // LU decomposition: modifies input to output
  zgetri_(&N,reinterpret_cast <__complex__ double*>(host_res),&N,IPIV,reinterpret_cast <__complex__ double*>(WORK),&LWORK,&info); // inverse
  cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nobs, 2*nlag+1, 2*nlag+1, &alpha, host_mem, nobs, host_res, 2*nlag+1, &beta, host_out, nobs ); 
  cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, 1, 2*nlag+1, nobs, &alpha, host_val, nobs, host_out, nobs, &beta, host_final, 1 ); */
