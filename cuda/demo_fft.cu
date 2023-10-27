// FFT et iFFT

// nvcc demo_fft.cu -o demo_fft -lcufft -lm
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <cufft.h>

#include <fftw3.h> // nvcc demo_fft.cu -o demo_fft -lcufft -lm -lfftw3

#include <sys/time.h>

#define PI 3.141592653589793

int main()
{
  //cublasHandle_t handle;
  int nobs=10000000;
  int k;
  double f=440.,fs=48000.,ph=0.;
  cufftDoubleComplex *dev_mem;     // .x, .y
  std::complex<double> *host_mem;  // real(), imag()
  struct timeval tv1,tv2;

  cufftHandle plan;
  cudaSetDevice(0);
  //cublasCreate(&handle);
  cudaMalloc((void **)&dev_mem, sizeof(cufftDoubleComplex) * nobs);
  host_mem=(std::complex<double>*)malloc(sizeof(std::complex<double>)*nobs);
  for (k=0;k<nobs;k++)
  {host_mem[k].real(cos(ph));
   host_mem[k].imag(sin(ph));
   ph+=2*PI*f/fs;if (ph>2*PI) ph-=2*PI;
   if (k<20) printf("%.2lf ",real(host_mem[k]));
  }
  printf("\n");
  cudaMemcpy(dev_mem, host_mem, sizeof(cufftDoubleComplex) * nobs, cudaMemcpyHostToDevice);
  gettimeofday(&tv1,NULL);
  cufftPlan1d(&plan, nobs, CUFFT_Z2Z, 1);
  gettimeofday(&tv2,NULL); printf("\ntime %d\n",tv2.tv_usec-tv1.tv_usec);
  gettimeofday(&tv1,NULL);
  cufftExecZ2Z(plan, dev_mem, dev_mem, CUFFT_FORWARD);
  gettimeofday(&tv2,NULL); printf("\ntime %d\n",tv2.tv_usec-tv1.tv_usec);
  cudaMemcpy(host_mem, dev_mem, sizeof(cufftDoubleComplex) * nobs, cudaMemcpyDeviceToHost);
  for (k=0;k<20;k++) printf("%.2lf ",abs(host_mem[k+(int)(f/fs*(float)nobs)-10])); // 440/fs*nobs
  gettimeofday(&tv1,NULL);
  cufftExecZ2Z(plan, dev_mem, dev_mem, CUFFT_INVERSE);
  gettimeofday(&tv2,NULL); printf("\ntime %d\n",tv2.tv_usec-tv1.tv_usec);
  cufftDestroy(plan);
  
  printf("\n");
  cudaMemcpy(host_mem, dev_mem, sizeof(cufftDoubleComplex) * nobs, cudaMemcpyDeviceToHost);
  for (k=0;k<20;k++) printf("%.2lf ",real(host_mem[k])/(double)nobs);
  printf("\n");

  ph=0.;
  for (k=0;k<nobs;k++)
  {host_mem[k].real(cos(ph));
   host_mem[k].imag(sin(ph));
   ph+=2*PI*f/fs;if (ph>2*PI) ph-=2*PI;
   if (k<20) printf("%.2lf ",real(host_mem[k]));
  }
  printf("\n");
  fftw_plan _plan_a_dx;
  fftw_plan _ifft_dx;
  _plan_a_dx = fftw_plan_dft_1d(nobs,
     reinterpret_cast<fftw_complex*>(host_mem), reinterpret_cast<fftw_complex*>(host_mem),
     FFTW_FORWARD, FFTW_ESTIMATE);
  _ifft_dx = fftw_plan_dft_1d(nobs,
     reinterpret_cast<fftw_complex*>(host_mem), reinterpret_cast<fftw_complex*>(host_mem),
     FFTW_BACKWARD, FFTW_ESTIMATE);
  gettimeofday(&tv1,NULL);
  fftw_execute(_plan_a_dx);
  gettimeofday(&tv2,NULL); printf("\ntime %d\n",tv2.tv_usec-tv1.tv_usec);
  for (k=0;k<20;k++) printf("%.2lf ",abs(host_mem[k+(int)(f/fs*(float)nobs)-10])); // 440/fs*nobs
  gettimeofday(&tv1,NULL);
  fftw_execute(_ifft_dx);
  gettimeofday(&tv2,NULL); printf("\ntime %d\n",tv2.tv_usec-tv1.tv_usec);
  for (k=0;k<20;k++) printf("%.2lf ",real(host_mem[k])/(double)nobs);
}
