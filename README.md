## Learning BLAS

The objective of this set of programs is to become familiar with the API and compare implementation, 
performance and results
of linear algebra algorithms implemented using BLAS (CPU), cuBLAS (CUDA on GPU), and GSL, and
for the Fast Fourier Transform (FFT) using <a href="https://www.fftw.org/">fftw3</a> (CPU) and 
cuFFT (CUDA on GPU).

All programs have been tested on either an Intel Xeon W-2295 CPU @ 3.00GHz with its 36 cores, or a
T400 GPU with <a href="https://github.com/NVIDIA/cuda-samples">deviceQuery</a> stating
```
Device 0: "NVIDIA T400 4GB"
  CUDA Driver Version / Runtime Version          12.0 / 12.0
  CUDA Capability Major/Minor version number:    7.5
  Total amount of global memory:                 3897 MBytes (4086366208 bytes)
  (006) Multiprocessors, (064) CUDA Cores/MP:    384 CUDA Cores
...
```

All subdirectories include a Makefile for compiling the programs and a description of the expected output.
The focus is on least square optimization so aimed at pseudo-inverse matrix calculation, reaching this
result step by step from handling square matrices to rectangular matrices to inverting such matrices.
The FFT are computed in the context of correlations assessed as inverse Fourier transform of the product 
of the Fourier transform of the signal and the pattern.

Start with ``octave`` to get the basics on signal processing involved in this demonstration and
the expected outcome.

Then ``blas`` demonstrates the use of the linear algebra BLAS package to achieve the same result
in C++ on a general purpose CPU.

Then ``cuda`` demonstrates the use of cuBLAS and cuFFT for processing on an Nvidia GPU.

Finally ``gsl`` aims at wrapping BLAS function calls in interfaces removing some of the FORTRAN
flavor of the original library.
