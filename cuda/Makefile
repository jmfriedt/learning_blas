all: demo_matrix demo_fft demo_pinv

CFLAGS=-O2 -arch=sm_50

demo_matrix: demo_matrix.cu
	nvcc demo_matrix.cu $(CFLAGS) -o demo_matrix -lcublas -lm

demo_pinv: demo_pinv.cu
	nvcc demo_pinv.cu $(CFLAGS) -o demo_pinv -lcublas  -lcusolver

demo_fft: demo_fft.cu
	nvcc demo_fft.cu $(CFLAGS) -o demo_fft -lcufft -lm -lfftw3

clean:
	rm demo_matrix demo_fft demo_pinv
