all: demo_matrix demo_fft demo_pinv

demo_matrix: demo_matrix.cu
	nvcc demo_matrix.cu -o demo_matrix -lcublas -lm

demo_pinv: demo_pinv.cu
	nvcc demo_pinv.cu -o demo_pinv -lcublas  -lcusolver

demo_fft: demo_fft.cu
	nvcc demo_fft.cu -o demo_fft -lcufft -lm -lfftw3

clean:
	rm demo_matrix demo_fft demo_pinv
