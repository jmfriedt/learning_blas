all: demo_matrix demo_fft

demo_matrix: demo_matrix.cu
	nvcc demo_matrix.cu -o demo_matrix -lcublas -lm

demo_fft: demo_fft.cu
	nvcc demo_fft.cu -o demo_fft -lcufft -lm

clean:
	rm demo_matrix demo_fft
