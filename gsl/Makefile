all: demo_complex demo_matrix

demo_matrix: demo_matrix.c
	gcc -Wall -o demo_matrix demo_matrix.c -lgsl

demo_complex: demo_complex.c
	gcc -Wall -o demo_complex demo_complex.c -lgsl

clean:
	rm demo_complex demo_matrix
