all: demo_pinv demo_matrix

demo_matrix: demo_matrix.c
	gcc -Wall -o demo_matrix demo_matrix.c -lgsl

demo_pinv: demo_pinv.c
	gcc -Wall -o demo_pinv demo_pinv.c -lgsl

clean:
	rm demo_pinv demo_matrix
