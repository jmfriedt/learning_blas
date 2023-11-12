all: demo_pinv demo_matrix

CFLAGS=-Wall -O2 

demo_matrix: demo_matrix.c
	gcc $(CFLAGS) -o demo_matrix demo_matrix.c -lgsl

demo_pinv: demo_pinv.c
	gcc $(CFLAGS) -o demo_pinv demo_pinv.c -lgsl

clean:
	rm demo_pinv demo_matrix
