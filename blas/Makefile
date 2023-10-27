all: demo1_matrix_square demo2_matrix_rectangle

demo1_matrix_square: demo1_matrix_square.c
	gcc -Wall -o demo1_matrix_square demo1_matrix_square.c -L/usr/local/lib -lopenblas

demo2_matrix_rectangle: demo2_matrix_rectangle.c
	gcc -Wall -o demo2_matrix_rectangle demo2_matrix_rectangle.c -L/usr/local/lib -lopenblas

clean: 
	rm demo1_matrix_square demo2_matrix_rectangle
