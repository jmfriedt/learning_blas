all: demo1_matrix_square demo2_matrix_rectangleCM demo3_matrix_rectangleRM demo4_matrix_rectangleRMcomplex demo5_matrix_rectangleCM demo_matrix_final 

demo1_matrix_square: demo1_matrix_square.c
	gcc -Wall -o demo1_matrix_square demo1_matrix_square.c -L/usr/local/lib -lopenblas

demo2_matrix_rectangleCM: demo2_matrix_rectangleCM.c
	gcc -Wall -o demo2_matrix_rectangleCM demo2_matrix_rectangleCM.c -L/usr/local/lib -lopenblas

demo3_matrix_rectangleRM: demo3_matrix_rectangleRM.c
	gcc -Wall -o demo3_matrix_rectangleRM demo3_matrix_rectangleRM.c -L/usr/local/lib -lopenblas

demo4_matrix_rectangleRMcomplex:demo4_matrix_rectangleRMcomplex.cpp
	gcc -Wall -o demo4_matrix_rectangleRMcomplex demo4_matrix_rectangleRMcomplex.cpp -L/usr/local/lib -lopenblas

demo5_matrix_rectangleCM:demo5_matrix_rectangleCM.c
	gcc -Wall -o demo5_matrix_rectangleCM demo5_matrix_rectangleCM.c -L/usr/local/lib -lopenblas

demo_matrix_final: demo_matrix_final.cpp
	gcc -Wall -o demo_matrix_final demo_matrix_final.cpp -L/usr/local/lib -lopenblas -lm

clean: 
	rm demo1_matrix_square demo2_matrix_rectangleCM demo3_matrix_rectangleRM
