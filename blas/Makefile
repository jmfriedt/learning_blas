all: demo1_matrix_square demo2_matrix_rectangleCM demo3_matrix_rectangleRM demo4_matrix_rectangleRMcomplex demo5_matrix_rectangleCM demo6_matrix_xcorr demo7_matrix_inv demo8_matrix_pinv demo9_matrix_pinv_zgetrs

FLAGS=-I/usr/include/clapack -lopenblas -llapacke64

demo1_matrix_square: demo1_matrix_square.c
	gcc -Wall -o demo1_matrix_square demo1_matrix_square.c $(FLAGS)

demo2_matrix_rectangleCM: demo2_matrix_rectangleCM.c
	gcc -Wall -o demo2_matrix_rectangleCM demo2_matrix_rectangleCM.c $(FLAGS)

demo3_matrix_rectangleRM: demo3_matrix_rectangleRM.c
	gcc -Wall -o demo3_matrix_rectangleRM demo3_matrix_rectangleRM.c $(FLAGS)

demo4_matrix_rectangleRMcomplex:demo4_matrix_rectangleRMcomplex.cpp
	gcc -Wall -o demo4_matrix_rectangleRMcomplex demo4_matrix_rectangleRMcomplex.cpp $(FLAGS)

demo5_matrix_rectangleCM:demo5_matrix_rectangleCM.c
	gcc -Wall -o demo5_matrix_rectangleCM demo5_matrix_rectangleCM.c $(FLAGS)

demo6_matrix_xcorr: demo6_matrix_xcorr.cpp
	gcc -Wall -o demo6_matrix_xcorr demo6_matrix_xcorr.cpp $(FLAGS) -lm

demo7_matrix_inv: demo7_matrix_inv.cpp
	gcc -Wall -o demo7_matrix_inv demo7_matrix_inv.cpp $(FLAGS) -lm

demo8_matrix_pinv: demo8_matrix_pinv.cpp
	gcc -Wall -o demo8_matrix_pinv demo8_matrix_pinv.cpp $(FLAGS) -lm

demo9_matrix_pinv_zgetrs: demo9_matrix_pinv_zgetrs.cpp
	gcc -Wall -o demo9_matrix_pinv_zgetrs demo9_matrix_pinv_zgetrs.cpp $(FLAGS) -lm

clean: 
	rm demo1_matrix_square demo2_matrix_rectangleCM demo3_matrix_rectangleRM demo4_matrix_rectangleRMcomplex demo5_matrix_rectangleCM demo6_matrix_xcorr demo7_matrix_inv demo8_matrix_pinv demo9_matrix_pinv_zgetrs
