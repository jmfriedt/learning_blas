# BLAS examples

Executing the C program ``demo1_matrix_square`` will provide the following outputs (left) and right the GNU/Octave output
```
       C                    Octave
0.00 2.00 -1.00 1.00     a=[0 -1 ; 2 1]
0.00 2.00 -1.00 1.00 
CM a*a                   a*a
-2.00 -1.00              -2  -1
2.00 -1.00                2  -1

CM a'*a'                 a'*a'
-2.00 2.00               -2   2
-1.00 -1.00              -1  -1

CM a*a'                  a*a'
1.00 -1.00                1  -1
-1.00 5.00               -1   5

CM a'*a                  a'*a
4.00 2.00                4   2
2.00 2.00                2   2
                         a=[0 2 ; -1 1]
RM a*a                   the display function must be
-2.00 -1.00              changed to Rwo Major to match
2.00 -1.00               the Octave output

RM a'*a'
-2.00 2.00 
-1.00 -1.00 

RM a*a'
4.00 2.00 
2.00 2.00 

RM a'*a
1.00 -1.00 
-1.00 5.00 
```

Executing the C program ``demo2_matrix_rectangleCM`` will provide the following outputs (left) and right the GNU/Octave output
```
             C                       Octave
0.00 2.00 -1.00 1.00 -2.00 0.00   a=[0 -1;-2 2 ;1 0]
0.00 2.00 -1.00 1.00 -2.00 0.00   a'*a
5.00 -4.00                        5  -4
-4.00 5.00                       -4   5
                                  a*a'
1.00 -2.00 0.00                   1  -2   0
-2.00 8.00 -2.00                 -2   8  -2
0.00 -2.00 1.00                   0  -2   1
```

Executing the C program ``demo3_matrix_rectangleRM`` will provide the following outputs (left) and right the GNU/Octave output
```
             C                      Octave
0.00 -1.00 -2.00 2.00 1.00 0.00   a=[0 -1 -2; 2 1 0]
                                  a*a'
5.00 -1.00                        5  -1
-1.00 5.00                       -1   5
                                 a'*a
4.00 2.00 0.00                   4   2   0
2.00 2.00 2.00                   2   2   2
0.00 2.00 4.00                   0   2   4
```
