// g++ -o demo_matrix demo_matrix.cpp -L/usr/local/lib -lopenblas
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#include <cblas.h>
#include <complex>

#undef debug

void affiche_matrice(std::complex<double>*mat,int x,int y)
{int l,m;
  for (m=0;m<x;m++) 
   {for (l=0;l<y;l++)
       printf("%.2lf ",abs(mat[l*x+m]));
       //printf("(%.5lf)+j*(%.5lf) ",real(mat[l*x+m]),imag(mat[l*x+m]));
    printf(" ;\n");
   }
  printf("\n");
}

int main(int argc, char *argv[])
{ int nobs=20000;
  int nlag=15;
  int l,m;
  struct timeval tv1,tv2;
  std::complex<double> *host_mem,*host_code,*host_val,*host_res;
  std::complex<double> alpha,beta,pwr; 
  host_val=(std::complex<double>*)malloc(sizeof(std::complex<double>) * nobs);
  alpha=1.;beta=0.;
  for (m=0;m<nobs;m++)
    {host_val[m].real((double((m)%8)));
     host_val[m].imag((double((m)%9)));
    }
  gettimeofday(&tv1,NULL);
  pwr=cblas_zdotc(nobs, host_val, 1, host_val, 1);
  gettimeofday(&tv2,NULL); printf("\ntime %ld\n",tv2.tv_usec-tv1.tv_usec);
  printf("power %lf\n",sqrt(pwr.real()*pwr.real()+pwr.imag()+pwr.imag()));

//https://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
  host_code=(std::complex<double>*)malloc(sizeof(std::complex<double>)*nobs);
  host_res=(std::complex<double>*)malloc(sizeof(std::complex<double>)*(nlag*2+1));
  host_mem=(std::complex<double>*)malloc(sizeof(std::complex<double>)*(2*nlag+1)*nobs);
  for (m=0;m<nobs;m++) 
      {host_val[m].real(0); // (double)(random()/pow(2,31))-0.5);
       host_val[m].imag(0); //(double)(random()/pow(2,31))-0.5);
       host_code[m].real((double)(random()/pow(2,31))-0.5);
       host_code[m].imag((double)(random()/pow(2,31))-0.5);
      }
  printf("\n");
  for (m=0;m<nobs-5;m++)          // time shifted copies of the code
      {host_val[m+2]+=0.5*host_code[m]; // m+4 = ...
       host_val[m+5]+=1.2*host_code[m];
       host_val[m]+=0.3*host_code[m+3];
      }
  for (m=0;m<nlag;m++)          // time shifted copies of the code
    {pwr=cblas_zdotc(nobs, host_val+m, 1, host_code, 1);
     printf("power%d %lf\n",m,abs(pwr));
    }
  
#ifdef debug
  printf("a=[\n");
  for (m=0;m<nobs;m++) printf("(%.2lf)+j*(%.2lf) ",real(host_val[m]),imag(host_val[m]));
  printf("];\n");
#endif
  for (m=0;m<(2*nlag+1)*nobs;m++) host_mem[m]=0.;
  for (l=-nlag;l<=nlag;l++)
    for (m=0;m<nobs-(l+nlag);m++)
       if (l<0) host_mem[(m)+nobs*(l+nlag)]=host_code[m-l];
          else  host_mem[(m+l)+nobs*(l+nlag)]=host_code[m];
//correlation et position du max
#ifdef debug
//  printf("b=[\n"); for (m=0;m<nlag*nobs;m++) printf("(%.2lf)+j*(%.2lf) ",real(host_mem[(m)]),imag(host_mem[(m)]));
//  printf("];\n");
  printf("b=[\n");
  affiche_matrice(host_mem,nobs,2*nlag+1);
  printf("];\n");
#endif
  // cblas_zgemm(CblasColMajor, CblasTrans, CblasNoTrans, 1 , 2*nlag+1, nobs, &alpha, host_val, nobs, host_mem, nobs, &beta, host_res, 1 ); 
  cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, 1 , 2*nlag+1, nobs, &alpha, host_val, nobs, host_mem, nobs, &beta, host_res, 1 ); 
  affiche_matrice(host_res,1, 2*nlag+1);
  // C := alpha*op( A )*op( B ) + beta*C,
  // alpha and beta are scalars, and A, B and C are matrices, with op( A )
  // an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  // (T/N,T/N,m,n,k,alpha, A,m si N ou k si T, B, k si N ou n si T, beta, C, m)
//identical result:
  cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, 2*nlag+1 ,1, nobs, &alpha, host_mem, nobs, host_val, nobs, &beta, host_res, 2*nlag+1 );
  affiche_matrice(host_res,1,2*nlag+1);
  cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, 1, 2*nlag+1, nobs, &alpha, host_val, nobs, host_mem, nobs, &beta, host_res, 1 );
  affiche_matrice(host_res,1,2*nlag+1);
  cblas_zgemv(CblasColMajor, CblasConjTrans, nobs ,2*nlag+1, &alpha, host_mem, nobs, host_val, 1, &beta, host_res, 1 );
  affiche_matrice(host_res,1,2*nlag+1);
  printf("\n");
}

// cudaMemcpy(host_res, dev_res, sizeof(cuDoubleComplex) * nlag , cudaMemcpyDeviceToHost); +  pk_idx -= 1;

/*
time 16
power 703.000000

power0 1.119355
power1 0.945368
power2 1.690374
power3 0.879536
power4 0.835104
power5 2.346042
power6 1.226513
a=[
(-0.16)+j*(0.27) (-0.22)+j*(0.05) (0.32)+j*(0.02) (0.15)+j*(0.31) (0.86)+j*(0.11) (0.31)+j*(0.38) (-0.30)+j*(0.46) (-0.09)+j*(-0.43) (-0.66)+j*(0.59) (-0.11)+j*(0.37) (-0.26)+j*(-0.05) (0.01)+j*(-0.16) (-0.02)+j*(0.50) (-0.11)+j*(0.32) (-0.56)+j*(0.03) (-0.85)+j*(-0.65) (0.14)+j*(0.02) (-0.34)+j*(-0.10) (-0.37)+j*(-0.39) (0.50)+j*(-0.28) ];
b=[
0.61 0.14 0.13 0.23 0.31 0.51 0.41 0.36 0.00 0.00 0.00 0.00 0.00 0.00 0.00  ;
0.26 0.61 0.14 0.13 0.23 0.31 0.51 0.41 0.36 0.00 0.00 0.00 0.00 0.00 0.00  ;
0.37 0.26 0.61 0.14 0.13 0.23 0.31 0.51 0.41 0.36 0.00 0.00 0.00 0.00 0.00  ;
0.55 0.37 0.26 0.61 0.14 0.13 0.23 0.31 0.51 0.41 0.36 0.00 0.00 0.00 0.00  ;
0.47 0.55 0.37 0.26 0.61 0.14 0.13 0.23 0.31 0.51 0.41 0.36 0.00 0.00 0.00  ;
0.36 0.47 0.55 0.37 0.26 0.61 0.14 0.13 0.23 0.31 0.51 0.41 0.36 0.00 0.00  ;
0.54 0.36 0.47 0.55 0.37 0.26 0.61 0.14 0.13 0.23 0.31 0.51 0.41 0.36 0.00  ;
0.57 0.54 0.36 0.47 0.55 0.37 0.26 0.61 0.14 0.13 0.23 0.31 0.51 0.41 0.36  ;
0.34 0.57 0.54 0.36 0.47 0.55 0.37 0.26 0.61 0.14 0.13 0.23 0.31 0.51 0.41  ;
0.23 0.34 0.57 0.54 0.36 0.47 0.55 0.37 0.26 0.61 0.14 0.13 0.23 0.31 0.51  ;
0.14 0.23 0.34 0.57 0.54 0.36 0.47 0.55 0.37 0.26 0.61 0.14 0.13 0.23 0.31  ;
0.47 0.14 0.23 0.34 0.57 0.54 0.36 0.47 0.55 0.37 0.26 0.61 0.14 0.13 0.23  ;
0.34 0.47 0.14 0.23 0.34 0.57 0.54 0.36 0.47 0.55 0.37 0.26 0.61 0.14 0.13  ;
0.00 0.34 0.47 0.14 0.23 0.34 0.57 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  ;
0.00 0.00 0.34 0.47 0.14 0.23 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  ;
0.00 0.00 0.00 0.34 0.47 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  ;
0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  ;
0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  ;
0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  ;
0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  ;

];
0.87 0.33 0.36 0.22 1.32 0.47 0.31 0.25 0.42 0.69 0.51 0.50 0.86 0.59 0.43  ;
*/
