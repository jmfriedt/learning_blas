#include <stdio.h>
#include <math.h>
#include <complex.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>

int main()
{
  int nobs=2000;
  int nlag=20;
  int i;
  int l,m;
  gsl_matrix_complex *host_mem;
  gsl_vector_complex_view tmp_vector_view1,tmp_vector_view2;
  gsl_vector_complex *host_val,*host_code,*host_res;
  gsl_complex pwr;
  gsl_complex z;
  gsl_complex alpha=1+0*I;
  gsl_complex beta=0.+0*I;

// 1st demo: power as dot product
  host_val=gsl_vector_complex_alloc(nobs);
  host_mem=gsl_matrix_complex_alloc(nobs,nlag);
  
  for (i=0;i<nobs;i++)
    {GSL_REAL(z)=1.; GSL_IMAG(z)=1.;
     gsl_vector_complex_set(host_val, i, z); // gsl_complex_conjugate(z));
    }
  gsl_blas_zdotc(host_val, host_val, &pwr); //cublasZdotc(handle, nobs, dev_mem, 1, dev_mem, 1, &pwr);
#ifdef debug
  gsl_vector_complex_fprintf(stdout, host_val, "%g");
#endif
  printf("%% power %lf vs a=ones(%d,1)+j*ones(%d,1);a'*a\n",gsl_complex_abs(pwr),nobs,nobs);
  gsl_vector_complex_free(host_val);
  gsl_matrix_complex_free(host_mem);

// 2nd demo: xcorr as matrix multiplication for detecting time delayed weighted copies of a known code
//https://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
  host_val=gsl_vector_complex_alloc(nobs);
  host_code=gsl_vector_complex_alloc(nobs);
  host_res=gsl_vector_complex_alloc(nlag*2+1);
  host_mem=gsl_matrix_complex_alloc(nobs,nlag*2+1);
  for (m=0;m<nobs;m++) 
      {GSL_REAL(z)=(double)(random()/pow(2,31))-0.5; 
       GSL_IMAG(z)=(double)(random()/pow(2,31))-0.5; ;
       gsl_vector_complex_set(host_val,m,z);
       GSL_REAL(z)=(double)(random()/pow(2,31))-0.5; 
       GSL_IMAG(z)=(double)(random()/pow(2,31))-0.5; ;
       gsl_vector_complex_set(host_code,m,z);
      }
//  gsl_vector_complex_add(host_val,host_code);
//  gsl_vector_complex_fprintf(stdout, host_val, "%g");
  tmp_vector_view1=gsl_vector_complex_subvector(host_val,0,nobs-12);
  tmp_vector_view2=gsl_vector_complex_subvector(host_code,12,nobs-12);
  gsl_vector_complex_scale(&tmp_vector_view2.vector, 0.3);
  gsl_vector_complex_add(&tmp_vector_view1.vector,&tmp_vector_view2.vector);
  tmp_vector_view1=gsl_vector_complex_subvector(host_val,0,nobs-3);
  tmp_vector_view2=gsl_vector_complex_subvector(host_code,3,nobs-3);
  gsl_vector_complex_scale(&tmp_vector_view2.vector, 1.3);
  gsl_vector_complex_add(&tmp_vector_view1.vector,&tmp_vector_view2.vector);
  tmp_vector_view1=gsl_vector_complex_subvector(host_val,10,nobs-10);
  tmp_vector_view2=gsl_vector_complex_subvector(host_code,0,nobs-10);
  gsl_vector_complex_scale(&tmp_vector_view2.vector, 0.8);
  gsl_vector_complex_add(&tmp_vector_view1.vector,&tmp_vector_view2.vector);
  tmp_vector_view1=gsl_vector_complex_subvector(host_val,5,nobs-5);
  tmp_vector_view2=gsl_vector_complex_subvector(host_code,0,nobs-5);
  gsl_vector_complex_add(&tmp_vector_view1.vector,&tmp_vector_view2.vector);
//  gsl_vector_complex_fprintf(stdout, host_val, "%g");
/*
  for (m=0;m<nobs-12;m++)          // time shifted copies of the code
      {host_val[m+12]+=host_code[m];
       host_val[m+3]+=host_code[m];
       host_val[m]+=host_code[m+12];
       host_val[m]+=host_code[m+3];
      }
*/
  gsl_matrix_complex_set_zero(host_mem);  // memset(host_mem,0x0,sizeof(std::complex<double>)*nobs*nlag);
  for (l=-nlag;l<=nlag;l++)
      {                      //    for (m=0;m<nobs-l;m++)
#ifdef debug
       printf("lag=%d\n",l);
#endif
       if (l<0)
         {tmp_vector_view1=gsl_matrix_complex_subcolumn(host_mem, l+nlag,0,nobs+l);
          tmp_vector_view2=gsl_vector_complex_subvector(host_code,-l,nobs+l);
          gsl_vector_complex_memcpy( &tmp_vector_view1.vector, &tmp_vector_view2.vector);       
         }
       else 
         {tmp_vector_view1=gsl_matrix_complex_subcolumn(host_mem, l+nlag,l,nobs-(l));
          tmp_vector_view2=gsl_vector_complex_subvector(host_code,0,nobs-(l));
          gsl_vector_complex_memcpy( &tmp_vector_view1.vector, &tmp_vector_view2.vector);       
         }
#ifdef debug
       tmp_vector_view1=gsl_matrix_complex_subcolumn(host_mem, l+nlag,0,nobs);
       gsl_vector_complex_fprintf(stdout, &tmp_vector_view1.vector, "%g");
#endif
//       host_mem[(l+nlag)*nobs+m+l+nlag]=host_code[m];
      }
//  gsl_matrix_complex_fprintf(stdout, host_mem, "%g");
//  printf("done\n");
  gsl_blas_zgemv(CblasConjTrans, alpha, host_mem, host_val,beta, host_res); // /!\ ConjTrans
  gsl_vector_complex_fprintf(stdout, host_res, "%g");
  return 0;
}
/*
https://stackoverflow.com/questions/46526274/assigning-complex-values-in-gsl
*/
