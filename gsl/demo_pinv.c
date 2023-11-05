// http://theochem.mercer.edu/pipermail/csc335/2013-November/000101.html
#include <stdio.h>
#include <math.h>
#include <complex.h>  // I definition

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
/* Solving complex linear system
  *      (3+i)x -    2y  = 3-4i
  *      -3x    + (1-2i) = -1+.5i
  */

#undef debug

void add_with_offset(gsl_vector_complex *code, gsl_vector_complex *inout, float scale, int total_len, int offset)
{ gsl_vector_complex_view tmp_vector_view1,tmp_vector_view2;
  gsl_vector_complex *tmp=gsl_vector_complex_alloc(total_len-abs(offset));
  if (offset>=0)
    {tmp_vector_view1=gsl_vector_complex_subvector(inout,0,total_len-abs(offset));
     tmp_vector_view2=gsl_vector_complex_subvector(code,offset,total_len-abs(offset));
    }
  else
    {tmp_vector_view1=gsl_vector_complex_subvector(inout,-offset,total_len-abs(offset));
     tmp_vector_view2=gsl_vector_complex_subvector(code,0,total_len-abs(offset));
    }
  gsl_vector_complex_memcpy( tmp, &tmp_vector_view2.vector);       
  gsl_vector_complex_memcpy( tmp, &tmp_vector_view2.vector);       
  gsl_vector_complex_scale(tmp, scale);
  gsl_vector_complex_add(&tmp_vector_view1.vector,tmp);
  gsl_vector_complex_free(tmp);
}

int main (void)
{double a_data[] = { 3,1, -2,0,
                    -3,0,  1,-2 };
 double b_data[] = { 3,-4,
                    -1,0.5 };
 int s;
 gsl_matrix_complex_view a = gsl_matrix_complex_view_array(a_data,2,2);
 gsl_vector_complex_view b = gsl_vector_complex_view_array(b_data, 2);
 gsl_vector_complex *x = gsl_vector_complex_alloc (2);
 gsl_permutation * p = gsl_permutation_alloc (2);
 gsl_linalg_complex_LU_decomp(&a.matrix, p, &s);
 gsl_linalg_complex_LU_solve (&a.matrix, p, &b.vector, x);
#ifdef debug
 gsl_vector_complex_fprintf(stdout, x, "%g");
#endif
 gsl_permutation_free (p);

 int nobs=5000;
 int nlag=25;
 int l,m;
 gsl_vector_complex *host_val=gsl_vector_complex_alloc(nobs);
//https://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
 gsl_vector_complex *host_code=gsl_vector_complex_alloc(nobs);
 gsl_matrix_complex *host_res=gsl_matrix_complex_alloc((nlag*2+1),(nlag*2+1));
 gsl_matrix_complex *host_mem=gsl_matrix_complex_alloc(nobs,(2*nlag+1));
 gsl_matrix_complex *host_out=gsl_matrix_complex_alloc(nobs,2*nlag+1);
 gsl_vector_complex *host_final=gsl_vector_complex_alloc(2*nlag+1);

// identique a` l'exemple matrix -- see https://www.gnu.org/software/gsl/doc/html/blas.html for BLAS functions in GSL
 gsl_vector_complex_view tmp_vector_view1,tmp_vector_view2;
 gsl_complex z;
 gsl_complex alpha=1+0*I;
 gsl_complex beta=0.+0*I;
 for (m=0;m<nobs;m++) 
     {GSL_REAL(z)=(double)(random()/pow(2,31))-0.5; 
      GSL_IMAG(z)=(double)(random()/pow(2,31))-0.5; ;
      gsl_vector_complex_set(host_val,m,z);
      GSL_REAL(z)=(double)(random()/pow(2,31))-0.5; 
      GSL_IMAG(z)=(double)(random()/pow(2,31))-0.5; ;
      gsl_vector_complex_set(host_code,m,z);
     }
  printf("\n");
  add_with_offset(host_code, host_val, 0.3, nobs, -12);
  add_with_offset(host_code, host_val, 1.3, nobs, -3);
  add_with_offset(host_code, host_val, 0.8, nobs, 10);
  add_with_offset(host_code, host_val, 1.0, nobs, 5);
/*
  gsl_vector_complex *tmp=gsl_vector_complex_alloc(nobs-12);
  tmp_vector_view1=gsl_vector_complex_subvector(host_val,0,nobs-12);
  tmp_vector_view2=gsl_vector_complex_subvector(host_code,12,nobs-12);
  gsl_vector_complex_memcpy( tmp, &tmp_vector_view2.vector);       
  gsl_vector_complex_scale(tmp, 0.3);
  gsl_vector_complex_add(&tmp_vector_view1.vector,tmp);
  gsl_vector_complex_free(tmp);

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
*/
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
// fin identique
  gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, alpha, host_mem, host_mem, beta, host_res); // /!\ ConjTrans
  p = gsl_permutation_alloc (2*nlag+1);
  gsl_linalg_complex_LU_decomp(host_res, p, &s);
  gsl_linalg_complex_LU_invert (host_res, p, host_res);
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, host_mem, host_res, beta, host_out);
  gsl_blas_zgemv(CblasConjTrans, alpha, host_out, host_val, beta, host_final);
  gsl_vector_complex_fprintf(stdout, host_final, "%g");
  printf("%% load t; plot([-%d:%d],abs(t(:,1)+j*t(:,2)))\n",nlag,nlag);
  return 0;
}
