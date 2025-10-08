#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>

#include "fisherinversion.h"
#include "integerPowers.h"

void fisherInversion(vector<double> dF, vector<double> err, int npar, int ndata, string filename, string header, string values, gsl_matrix** covmatrix, gsl_matrix** covinverse)
{

  ofstream out;
  int output=filename.length();
  out << std::scientific;

  if(output)
    {
      out.open(filename.c_str());

      if(!out)
		{
		  cerr << __FILE__ << ": Error: could not open output file (" << filename << ")" << endl;
		  exit(1);
		}
    }
  
  gsl_set_error_handler_off();

  gsl_matrix* bij = gsl_matrix_calloc(npar,npar); //inverse of cov
  gsl_matrix* bij_copy = gsl_matrix_calloc(npar,npar); //copy of the inverse of cov
  gsl_matrix* cij = gsl_matrix_alloc(npar,npar); //cov
  gsl_matrix* identity_test = gsl_matrix_alloc(npar,npar);
  gsl_permutation* permutation = gsl_permutation_alloc(npar);
  int signum;

  //Check allocation proceeded
  if(bij==NULL || cij==NULL || permutation==NULL)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Failed to allocate sufficient matrix memory (npar = " << npar << ")" << endl;
      exit(1);
    }

  out.precision(16);

  //Calculate the bij matrix
  if(output) out << "bij:\n";
  for(int i=0;i<npar;i++)
    {
      int ishift = i*ndata;
      for(int j=0;j<npar;j++)
		{
		  int jshift = j*ndata;

		  double sum=0;
	  
		  for(int idx=0;idx<ndata;idx++)
			{
			  sum += dF[ishift+idx] * dF[jshift+idx] / sqr(err[idx]);
			}
	      
		  gsl_matrix_set(bij,i,j,sum);
		  gsl_matrix_set(bij_copy,i,j,sum);
		  if(output) out << sum << "\t";
	  
		}
      if(output) out << "\n";
    }
  if(output) out << "\n";

  

  //compute the LU decomposition
  gsl_linalg_LU_decomp(bij,permutation,&signum);
  //compute the inverse
  gsl_linalg_LU_invert(bij,permutation,cij);

  //Compute the identity matrix from the original and its inverse
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,bij_copy,cij,0.0,identity_test);


  if(output)
    {
      out << "LU:\n";
      for(int i=0;i<npar;i++)
		{
		  for(int j=0;j<npar;j++)
			{
			  out << gsl_matrix_get(bij,i,j) << "\t";
			}
		  out << "\n";
		}
      out << "\n";

      //print out the matrix
      out << "cij:\n";
      
      for(int i=0;i<npar;i++)
		{
		  for(int j=0;j<npar;j++)
			{
			  out << gsl_matrix_get(cij,i,j) << "\t";
			}
		  out << "\n";
		}
      out << "\n";

	  //print out the matrix
      out << "bij x cij:\n";
      
      for(int i=0;i<npar;i++)
		{
		  for(int j=0;j<npar;j++)
			{
			  out << gsl_matrix_get(identity_test,i,j) << "\t";
			}
		  out << "\n";
		}
      out << "\n";
	  out << "sqrt(diag)-1\n";
	  for(int i=0;i<npar;i++)
		{
		  out << sqrt(gsl_matrix_get(identity_test,i,i))-1 << "\t";
		}
	  out << "\n";
	  out << "max(abs(off-diag))\n";
      for(int i=0;i<npar;i++)
		{
		  double maxel=0;
		  for(int j=0;j<npar;j++)
			{
			  double tmp = gsl_matrix_get(identity_test,i,j);
			  if(i!=j&&abs(tmp)>maxel) maxel=tmp; 
			}
		  out << maxel << " ";
		}
	  out << "\n";


	  

	  
      //Compute the errors
      out << "errors:\n";
      out << header << "\n";
      out << values << "\n";

      for(int i=0;i<npar;i++)
		{
		  out << sqrt(gsl_matrix_get(cij,i,i)) << "\t";
		}
      out << endl;
      out.close();
    }

  if(covmatrix!=NULL)
    {
      //gsl_matrix_free(*covmatrix);
      //*covmatrix = gsl_matrix_alloc(npar,npar);
      gsl_matrix_memcpy(*covmatrix,cij);
    }

  if(covinverse!=NULL)
    {
      //gsl_matrix_free(*covinverse);
      //*covinverse = gsl_matrix_alloc(npar,npar);
      gsl_matrix_memcpy(*covinverse,bij);
    }

  gsl_matrix_free(bij);
  gsl_matrix_free(bij_copy);
  gsl_matrix_free(cij);
  gsl_matrix_free(identity_test);
  gsl_permutation_free(permutation);

}
