#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

using namespace std;
void eig_decom (double data[], int dim, vector<double> & pEigVals, vector< vector<double> > & pEigVecs)
{
  
  //compute eig decomposition 
  gsl_matrix_view m 
    = gsl_matrix_view_array (data, dim, dim);

  gsl_vector_complex *eval = gsl_vector_complex_alloc (dim);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (dim, dim);

  gsl_eigen_nonsymmv_workspace * w = 
    gsl_eigen_nonsymmv_alloc (dim);
  
  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);

  gsl_eigen_nonsymmv_free (w);

  gsl_eigen_nonsymmv_sort (eval, evec, 
                           GSL_EIGEN_SORT_ABS_DESC);
  
  //store the real part in vectors


  {
    int i, j;

    for (i = 0; i < dim; i++)
      {
        gsl_complex eval_i 
           = gsl_vector_complex_get (eval, i);
        gsl_vector_complex_view evec_i 
           = gsl_matrix_complex_column (evec, i);
        
        pEigVals.push_back(GSL_REAL(eval_i));
        
        vector<double> lOneEigVec;
        for (j = 0; j < dim; ++j)
          {
            gsl_complex z = 
              gsl_vector_complex_get(&evec_i.vector, j);

            lOneEigVec.push_back(GSL_REAL(z));
          }
          pEigVecs.push_back(lOneEigVec);
      }
  }

  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);

}

// the returned samples index from 1
vector<int> sample_k(const vector<double> &pEigVals, const int K){

    // compute the elementary symmetric polyonomials 
    vector<vector<double> > pE; //

    int dim=pEigVals.size();

    vector<double> firstrow(dim+1, 1.0);
    pE.push_back(firstrow);
    for(int i=0; i<K; i++){
        std::vector<double> onerow(dim+1, 0.0);
        for(int n=0; n<dim; n++){
           onerow[n+1] = onerow[n] + pEigVals[n]*pE[i][n];
        }
        pE.push_back(onerow);
    }

    //iterate 
    int remaining = K;
    int i=dim;
    double marg=0.0;

    vector<int> Samples(K,0);
    while (remaining){
        if( remaining == i)
            marg=1.0;
        else
            marg = pEigVals[i-1] * pE[remaining-1][i-1]/pE[remaining][i];
        
        srand(time(0));
        double randN= (double)rand()/RAND_MAX;
        if(randN < marg) {
            Samples[remaining-1]=i;
            remaining--;
        }
        i--; 
    }
    return Samples;
}

