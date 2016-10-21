#include "dpp.h"

DPP::DPP( double pData[], int pDim)
: mDim(pDim){
    mData= new double[pDim*pDim];
   memcpy(mData, pData, pDim*pDim*sizeof(double));
}

void DPP::eig_decom ()
{
  
  //compute eig decomposition 
  gsl_matrix_view m 
    = gsl_matrix_view_array (mData, mDim, mDim);

  gsl_vector_complex *eval = gsl_vector_complex_alloc (mDim);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (mDim, mDim);

  gsl_eigen_nonsymmv_workspace * w = 
    gsl_eigen_nonsymmv_alloc (mDim);
  
  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);

  gsl_eigen_nonsymmv_free (w);

  gsl_eigen_nonsymmv_sort (eval, evec, 
                           GSL_EIGEN_SORT_ABS_DESC);
  
  //store the real part in vectors


  {
    int i, j;

    for (i = 0; i < mDim; i++)
      {
        gsl_complex eval_i 
           = gsl_vector_complex_get (eval, i);
        gsl_vector_complex_view evec_i 
           = gsl_matrix_complex_column (evec, i);
        
        mEigVals.push_back(GSL_REAL(eval_i));
        
        vector<double> lOneEigVec;
        for (j = 0; j < mDim; ++j)
          {
            gsl_complex z = 
              gsl_vector_complex_get(&evec_i.vector, j);

            lOneEigVec.push_back(GSL_REAL(z));
          }
         mEigVecs.push_back(lOneEigVec);
      }
  }

  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);

}

// the returned samples index from 1
vector<int> DPP::sample_k( const int K){

    // compute the elementary symmetric polyonomials 
    vector<vector<double> > pE; //
    vector<double> firstrow(mDim+1, 1.0);
    pE.push_back(firstrow);
    for(int i=0; i<K; i++){
        std::vector<double> onerow(mDim+1, 0.0);
        for(int n=0; n<mDim; n++){
           onerow[n+1] = onerow[n] + mEigVals[n]*pE[i][n];
        }
        pE.push_back(onerow);
    }

    //iterate 
    int remaining = K;
    int i=mDim;
    double marg=0.0;

    vector<int> Samples(K,0);
    while (remaining){
        if( remaining == i)
            marg=1.0;
        else
            marg = mEigVals[i-1] * pE[remaining-1][i-1]/pE[remaining][i];
        
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

vector<int> DPP::sample_dpp(const int K){
  
  vector<int> Y(K, 0);
  //Step 1 select a subseti (k) of eignvectors, where the probability of selecting each eignvector depends on its associated eignvalue
  vector<int> Samples=sample_k(K);
  
  vector<vector<double> > V;
  for(int k=0; k<K; k++){
     V.push_back(mEigVecs[Samples[k]-1]);
  }

  //step 2: produce samples based on the selected eigen vectors

  for (int i=K-1; i>=0; i--){
    vector<double> P(mDim,0.0);
    double sumP=0.0;
    for(int dd=0; dd< mDim; dd++){
        for(int k=0; k<V.size(); k++){
            P[dd]+=pow(V[k][dd] , 2.0);
        }
        sumP+=P[dd];
    }
    
    vector<double> CumSumP(mDim,0.0);
    for(int dd=0; dd< mDim; dd++){
        P[dd]/=sumP;
        if(dd==0)
            CumSumP[dd]=P[dd];
        else
            CumSumP[dd]=CumSumP[dd-1]+P[dd];
    }

    srand(time(0));
    double randI= (double) rand()/RAND_MAX;
    for(int dd=0; dd< mDim; dd++){
        if(randI<CumSumP[dd]){
            Y[i]=dd;
            break;
        }
    }

    
    //eliminate the last vector and update V
     vector<vector<double> > Vnew;
     int currentR=V.size();
     for(int k=0; k<currentR-1; k++){
         vector<double> Vnew_one(V[0].size(),0.0);
         for(int d=0; d<V[0].size(); d++){
            Vnew_one[d]=V[k][d]-V.back().at(d)*V[k][Y[i]]/V.back().at(Y[i]);
         }
         Vnew.push_back(Vnew_one);

     }
    
     V=Vnew;

     //orthogonalize
     for(int a=0; a< i-1; a++){
         for(int  b= 0; b<a-1; b++){
             double Stmp=0.0; //V(:,a)'*V(:,b)
             for(int dd=0; dd<mDim; dd++){
                 Stmp+=V[a][dd]*V[b][dd];
             }
             for(int dd=0; dd<mDim; dd++){
                 V[a][dd]-=Stmp*V[b][dd];
             }

         }
         double S_norm=0.0;
         for(int dd=0; dd<mDim; dd++){
             S_norm+=pow(V[a][dd],2.0);
         }
         S_norm=sqrt(S_norm);
         for(int dd=0; dd<mDim; dd++){
             V[a][dd]/=S_norm;
         }
     }
  }

   return Y;
}
