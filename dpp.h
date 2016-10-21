#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iterator>
#include <algorithm>
#include <memory.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h> //matrix multiplications

// g++ -o main main.cc -lgsl -lgslcblas -lm

using namespace std;

class DPP{
    public:
        DPP(double pData[], int pDim);
        void eig_decom ();
        vector<int> sample_k(const int K);
        vector<int> sample_dpp(const int K);
     private:
     double *mData;
     int mDim; 
     vector<double> mEigVals;
     vector<vector<double> > mEigVecs;
};
