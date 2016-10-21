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
        void eig_decom (double data[], int dim, vector<double> & pEigVals, vector< vector<double> > & pEigVecs);

        vector<int> sample_k(const vector<double> &pEigVals, const int K);

        vector<int> sample_dpp(const vector<double> &pEigVals, const vector<vector<double> > & pEigVecs, const int K);
    
};
