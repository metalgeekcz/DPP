#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
// g++ -o main main.cc -lgsl -lgslcblas -lm

using namespace std;

void eig_decom (double data[], int dim, vector<double> & pEigVals, vector< vector<double> > & pEigVecs);

vector<int> sample_k(const vector<double> &pEigVals, const int K);

vector<int> sample_dpp(const vector<double> &pEigVals, const vector<vector<double> > & pEigVecs, const int K);

