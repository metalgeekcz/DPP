#include "dpp.h"

int main (void)
{
 /* double data[] = { -1.0, 1.0, -1.0, 1.0,
                    -8.0, 4.0, -2.0, 1.0,
                    27.0, 9.0, 3.0, 1.0,
                    64.0, 16.0, 4.0, 1.0 };*/
   double data[] ={  5.0 ,    7.0  ,  10.0  ,   5.0  ,   9.0,
                     7.0 ,   26.0  ,  14.0  ,   7.0  ,  27.0,
                    10.0 ,   14.0  ,  20.0  ,  10.0  ,  18.0,
                     5.0 ,    7.0  ,  10.0  ,   5.0  ,   9.0,
                     9.0 ,   27.0  ,  18.0  ,   9.0  ,  29.0};
  int dim=5;
  
  vector<double> lEigVals;
  vector<vector<double> > lEigVecs;
  eig_decom (data, dim, lEigVals, lEigVecs);
  //sampling from a DPP
  int K=3;
  vector<int> Sample_points;
  Sample_points=sample_dpp(lEigVals, lEigVecs, K);

  cout<<"Sample_points \n";
  for (int ss=0; ss<Sample_points.size(); ss++){
      cout<<Sample_points[ss]<<" ";
  }
  
  return 0;
}

