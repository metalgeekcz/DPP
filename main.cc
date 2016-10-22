#include "dpp.h"

int main (void)
{
  //the data is the similarity measure matrix which can be obtained such as B*B' where B is doc*feature matrix
   double data[] ={    
     
    5.2100,    5.3000,    5.2000,   12.2110,   12.2550,   12.3000,    9.5550,    9.5110,    9.5550,
    5.3000,    5.4100,    5.3050,   12.5100,   12.5500,   12.6050,    9.2500,    9.2100,    9.2500,
    5.2000,    5.3050,    5.2025,   12.2600,   12.3000,   12.3525,    9.1500,    9.1100,    9.1500,
   12.2110,   12.5100,   12.2600,   29.0401,   29.1205,   29.2700,   20.1505,   20.0701,   20.1505,
   12.2550,   12.5500,   12.3000,   29.1205,   29.2025,   29.3500,   20.3525,   20.2705,   20.3525,
   12.3000,   12.6050,   12.3525,   29.2700,   29.3500,   29.5025,   20.2000,   20.1200,   20.2000,
    9.5550,    9.2500,    9.1500,   20.1505,   20.3525,   20.2000,   29.5025,   29.3005,   29.5025,
    9.5110,    9.2100,    9.1100,   20.0701,   20.2705,   20.1200,   29.3005,   29.1001,   29.3005,
    9.5550,    9.2500,    9.1500,   20.1505,   20.3525,   20.2000,   29.5025,   29.3005,   29.5025
   };
  int dim=9;
  int K=3;
  DPP myDPP(data, dim); 
  myDPP.eig_decom ();

  //sampling from a DPP
  srand(time(0));
  for(int rr=0; rr<5; rr++){
      vector<int> Sample_points;
      Sample_points=myDPP.sample_dpp(K);

      cout<<"Sample_points \n";
      for (int ss=0; ss<Sample_points.size(); ss++){
          cout<<Sample_points[ss]<<" ";
      }
      cout<<endl;
  }
  return 0;
}

