#include "dpp.h"

vector<int> sample_dpp(const vector<double> &pEigVals, const vector<vector<double> > & pEigVecs, const int K){
  
  vector<int> Y(K, 0);
  int Dim=pEigVals.size();
  //Step 1 select a subseti (k) of eignvectors, where the probability of selecting each eignvector depends on its associated eignvalue
  vector<int> Samples=sample_k(pEigVals, K);
  
  vector<vector<double> > V;
  for(int k=0; k<K; k++){
     V.push_back(pEigVecs[Samples[k]-1]);
  }

  //step 2: produce samples based on the selected eigen vectors

  for (int i=K-1; i>=0; i--){
    vector<double> P(Dim,0.0);
    double sumP=0.0;
    for(int dd=0; dd< Dim; dd++){
        for(int k=0; k<V.size(); k++){
            P[dd]+=pow(V[k][dd] , 2.0);
        }
        sumP+=P[dd];
    }
    
    vector<double> CumSumP(Dim,0.0);
    for(int dd=0; dd< Dim; dd++){
        P[dd]/=sumP;
        if(dd==0)
            CumSumP[dd]=P[dd];
        else
            CumSumP[dd]=CumSumP[dd-1]+P[dd];
    }

    srand(time(0));
    double randI= (double) rand()/RAND_MAX;
    for(int dd=0; dd< Dim; dd++){
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
             for(int dd=0; dd<Dim; dd++){
                 Stmp+=V[a][dd]*V[b][dd];
             }
             for(int dd=0; dd<Dim; dd++){
                 V[a][dd]-=Stmp*V[b][dd];
             }

         }
         double S_norm=0.0;
         for(int dd=0; dd<Dim; dd++){
             S_norm+=pow(V[a][dd],2.0);
         }
         S_norm=sqrt(S_norm);
         for(int dd=0; dd<Dim; dd++){
             V[a][dd]/=S_norm;
         }
     }
  }

   return Y;
}
