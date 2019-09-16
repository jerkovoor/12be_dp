using namespace std;

#include "TFile.h" 
#include "TMath.h"
#include <cmath>
#include "iostream"
#include "fstream"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"

void BeamOffset(){
    
    double A, B;

	double OffsetInitialx=-0.5;
    double OffsetFinalx=0.5;
    double OffsetInitialy=0.5;
    double OffsetFinaly=-0.5;
    double interval=0.2;
    int OffsetLengthx=1+(OffsetFinalx-OffsetInitialx)/interval;
    int OffsetLengthy=1+fabs(OffsetFinaly-OffsetInitialy)/interval;
    int Length=OffsetLengthx*OffsetLengthy;
    
    double x[OffsetLengthx], y[OffsetLengthy];
    
    for (int i=0;i<OffsetLengthx;i++){
        x[i]=OffsetInitialx+i*interval;
    }
    
    for (int i=0;i<OffsetLengthy;i++){
        y[i]=OffsetInitialy+i*interval;
    }
    
    double phi[8], r[16], R[8][16][OffsetLengthx][OffsetLengthy], DetAngle[8][16][OffsetLengthx][OffsetLengthy];
    phi[0]=22.5;
    for (int i=1;i<8;i++){
        phi[i]=phi[0]+45*i;
        //cout << phi[i] << endl;
    }
    
    for (int i=0;i<16;i++){
        r[i]=50+((129.-50.)/16)*(0.5+i);
        //cout << r[i] << endl;
    }
    
    //Offsetting the beam
        for (int xindex=0;xindex<OffsetLengthx;xindex++){
            for (int yindex=0;yindex<OffsetLengthy;yindex++){
                //Calculating the angles
                for (int i=0;i<8;i++){
                    for (int j=0;j<16;j++){
                        A=r[j]*TMath::Sin(phi[i]);
                        B=r[j]*TMath::Cos(phi[i]);
                        R[i][j][xindex][yindex]=TMath::Sqrt(pow((x[xindex]-A),2)+pow((y[yindex]-B),2));
                        DetAngle[i][j][xindex][yindex]=180-TMath::RadToDeg()*TMath::ATan(R[i][j][xindex][yindex]/85);
                        cout << DetAngle[i][j][xindex][yindex] << endl;
                    }
                }
            }
        }

}
