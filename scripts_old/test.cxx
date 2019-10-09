using namespace std;

#include <TFile.h>
#include "TString.h"
#include "string"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "iostream"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TCutG.h"
#include <vector>
#include <TObject.h>
#include <TClass.h>
#include "TF1.h"
#include "TSpectrum.h"
#include "fstream"
#include "TGraph.h"
#include "TCanvas.h"



void test()
{
  TGraph* gr[128];
TGraph* gr40[128];
TGraph* gr3[128];
 
  
   //open the parameter file
 // Double_t peakX[127][7];
  //double strip[128], peakX1[128],peakX2[128],peakX3[128],peakX4[128],peakX5[128],peakX6[128];
  double strip[128], ped[128], peakX1[128], peakX2[128], peakX3[128];
  
  ifstream Alpha_peaks;
  Alpha_peaks.open("/home/jerome/12Be_exp/scripts/Yd_AlphaPeaks.txt");
  
  
  if(Alpha_peaks.is_open())
  {
    for(int i=0;i<128;i++)
    {
      //  Alpha_peaks >> strip[i] >> ped[i] >> peakX1[i] >> peakX2[i] >> peakX3[i];
      Alpha_peaks >> strip[i] >> peakX1[i] >> peakX2[i] >> peakX3[i];
    }
  }
  
  else{cout << "No file found " << endl;}
 Alpha_peaks.close();
 
 
 double ring[16], comb_loss5[16];
ifstream combloss;
combloss.open("/home/jerome/Software/Energy_Loss/Combined_loss");
if(combloss.is_open())
{for(int i=0; i<16; i++)
                
    {
        combloss >> strip[i] >> comb_loss5[i];
        //cout << comb_loss5[i] << endl;
    }
}
 
//cout << peakX[10][0] << " " << peakX[10][1] << " " << peakX[10][2] << " " << peakX[10][3] << " " << peakX[10][4] << " " << peakX[10][5] << " " << peakX[10][6] << endl;
 //cout <<  strip[10] << " " << peakX1[10] << " " << peakX2[10] << " " << peakX3[10] << " " << peakX4[10] << " " << peakX5[10] << " " << peakX6[10] << endl;
 
    double realEn[3]= {5156,5486,5805};
    double realE[3][128]={0};
    for(int i=0;i<128;i++){
        realE[0][i]=5156;
        realE[1][i]=5486;
        realE[2][i]=5805;}
    
    double diff[16]={0};
    for(int i=0;i<16;i++)//converts MeV/u to KeV for an alpha particle
    {diff[i]=4000*comb_loss5[i];
        //cout << diff[i] << endl;
    }
    
    double Ethe[3][128];
    
   for(int i=0; i<128;i++) 
   { for(int j=0; j<8; j++)
{if(i==j*16){Ethe[0][i]=realE[0][i]-diff[0];}
else if ( i==16*j+1){Ethe[0][i]=realE[0][i]-diff[1];}
else if ( i==16*j+2){Ethe[0][i]=realE[0][i]-diff[2];}
else if ( i==16*j+3){Ethe[0][i]=realE[0][i]-diff[3];}
else if ( i==16*j+4){Ethe[0][i]=realE[0][i]-diff[4];}
else if ( i==16*j+5){Ethe[0][i]=realE[0][i]-diff[5];}
else if ( i==16*j+6){Ethe[0][i]=realE[0][i]-diff[6];}
else if ( i==16*j+7){Ethe[0][i]=realE[0][i]-diff[7];}
else if ( i==16*j+8){Ethe[0][i]=realE[0][i]-diff[8];}
else if ( i==16*j+9){Ethe[0][i]=realE[0][i]-diff[9];}
else if ( i==16*j+10){Ethe[0][i]=realE[0][i]-diff[10];}
else if ( i==16*j+11){Ethe[0][i]=realE[0][i]-diff[11];}
else if ( i==16*j+12){Ethe[0][i]=realE[0][i]-diff[12];}
else if ( i==16*j+13){Ethe[0][i]=realE[0][i]-diff[13];}
else if ( i==16*j+14){Ethe[0][i]=realE[0][i]-diff[14];}
else if ( i==16*j+15){Ethe[0][i]=realE[0][i]-diff[15];}
}
  cout << Ethe[0][i] << endl;     
}

for(int i=0; i<128;i++) 
   { for(int j=0; j<8; j++)
{if(i==j*16){Ethe[1][i]=realE[1][i]-diff[0];}
else if ( i==16*j+1){Ethe[1][i]=realE[1][i]-diff[1];}
else if ( i==16*j+2){Ethe[1][i]=realE[1][i]-diff[2];}
else if ( i==16*j+3){Ethe[1][i]=realE[1][i]-diff[3];}
else if ( i==16*j+4){Ethe[1][i]=realE[1][i]-diff[4];}
else if ( i==16*j+5){Ethe[1][i]=realE[1][i]-diff[5];}
else if ( i==16*j+6){Ethe[1][i]=realE[1][i]-diff[6];}
else if ( i==16*j+7){Ethe[1][i]=realE[1][i]-diff[7];}
else if ( i==16*j+8){Ethe[1][i]=realE[1][i]-diff[8];}
else if ( i==16*j+9){Ethe[1][i]=realE[1][i]-diff[9];}
else if ( i==16*j+10){Ethe[1][i]=realE[1][i]-diff[10];}
else if ( i==16*j+11){Ethe[1][i]=realE[1][i]-diff[11];}
else if ( i==16*j+12){Ethe[1][i]=realE[1][i]-diff[12];}
else if ( i==16*j+13){Ethe[1][i]=realE[1][i]-diff[13];}
else if ( i==16*j+14){Ethe[1][i]=realE[1][i]-diff[14];}
else if ( i==16*j+15){Ethe[1][i]=realE[1][i]-diff[15];}
}}

for(int i=0; i<128;i++) 
   { for(int j=0; j<8; j++)
{if(i==j*16){Ethe[2][i]=realE[2][i]-diff[0];}
else if ( i==16*j+1){Ethe[2][i]=realE[2][i]-diff[1];}
else if ( i==16*j+2){Ethe[2][i]=realE[2][i]-diff[2];}
else if ( i==16*j+3){Ethe[2][i]=realE[2][i]-diff[3];}
else if ( i==16*j+4){Ethe[2][i]=realE[2][i]-diff[4];}
else if ( i==16*j+5){Ethe[2][i]=realE[2][i]-diff[5];}
else if ( i==16*j+6){Ethe[2][i]=realE[2][i]-diff[6];}
else if ( i==16*j+7){Ethe[2][i]=realE[2][i]-diff[7];}
else if ( i==16*j+8){Ethe[2][i]=realE[2][i]-diff[8];}
else if ( i==16*j+9){Ethe[2][i]=realE[2][i]-diff[9];}
else if ( i==16*j+10){Ethe[2][i]=realE[2][i]-diff[10];}
else if ( i==16*j+11){Ethe[2][i]=realE[2][i]-diff[11];}
else if ( i==16*j+12){Ethe[2][i]=realE[2][i]-diff[12];}
else if ( i==16*j+13){Ethe[2][i]=realE[2][i]-diff[13];}
else if ( i==16*j+14){Ethe[2][i]=realE[2][i]-diff[14];}
else if ( i==16*j+15){Ethe[2][i]=realE[2][i]-diff[15];}
}}

double pedYu[128]= {22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,47,22,47,22,22,22,22,22,22,22,22,22,22,22,47,22,22,22,47,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,47,22,47,22,22,22,0,22,22,22,22,22,22,22,22,22,47,47,22,47,0,22,22,22,22,22,22,22,22,22,0,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22};

for(int i=0; i<128;i++)
 {
   
  //double Eexp[3]={peakX2[i],peakX4[i],peakX6[i]};
   double Eexp[3]={peakX1[i],peakX2[i],peakX3[i]}; //3 peak fit
   //double Eexp4[4]={ped2[l][i],peakX1[i],peakX2[i],peakX3[i]}; //3 peak + pedestal
   double Eexp40[4]={ped[i],peakX1[i],peakX2[i],peakX3[i]};
double Eexp3[4]={pedYu[i],peakX1[i],peakX2[i],peakX3[i]};
double Etheo[4]={0,Ethe[0][i],Ethe[1][i],Ethe[2][i]};
gr[i] = new TGraph(3,Eexp,realEn);
  //gr40[i] = new TGraph(4,Eexp40,Etheo);
  gr3[i] = new TGraph(4,Eexp3,Etheo); 
  cout << Etheo[0] << " " << Etheo[1] << " " << Etheo[2] << " " << Etheo[3] << endl;
 }
}
