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



double calAl()
{
  TGraph* gr[128];
TGraph* gr40[128];
TGraph* gr3[128];
 
  
   //open the parameter file
 // Double_t peakX[127][7];
  //double strip[128], peakX1[128],peakX2[128],peakX3[128],peakX4[128],peakX5[128],peakX6[128];
  double strip[128], ped[128], peakX1[128], peakX2[128], peakX3[128];
  
  ifstream Alpha_peaks;
  Alpha_peaks.open("/home/jerome/12Be_exp/scripts/Yu_AlphaPeaks5225.txt");
  
  
  if(Alpha_peaks.is_open())
  {
    for(int i=0;i<128;i++)
    {
      Alpha_peaks >> strip[i] >> ped[i] >> peakX1[i] >> peakX2[i] >> peakX3[i];
      //Alpha_peaks >> strip[i] >> peakX1[i] >> peakX2[i] >> peakX3[i];
    }
  }
  
  else{cout << "No file found " << endl;}
 Alpha_peaks.close();
 
 
 double alpha_loss1[16];
ifstream loss1;
loss1.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Alpha_calibration/alpha1.txt");
if(loss1.is_open())
{for(int i=0; i<16; i++)
                
    {
        loss1 >> alpha_loss1[i];
    }
}

double alpha_loss2[16];
ifstream loss2;
loss2.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Alpha_calibration/alpha2.txt");
if(loss2.is_open())
{for(int i=0; i<16; i++)
                
    {
        loss2 >> alpha_loss2[i];
    }
}

double alpha_loss3[16];
ifstream loss3;
loss3.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Alpha_calibration/alpha3.txt");
if(loss3.is_open())
{for(int i=0; i<16; i++)
                
    {
        loss3 >> alpha_loss3[i];
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
    {diff[i]=4000*alpha_loss1[i];
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
  //cout << Ethe[0][i] << endl;     
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

 
/*else if ( i==16*j+1){Ethe[2][i]=realE[2][i]-diff[1];}
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
	}
}}*/
    
    //double Etheo[4]= {0,5156,5486,5805};
    double pedYu[128]= {22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,47,22,47,22,22,22,22,22,22,22,22,22,22,22,47,22,22,22,47,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,47,22,47,22,22,22,0,22,22,22,22,22,22,22,22,22,47,47,22,47,0,22,22,22,22,22,22,22,22,22,0,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22};
    
    double pedYd[128]={0,34,0,34,34,34,34,38,38,34,34,34,38,34,34,34,36,34,34,34,0,34,34,34,34,34,38,34,34,36,34,38,34,34,34,34,34,34,34,34,34,34,38,38,38,34,34,34,34,34,34,34,34,36,38,34,34,38,38,38,38,38,38,34,0,34,34,34,34,34,34,34,34,34,38,38,38,38,34,32,34,34,34,34,34,38,34,38,34,34,34,38,38,38,38,0,34,34,34,34,34,34,34,34,34,36,38,34,34,34,38,34,38,34,0,34,38,34,34,38,34,34,38,38,38,34,34,0};
    
    double ped2[5][128];
    
  for(int k=0; k<128; k++){ for(int p=0; p<5; p++){ped2[p][k]={20.+p};}}
  
 
 for(int i=0; i<128;i++)
 {
   
  //double Eexp[3]={peakX2[i],peakX4[i],peakX6[i]};
   double Eexp[3]={peakX1[i],peakX2[i],peakX3[i]}; //3 peak fit
   //double Eexp4[4]={ped2[l][i],peakX1[i],peakX2[i],peakX3[i]}; //3 peak + pedestal
   double Eexp40[4]={ped[i],peakX1[i],peakX2[i],peakX3[i]};
double Eexp3[4]={pedYu[i],peakX1[i],peakX2[i],peakX3[i]};
double Etheo[4]={0,Ethe[0][i],Ethe[1][i],Ethe[2][i]};

  //gr[i] = new TGraph(3,Eexp,realE);
  gr3[i] = new TGraph(4,Eexp40,Etheo);
  //gr3[i] = new TGraph(4,Eexp3,Etheo);
 }
 
  //open a txt file to store the positions of the peaks
  ofstream Yd_alpha3;
  Yd_alpha3.open("Yu5225_pedestal_loss_AllPar.txt");
  Yd_alpha3 << " strip / gain / offset  \n";
  
  /*ofstream Yd_alpha4;
  Yd_alpha4.open("Yu5226_AllPar.txt");
  Yd_alpha4 << " strip / gain / offset  \n";*/
 
 
    TF1 *fit_lin = new TF1 ( "fit_lin","[0]*x+[1]",0,4000 );
    TF1 *fit_lin1 = new TF1 ( "fit_lin1","[0]*x+[1]",0,4000 );
   //  TF1 *fit_lin1 = new TF1 ( "fit_lin1","[0]*x*x+[1]*x+[2]",0,4000 );
    
  /*  for(int m=0; m<5; m++)
    {
       gr4[0][m]->Draw ( "AP*" );
       gr4[0][m]->Fit ( "fit_lin1","","",0,4000 );
    }
    
     gr40[0]->Draw ( "AP*" );
     gr40[0]->Fit ( "fit_lin1","","",0,4000 );*/

    TCanvas *c1 = new TCanvas ( "c1" );
    c1->Divide ( 4,4 );
    
    
    for ( int i=0; i<16; i++ )
    {
        c1->cd ( i+1 );
       /*gr[i]->Draw ( "AP*" );
       gr[i]->Fit ( "fit_lin1","","",0,4000 );
       Yd_alpha4 << i << " "  << fit_lin1->GetParameter ( 0 ) << " " << fit_lin1->GetParameter ( 1 ) << endl;*/
	//Yd_alpha4 << i << " "  << fit_lin1->GetParameter ( 0 ) << " " << fit_lin1->GetParameter ( 1 ) << " " << fit_lin1->GetParameter ( 2 ) << endl;

         gr3[i]->Draw ( "AP*" );
         gr3[i]->Fit ( "fit_lin","","",0,4000 );
         Yd_alpha3 << i << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;
    
    }

    TCanvas *c2 = new TCanvas ( "c2" );
    c2->Divide ( 4,4 );
    for ( int i=0; i<16; i++ )
    {
        c2->cd ( i+1 );
       /*gr[i+16]->Draw ( "AP*" );
       gr[i+16]->Fit ( "fit_lin1","","",0,4000 );
       Yd_alpha4 << i+16 << " "  << fit_lin1->GetParameter ( 0 ) << " " << fit_lin1->GetParameter ( 1 ) << endl;*/
	//Yd_alpha4 << i+16 << " "  << fit_lin1->GetParameter ( 0 ) << " " << fit_lin1->GetParameter ( 1 ) << " " << fit_lin1->GetParameter ( 2 ) << endl;
       
         gr3[i+16]->Draw ( "AP*" );
         gr3[i+16]->Fit ( "fit_lin","","",0,4000 );
         Yd_alpha3 << i+16 << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;  
    }

    TCanvas *c3 = new TCanvas ( "c3" );
    c3->Divide ( 4,4 );
    for ( int i=0; i<16; i++ )
    {
        c3->cd ( i+1 );
      /* gr[i+32]->Draw ( "AP*" );
       gr[i+32]->Fit ( "fit_lin1","","",0,4000 );
       Yd_alpha4 << i+32 << " "  << fit_lin1->GetParameter ( 0 ) << " " << fit_lin1->GetParameter ( 1 ) << endl;*/
	//Yd_alpha4 << i+32 << " "  << fit_lin1->GetParameter ( 0 ) << " " << fit_lin1->GetParameter ( 1 ) << " " << fit_lin1->GetParameter ( 2 ) << endl;

         gr3[i+32]->Draw ( "AP*" );
         gr3[i+32]->Fit ( "fit_lin","","",0,4000 );
         Yd_alpha3 << i+32 << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;  
    }

    TCanvas *c4 = new TCanvas ( "c4" );
    c4->Divide ( 4,4 );
    for ( int i=0; i<16; i++ )
    {
        c4->cd ( i+1 );
      /*gr[i+48]->Draw ( "AP*" );
       gr[i+48]->Fit ( "fit_lin1","","",0,4000 );
      Yd_alpha4 << i+48 << " "  << fit_lin1->GetParameter ( 0 ) << " " << fit_lin1->GetParameter ( 1 ) << endl;*/
	//Yd_alpha4 << i+48 << " "  << fit_lin1->GetParameter ( 0 ) << " " << fit_lin1->GetParameter ( 1 ) << " " << fit_lin1->GetParameter ( 2 ) << endl;

        gr3[i+48]->Draw ( "AP*" );
        gr3[i+48]->Fit ( "fit_lin","","",0,4000 );
        Yd_alpha3 << i+48 << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;  
    }

    TCanvas *c5 = new TCanvas ( "c5" );
    c5->Divide ( 4,4 );
    for ( int i=0; i<16; i++ )
    {
        c5->cd ( i+1 );
      /* gr[i+64]->Draw ( "AP*" );
      gr[i+64]->Fit ( "fit_lin1","","",0,4000 );
    Yd_alpha4 << i+64 << " "  << fit_lin1->GetParameter ( 0 ) << " " << fit_lin1->GetParameter ( 1 ) << endl;*/
	//Yd_alpha4 << i+64 << " "  << fit_lin1->GetParameter ( 0 ) << " " << fit_lin1->GetParameter ( 1 ) << " " << fit_lin1->GetParameter ( 2 ) << endl;

         gr3[i+64]->Draw ( "AP*" );
        gr3[i+64]->Fit ( "fit_lin","","",0,4000 );
        Yd_alpha3 << i+64 << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;  
    }

    TCanvas *c6 = new TCanvas ( "c6" );
    c6->Divide ( 4,4 );
    for ( int i=0; i<16; i++ )
    {
        c6->cd ( i+1 );
       /*gr[i+80]->Draw("AP*" );
       gr[i+80]->Fit( "fit_lin1","","",0,4000 );
      Yd_alpha4 << i+80 << " "  << fit_lin1->GetParameter ( 0 ) << " " << fit_lin1->GetParameter ( 1 ) << endl;*/
	//Yd_alpha4 << i+80 << " "  << fit_lin1->GetParameter ( 0 ) << " " << fit_lin1->GetParameter ( 1 ) << " " << fit_lin1->GetParameter ( 2 ) << endl;

         gr3[i+80]->Draw ( "AP*" );
         gr3[i+80]->Fit ( "fit_lin","","",0,4000 );
         Yd_alpha3 << i+80 << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;  
    }

    TCanvas *c7 = new TCanvas ( "c7" );
    c7->Divide ( 4,4 );
    for ( int i=0; i<16; i++ )
    {
        c7->cd ( i+1 );
       /*gr[i+96]->Draw( "AP*" );
       gr[i+96]->Fit( "fit_lin1","","",0,4000 );
       Yd_alpha4 << i+96 << " "  << fit_lin1->GetParameter ( 0 ) << " " << fit_lin1->GetParameter ( 1 ) << endl;*/
	//Yd_alpha4 << i+96 << " "  << fit_lin1->GetParameter ( 0 ) << " " << fit_lin1->GetParameter ( 1 ) << " " << fit_lin1->GetParameter ( 2 ) << endl;

         gr3[i+96]->Draw ( "AP*" );
         gr3[i+96]->Fit ( "fit_lin","","",0,4000 );
         Yd_alpha3 << i+96 << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;  
    }

    TCanvas *c8 = new TCanvas ( "c8" );
    c8->Divide ( 4,4 );
    for ( int i=0; i<16; i++ )
    {
        c8->cd ( i+1 );
       /*gr[i+112]->Draw( "AP*" );
       gr[i+112]->Fit( "fit_lin1","","",0,4000 );
       Yd_alpha4 << i+112 << " "  << fit_lin1->GetParameter ( 0 ) << " " << fit_lin1->GetParameter ( 1 ) << endl;*/
	//Yd_alpha4 << i+112 << " "  << fit_lin1->GetParameter ( 0 ) << " " << fit_lin1->GetParameter ( 1 ) << " " << fit_lin1->GetParameter ( 2 ) << endl;

         gr3[i+112]->Draw ( "AP*" );
         gr3[i+112]->Fit ( "fit_lin","","",0,4000 );
         Yd_alpha3 << i+112 << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;  
    }
    
    
    TFile *f_out = new TFile("AlphaYu5225_pedestal_loss.root","RECREATE");
    
    f_out->cd();
    
    for ( int i=0; i<128; i++ )
   {
         //  string name1 = Form ( "4points_%i",i );
	string name2 = Form ( "3points_%i",i );
       // gr[i]->Write(name1.c_str());
        gr3[i]->Write(name2.c_str());
      
    }
      
    f_out->Close();
 
return 0;
} //end of script
