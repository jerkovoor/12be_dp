//Alpha calibration of Sd1 detector

using namespace std;

#include <TFile.h>
#include "TString.h"
#include "string"
#include "TTree.h"
#include "TBranch.h"
#include "iostream"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include <vector>
#include <TObject.h>
#include <TClass.h>
#include "TSpectrum.h"
#include <algorithm>
#include "TF1.h"
#include "fstream"

double Sd_alpha_peak_read()
{
    
    TFile *f = new TFile ( "/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1_AlphaPeaks.root","READ" );
    
    TH1D* hSd1s[32];
    TH1D* hSd1r[24];
    
    for (int i=0;i<24;i++){
        hSd1r[i] = (TH1D*)f->Get(Form("Sd1r_Ch%i",i));
    }
    
    for (int i=0;i<32;i++){
        hSd1s[i] = (TH1D*)f->Get(Form("Sd1s_Ch%i",i));
    }
    
    
        
    TF1* fit_func1 = new TF1 ( "fit_func1","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",500,5000 ); //removed the liner part for comparison
    TF1* fit_func3 = new TF1 ( "fit_func3","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))",0,5000 );
    TF1* fit_func9 = new TF1 ( "fit_func9","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp(-pow((x-[10]),2)/(2*[2]*[2]))+[11]*TMath::Exp(-pow((x-[12]),2)/(2*[2]*[2]))+[13]*TMath::Exp(-pow((x-[14]),2)/(2*[2]*[2]))+[15]*TMath::Exp(-pow((x-[16]),2)/(2*[2]*[2]))+[17]*TMath::Exp(-pow((x-[18]),2)/(2*[2]*[2]))+[19]+[20]*x",2000,5000 );
    
    //Search for peaks in Sd strips
    int npeaks=3;
    TSpectrum* s = new TSpectrum ( 2*npeaks );
    int nfound=0;
    double* txp = s->GetPositionX();
    
    //find peaks and fir the Sd1r detector
    double Sdrxp[3];
    double Sdra[24],Sdrb[24],Sdrc[24],Sdrped[24];

    for ( int j=0; j<24; j++ )
    {
        nfound = s->Search ( hSd1r[j],2,"",0.1 );
       
        for ( int p=0; p<nfound; p++ ){Sdrxp[p] = txp[p];}
       
        Sdra[j] = min ( min ( Sdrxp[0],Sdrxp[1] ), Sdrxp[2] );
        Sdrc[j] = max ( max ( Sdrxp[0],Sdrxp[1] ), Sdrxp[2] );

    if(Sdrxp[0]!=Sdra[j] && Sdrxp[0]!=Sdrc[j]) {Sdrb[j] = Sdrxp[0];}
    else if(Sdrxp[1]!=Sdra[j] && Sdrxp[1]!=Sdrc[j]) {Sdrb[j] = Sdrxp[1];}
    else if(Sdrxp[2]!=Sdra[j] && Sdrxp[2]!=Sdrc[j]) {Sdrb[j] = Sdrxp[2];}            
    
    }
  
    //open a txt file to store the positions of the peaks
    ofstream Sd1r_alpha;
    Sd1r_alpha.open ( "/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1r_AlphaPeaks.txt" );
    Sd1r_alpha << " strip / 1 / 2 / 3  \n";
  
    ofstream Sd1r_alphaCounts;
    Sd1r_alphaCounts.open ( "/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1rCounts.txt" );
    Sd1r_alphaCounts << " strip / Counts 1 / 2 / 3  \n";
    for ( int i=0; i<24; i++ )
    {
      //cout << i << " a " << Sdra[i] << " b " << Sdrb[i] << " c " << Sdrc[i] << endl;
       if(i==0)
       {fit_func3->SetParLimits ( 1,Sdra[i]-5,Sdra[i]+5 );
        fit_func3->SetParLimits ( 4,Sdrb[i]-5,Sdrb[i]+5 );
        fit_func3->SetParLimits ( 6,Sdrc[i]-5,Sdrc[i]+5 );
        fit_func3->SetParLimits ( 2,0,20 );
       }
       else if(i>0)
       {fit_func3->SetParLimits ( 1,Sdra[i]-5,Sdra[i]+5 );
        fit_func3->SetParLimits ( 4,Sdrb[i]-5,Sdrb[i]+5 );
        fit_func3->SetParLimits ( 6,Sdrc[i]-5,Sdrc[i]+5 );
        fit_func3->SetParLimits ( 2,0,10 );
        }
        hSd1r[i]->Fit ( "fit_func3","","",Sdra[i]-50,Sdrc[i]+50 );

        Sd1r_alpha << i << " "  << fit_func3->GetParameter ( 1 ) << " " << fit_func3->GetParameter ( 4 ) << " " << fit_func3->GetParameter ( 6 ) << endl;
        Sd1r_alphaCounts << i << " " << fit_func3->GetParameter ( 0 ) << " " <<  fit_func3->GetParameter ( 3 ) << " " <<  fit_func3->GetParameter ( 5 ) << endl;
    }//end of fit on Sd1r detector
    
    //find peaks and fit the Sd1s detector

    double Sdsxp[3];
    double Sdsa[32],Sdsb[32],Sdsc[32],Sdsped[32];

    for ( int j=0; j<32; j++ ){
        nfound = s->Search ( hSd1s[j],2,"",0.1 );
       
        for ( int p=0; p<nfound; p++ ){Sdsxp[p] = txp[p];}
        
        Sdsa[j] = min ( min ( Sdsxp[0],Sdsxp[1] ), Sdsxp[2] );
        Sdsc[j] = max ( max ( Sdsxp[0],Sdsxp[1] ), Sdsxp[2] );

        if(Sdsxp[0]!=Sdsa[j] && Sdsxp[0]!=Sdsc[j]) {Sdsb[j] = Sdsxp[0];}
        else if(Sdsxp[1]!=Sdsa[j] && Sdsxp[1]!=Sdsc[j]) {Sdsb[j] = Sdsxp[1];}
        else if(Sdsxp[2]!=Sdsa[j] && Sdsxp[2]!=Sdsc[j]) {Sdsb[j] = Sdsxp[2];}       
                           
    }
  
    //open a txt file to store the positions of the peaks
    ofstream Sd1s_alpha;
    Sd1s_alpha.open ( "/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1s_AlphaPeaks.txt" );
    Sd1s_alpha << " strip / 1 / 2 / 3  \n";
  
    ofstream Sd1s_alphaCounts;
    Sd1s_alphaCounts.open ( "/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1sCounts.txt" );
    Sd1s_alphaCounts << " strip / Counts 1 / 2 / 3  \n";
    for ( int i=0; i<32; i++ ){
        fit_func3->SetParLimits ( 1,Sdsa[i]-10,Sdsa[i]+10 );
        fit_func3->SetParLimits ( 4,Sdsb[i]-10,Sdsb[i]+10 );
        fit_func3->SetParLimits ( 6,Sdsc[i]-10,Sdsc[i]+10 );
        fit_func3->SetParLimits ( 2,0,10 );

        
        hSd1s[i]->Fit ( "fit_func3","","",Sdsa[i]-50,Sdsc[i]+50 );

        Sd1s_alpha << i << " "  << fit_func3->GetParameter ( 1 ) << " " << fit_func3->GetParameter ( 4 ) << " " << fit_func3->GetParameter ( 6 ) << endl;
        Sd1s_alphaCounts << i << " " << fit_func3->GetParameter ( 0 ) << " " <<  fit_func3->GetParameter ( 3 ) << " " <<  fit_func3->GetParameter ( 5 ) << endl;
    }//end of fit on Sd1s detector
        
    TFile *f_out = new TFile ( "/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1_AlphaPeaks_withFit.root","RECREATE" );
    
    f_out->cd();

    for ( int i=0; i<24; i++ ){
        hSd1r[i]->GetXaxis()->SetRangeUser(600,1000);
        hSd1r[i]->Write();
    }
    
    for ( int i=0; i<32; i++ ){
        hSd1s[i]->GetXaxis()->SetRangeUser(600,1000);
        hSd1s[i]->Write();   
    }
    
    f_out->Close();

    return 0;
}
