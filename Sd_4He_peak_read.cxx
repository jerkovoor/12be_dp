//4He in-beam calibration of Sd detectors

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

double Sd_4He_peak_read()
{
    
    TFile *f = new TFile ( "/home/jerome/12Be_exp/Analysis/Sd_Calibration/C_target.root","READ" );
    
    TH1D* hSd1r[24];
    TH1D* hSd1s[32];
    TH1D* hSd2r[24];
    TH1D* hSd2s[32];
    
    
    for (int i=0;i<24;i++){
        hSd1r[i] = (TH1D*)f->Get(Form("Sd1r_Ch%i",i));
        hSd2r[i] = (TH1D*)f->Get(Form("Sd2r_Ch%i",i));
    }
    
    for (int i=0;i<32;i++){
        hSd1s[i] = (TH1D*)f->Get(Form("Sd1s_Ch%i",i));
        hSd2s[i] = (TH1D*)f->Get(Form("Sd2s_Ch%i",i));
    }
    
    
        
    TF1* fit_func1 = new TF1 ( "fit_func1","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",500,5000 ); //removed the liner part for comparison
    TF1* fit_func3 = new TF1 ( "fit_func3","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))",0,5000 );
    TF1* fit_func9 = new TF1 ( "fit_func9","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp(-pow((x-[10]),2)/(2*[2]*[2]))+[11]*TMath::Exp(-pow((x-[12]),2)/(2*[2]*[2]))+[13]*TMath::Exp(-pow((x-[14]),2)/(2*[2]*[2]))+[15]*TMath::Exp(-pow((x-[16]),2)/(2*[2]*[2]))+[17]*TMath::Exp(-pow((x-[18]),2)/(2*[2]*[2]))+[19]+[20]*x",2000,5000 );
    
    //Search for peaks in Sd strips
    int npeaks=1;
    TSpectrum* s = new TSpectrum ( 2*npeaks );
    int nfound=0;
    double* txp = s->GetPositionX();
    
    
    //////////////////////////
    ///        Sd1r        ///
    //////////////////////////
    
    //find peaks and fit the Sd1r detector
    double Sd1rxp[3];
    double Sd1ra[24],Sd1rped[24];

    for ( int j=0; j<24; j++ ){
        nfound = s->Search ( hSd1r[j],2,"",0.1 );
       
        for ( int p=0; p<nfound; p++ ){Sd1rxp[p] = txp[p];}
       
        Sd1ra[j] = Sd1rxp[0];
    
    }

  
    //open a txt file to store the positions of the peaks
    ofstream Sd1r_He;
    Sd1r_He.open ( "/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1r_4He_Peaks.txt" );
    Sd1r_He << " strip / 1  \n";
  
    ofstream Sd1r_HeCounts;
    Sd1r_HeCounts.open ( "/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1r_4He_Counts.txt" );
    Sd1r_HeCounts << " strip / Counts  \n";
    
    //fit_func1->SetParLimits ( 1,135,145 );

    //fit_func1->SetParLimits ( 2,0,40);
    //hSd1r[0]->Fit ( "fit_func1","Q","",100,180);
    
    for ( int i=0; i<24; i++ ){
        fit_func1->SetParLimits ( 1,Sd1ra[i]-5,Sd1ra[i]+5 );

        fit_func1->SetParLimits ( 2,0,30);
        hSd1r[i]->Fit ( "fit_func1","Q","",Sd1ra[i]-50,Sd1ra[i]+50 );

        Sd1r_He << i << " "  << fit_func1->GetParameter ( 1 ) << endl;
        Sd1r_HeCounts << i << " " << fit_func1->GetParameter ( 0 ) << endl;
    }//end of fit on Sd1r detector
    
    
    //////////////////////////
    ///        Sd1s        ///
    //////////////////////////
    

    double Sd1sxp[3];
    double Sd1sa[32],Sd1sped[32];

    for ( int j=0; j<32; j++ ){
        nfound = s->Search ( hSd1s[j],2,"",0.1 );
       
        for ( int p=0; p<nfound; p++ ){Sd1sxp[p] = txp[p];}
        
        Sd1sa[j] = Sd1sxp[0];
                           
    }
  
    //open a txt file to store the positions of the peaks
    ofstream Sd1s_alpha;
    Sd1s_alpha.open ( "/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1s_4He_Peaks.txt" );
    Sd1s_alpha << " strip / 1   \n";
  
    ofstream Sd1s_alphaCounts;
    Sd1s_alphaCounts.open ( "/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1s_4He_Counts.txt" );
    Sd1s_alphaCounts << " strip / Counts  \n";
    for ( int i=0; i<32; i++ ){
        fit_func1->SetParLimits ( 1,Sd1sa[i]-5,Sd1sa[i]+5 );
    
        fit_func1->SetParLimits ( 2,0,30 );

        
        hSd1s[i]->Fit ( "fit_func1","Q","",Sd1sa[i]-50,Sd1sa[i]+50 );

        Sd1s_alpha << i << " "  << fit_func1->GetParameter ( 1 ) << endl;
        Sd1s_alphaCounts << i << " " << fit_func1->GetParameter ( 0 ) << endl;
    }//end of fit on Sd1s detector
    
    
    //////////////////////////
    ///        Sd2r        ///
    //////////////////////////
    
    //find peaks and fir the Sd1r detector
    double Sd2rxp[3];
    double Sd2ra[24],Sd2rped[24];

    for ( int j=0; j<24; j++ ){
        nfound = s->Search ( hSd2r[j],2,"",0.1 );
       
        for ( int p=0; p<nfound; p++ ){Sd2rxp[p] = txp[p];}
       
        Sd2ra[j] = Sd2rxp[0];
    
    }

  
    //open a txt file to store the positions of the peaks
    ofstream Sd2r_He;
    Sd2r_He.open ( "/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd2r_4He_Peaks.txt" );
    Sd2r_He << " strip / 1  \n";
  
    ofstream Sd2r_HeCounts;
    Sd2r_HeCounts.open ( "/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd2r_4He_Counts.txt" );
    Sd2r_HeCounts << " strip / Counts  \n";
    
    //fit_func1->SetParLimits ( 1,135,145 );

    //fit_func1->SetParLimits ( 2,0,40);
    //hSd1r[0]->Fit ( "fit_func1","Q","",100,180);
    
    for ( int i=0; i<24; i++ ){
        fit_func1->SetParLimits ( 0,0,10 );
        fit_func1->SetParLimits ( 1,Sd2ra[i]-5,Sd2ra[i]+5 );

        fit_func1->SetParLimits ( 2,0,30);
        hSd2r[i]->Fit ( "fit_func1","Q","",Sd2ra[i]-50,Sd2ra[i]+50 );

        Sd2r_He << i << " "  << fit_func1->GetParameter ( 1 ) << endl;
        Sd2r_HeCounts << i << " " << fit_func1->GetParameter ( 0 ) << endl;
    }//end of fit on Sd2r detector
    
    
    //////////////////////////
    ///        Sd2s        ///
    //////////////////////////
    
    double Sd2sxp[3];
    double Sd2sa[32],Sd2sped[32];

    for ( int j=0; j<32; j++ ){
        nfound = s->Search ( hSd2s[j],2,"",0.1 );
       
        for ( int p=0; p<nfound; p++ ){Sd2sxp[p] = txp[p];}
        
        Sd2sa[j] = Sd2sxp[0];
                           
    }
  
    //open a txt file to store the positions of the peaks
    ofstream Sd2s_He;
    Sd2s_He.open ( "/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd2s_4He_Peaks.txt" );
    Sd2s_He << " strip / 1   \n";
  
    ofstream Sd2s_HeCounts;
    Sd2s_HeCounts.open ( "/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd2s_4He_Counts.txt" );
    Sd2s_HeCounts << " strip / Counts  \n";
    for ( int i=0; i<32; i++ ){
        fit_func1->SetParLimits ( 1,Sd2sa[i]-5,Sd2sa[i]+5 );
    
        fit_func1->SetParLimits ( 2,0,30 );

        
        hSd2s[i]->Fit ( "fit_func1","Q","",Sd2sa[i]-50,Sd2sa[i]+50 );

        Sd2s_He << i << " "  << fit_func1->GetParameter ( 1 ) << endl;
        Sd2s_HeCounts << i << " " << fit_func1->GetParameter ( 0 ) << endl;
    }//end of fit on Sd2s detector
    
    
        
    TFile *f_out = new TFile ( "/home/jerome/12Be_exp/Analysis/Sd_Calibration/C_target_withFit.root","RECREATE" );
    
    f_out->cd();

    for ( int i=0; i<24; i++ ){
        hSd1r[i]->GetXaxis()->SetRangeUser(0,300);
        hSd1r[i]->Write();
    }
    
    for ( int i=0; i<32; i++ ){
        //hSd1s[i]->GetXaxis()->SetRangeUser(800,1600);
        hSd1s[i]->Write();   
    }
    
    for ( int i=0; i<24; i++ ){
        //hSd1r[i]->GetXaxis()->SetRangeUser(0,300);
        hSd2r[i]->Write();
    }
    
    for ( int i=0; i<32; i++ ){
        //hSd1s[i]->GetXaxis()->SetRangeUser(800,1600);
        hSd2s[i]->Write();   
    }
    
    f_out->Close();

    return 0;
}
