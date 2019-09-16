using namespace std;

#include <TFile.h>
#include <TString.h>
#include <string>
#include <iostream>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TCutG.h>
#include <vector>
#include <TObject.h>
#include <TClass.h>
#include <TF1.h>
#include <fstream>
#include <TGraph.h>
#include <TCanvas.h>

void dataplot(){

    double angleGS[181], ampGS[181], errorGS[181], angleFES[181], ampFES[181], errorFES[181], angleSES[181], ampSES[181], errorSES[181], angleTES[181], ampTES[181], errorTES[181];
    
    ifstream CdpGS;
    CdpGS.open("/home/jerome/12Be_exp/scripts/TWOFNR/24.12CdpGS");
    
    if(CdpGS.is_open()){
        for(int i=0;i<181;i++){
        CdpGS >> angleGS[i] >> ampGS[i];// >> errorGS[i];
        cout << angleGS[i] << endl;
        }
    }
    
    ifstream CdpFES;
    CdpFES.open("/home/jerome/12Be_exp/scripts/TWOFNR/24.12CdpFES");
    
    if(CdpFES.is_open()){
        for(int i=0;i<181;i++){
        CdpFES >> angleFES[180-i] >> ampFES[180-i];// >> errorFES[180-i];
        }
    }
    
    ifstream CdpSES;
    CdpSES.open("/home/jerome/12Be_exp/scripts/TWOFNR/24.12CdpSES");
    
    if(CdpSES.is_open()){
        for(int i=0;i<181;i++){
        CdpSES >> angleSES[180-i] >> ampSES[180-i];// >> errorSES[180-i];
        }
    }
    
    ifstream CdpTES;
    CdpTES.open("/home/jerome/12Be_exp/scripts/TWOFNR/24.12CdpTES");
    
    if(CdpTES.is_open()){
        for(int i=0;i<181;i++){
        CdpTES >> angleTES[180-i] >> ampTES[180-i];// >> errorTES[180-i];
        }
    }
    
    
    
    TGraph *gr1;
    TGraph *gr2;
    TGraph *gr3;
    TGraph *gr4;
    
    
    TCanvas *c4 = new TCanvas ( "c4" );

    gr1 = new TGraph(181,angleGS,ampGS);
    gr1->SetLineColor(kYellow);
    gr1->SetLineWidth(3);
    gr1->Draw();
    
    gr2 = new TGraph(181,angleFES,ampFES);
    gr2->SetLineColor(kGreen);
    gr2->SetLineWidth(3);
    gr2->Draw("same");
    
    gr3 = new TGraph(181,angleSES,ampSES);
    gr3->SetLineColor(kBlue);
    gr3->SetLineWidth(3);
    gr3->Draw("same");
    
    gr4 = new TGraph(181,angleTES,ampTES);
    gr4->SetLineColor(kBlack);
    gr4->SetLineWidth(3);
    gr4->Draw("same");
    
    auto legend2 = new TLegend(0.7,0.7,0.9,0.9);
	legend2->AddEntry(gr1,"GS (1/2-)","l");
	legend2->AddEntry(gr2,"1st (1/2+)","l");
	legend2->AddEntry(gr3,"2nd (3/2-)","l");
	legend2->AddEntry(gr4,"3rd (5/2+)","l");
    legend2->Draw();
    
    
    //gr1->GetXaxis()->SetRangeUser(122,150);
    //gr1->GetYaxis()->SetRangeUser(0,1.5);
}
