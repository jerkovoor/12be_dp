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

void FitBeamOffset(){
	TFile *f = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/C_pedestal_TRIUMF_DL_BeamOffset5_5_shift0_091_TargetDistance81_2.root","READ");
	TF1* fit_func1 = new TF1 ( "fit_func1","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",-10,10);
	
    TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/C_pedestal_TRIUMF_DL_BeamOffset5_5_shift0_091_TargetDistance81_2_withFit.root","RECREATE");
    
	ofstream qvalfitpar;
	//qvalfitpar.open ( "/home/jerome/12Be_exp/Analysis/CarbonGain/TES_QValues.txt" ); //change the path and the name of the file acordingly
	qvalfitpar.open ( "/home/jerome/12Be_exp/Analysis/BeamOffset/GS_QValues_BeamOffset5_5_shift0_091_TargetDistance81_2.txt" ); //change the path and the name of the file acordingly
	//qvalfitpar << " x      y    	 Amplitude 	 Mean 	 SD  	Deviation\n";
	double OffsetInitialx=-5;
    double OffsetFinalx=5;
    double OffsetInitialy=-5;
    double OffsetFinaly=5;
    double interval=1;
    int OffsetLengthx=1+(OffsetFinalx-OffsetInitialx)/interval;
    int OffsetLengthy=1+(OffsetFinaly-OffsetInitialy)/interval;
    int Length=OffsetLengthx*OffsetLengthy;

	TCanvas *c1 = new TCanvas ( "c1" );
    c1->Divide ( OffsetLengthy,OffsetLengthx);
    
    f_out->cd();
	
	for(int i=0;i<Length;i++){
		TH1D *hQ = (TH1D*)f->Get(Form("hQval_%d",i));
		c1->cd(i+1);
		//hQ->Draw ();
		hQ->Rebin(2);
		//GS
		fit_func1->SetParLimits ( 2,0.05,1.5);
		fit_func1->SetParLimits ( 0,1,25);
		hQ->Fit ( "fit_func1","QE","",1,4);
		hQ->GetXaxis()->SetRangeUser(1.5,4.5);
        string name2 = Form ( "hQval_%i",i );
        hQ->Write(name2.c_str());
		
		//TES
		/*fit_func1->SetParLimits ( 2,0.1,0.5);
		if(i>1){
			fit_func1->SetParLimits ( 0,25,80);
		}
		hQ->Fit ( "fit_func1","","",-3,1);
		
		hQ->GetXaxis()->SetRangeUser(-3,1);
		deviation[i]=abs(fit_func1->GetParameter ( 1 )-C_TESQval);*/
        
        
		
		//qvalfitpar << OffsetInitialx+((i-i%OffsetLengthx)/OffsetLengthx)*interval << "	 " << OffsetInitialy-(i%OffsetLengthy)*interval << "	 " << fit_func1->GetParameter ( 0 ) << " 	" <<  fit_func1->GetParameter ( 1 ) << "	 " <<  fit_func1->GetParameter ( 2 ) << endl;
		qvalfitpar << i << "	 " << OffsetInitialx+((i-i%OffsetLengthy)/OffsetLengthy)*interval << "	 " << OffsetInitialy+(i%OffsetLengthy)*interval << "	 " <<  fit_func1->GetParameter ( 2 ) << endl;
	}
}
