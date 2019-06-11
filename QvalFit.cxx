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

void QvalFit(){
	TFile *f = new TFile("/home/jerome/12Be_exp/Analysis/CarbonGain/C_pedestal_TRIUMF_DL_CarbonGain_NewDecode_NewCut_Yu_Random.root","READ");
	TF1* fit_func1 = new TF1 ( "fit_func1","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",-10,10);
	
	ofstream qvalfitpar;
	qvalfitpar.open ( "/home/jerome/12Be_exp/Analysis/CarbonGain/TES_QValues.txt" ); //change the path and the name of the file acordingly
	//qvalfitpar.open ( "/home/jerome/12Be_exp/Analysis/CarbonGain/GS_QValues.txt" ); //change the path and the name of the file acordingly
	qvalfitpar << " No. 	 Amplitude 	 Mean 	 SD  	Deviation\n";
	double C_GSQval=2.722;
	double C_TESQval=-1.132;
	int shiftLength=21;
	double shift[shiftLength];
	double shiftInitial=-0.4;
	double deviation[shiftLength];
	for(int i=0;i<shiftLength;i++){
		shift[i]=shiftInitial;
		shiftInitial=shiftInitial+0.02;//-0.4:0.02:0
	}
		
	TF1 *fit_quad = new TF1 ("fit_quad","[0]+[1]*x+[2]*pow(x,2)",0,200);
	TGraphErrors *gr1 = new TGraphErrors();
	
	TCanvas *c1 = new TCanvas ( "c1" );
    c1->Divide ( 7,3 );
	
	for(int i=0;i<shiftLength;i++){
		TH1D *hQ = (TH1D*)f->Get(Form("hQval_%d",i));
		c1->cd(i+1);
		hQ->Draw ();
		hQ->Rebin(2);
		//GS
		/*fit_func1->SetParLimits ( 2,0.05,0.3);
		fit_func1->SetParLimits ( 0,9,10);
		hQ->Fit ( "fit_func1","","",1,4);
		hQ->GetXaxis()->SetRangeUser(1,4.5);
		deviation[i]=abs(fit_func1->GetParameter ( 1 )-C_GSQval);*/
		
		//TES
		fit_func1->SetParLimits ( 2,0.1,0.4);
		if(i>1){
			fit_func1->SetParLimits ( 0,25,40);
		}
		hQ->Fit ( "fit_func1","","",-3,1);
		
		hQ->GetXaxis()->SetRangeUser(-3,1);
		deviation[i]=abs(fit_func1->GetParameter ( 1 )-C_TESQval);
		
		qvalfitpar << i << "	 " << fit_func1->GetParameter ( 0 ) << " 	" <<  fit_func1->GetParameter ( 1 ) << "	 " <<  fit_func1->GetParameter ( 2 ) << " 	" << deviation[i] << endl;
	}
	int Number=0;
	for(int i=11;i<shiftLength;i++){
		gr1->SetPoint(Number,shift[i],deviation[i]);
		Number++;
	}
	TCanvas *c2 = new TCanvas ( "c2" );
	gr1->Draw("AP*");
	gr1->Fit(fit_quad,"QE","",shift[0],shift[20]);
	cout << -fit_quad->GetParameter(1)/(2*fit_quad->GetParameter(2)) << endl;
}