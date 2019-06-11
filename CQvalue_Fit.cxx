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

void CQvalue_Fit(){
	TFile *f = new TFile("/home/jerome/12Be_exp/Analysis/CarbonGain/C_pedestal_TRIUMF_DL_QvalMinGS_Yu.root","READ");
	//TFile *f = new TFile("/home/jerome/12Be_exp/Analysis/Be/CarbonGain/Be_pedestal_TRIUMF_DL_CarbonGain_QvalMinGS_Random_NewDecode_Yu.root","READ");
	TF1* fit_func1a = new TF1 ( "fit_func1a","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",-2,-0.5);
	TF1* fit_func1b = new TF1 ( "fit_func1b","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",-0.6,0);
	TF1* fit_func1c = new TF1 ( "fit_func1c","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",2,3.5);
	//TF1* fit_func2 = new TF1 ( "fit_func2","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))",-10,10);
	TF1* fit_func3 = new TF1 ( "fit_func3","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))",-10,10);
	TH1D *hQ = (TH1D*)f->Get("hQval");
	hQ->Draw ();
	//hQ->Rebin(2);
	
	// 3-Gaussian Fit for Carbon
	fit_func3->SetParLimits ( 0,25,40);
	fit_func3->SetParLimits ( 1,-1.3,-0.5);
	//fit_func3->SetParLimits ( 3,8,14);
	fit_func3->SetParLimits ( 4,-0.5,0);
	fit_func3->SetParLimits ( 5,5,20);
	fit_func3->SetParLimits ( 6,2.1,3.2);
	fit_func3->SetParLimits ( 2,0,0.4);
	
	fit_func1a->SetParameters(35,-1.1,0.25);

	fit_func1b->SetParameters(8.5,-0.36,0.25);
	
	fit_func1c->SetParameters(10,2.7,0.25);
		
	hQ->Fit ( fit_func1a, "R");
	hQ->Fit ( fit_func1b, "R+");
	hQ->Fit ( fit_func1c, "R+");
	hQ->Fit ( "fit_func3","R+","",-3,4);
	hQ->GetXaxis()->SetRangeUser(-3,4);
	
	// 2-Gaussian Fit for Beryllium
	/*fit_func2->SetParLimits ( 0,20,24);
	fit_func2->SetParLimits ( 1,-4.2,-3.8);
	fit_func2->SetParLimits ( 3,8,14);
	fit_func2->SetParLimits ( 4,-3,-2.5);
	fit_func2->SetParLimits ( 2,0.15,0.25);
	
	hQ->Fit ( "fit_func2","","",-5,-1);
	hQ->GetXaxis()->SetRangeUser(-5,-1);*/
	
	
	// 3-Gaussian Fit for Beryllium
	/*fit_func3->SetParLimits ( 0,20,24);
	fit_func3->SetParLimits ( 1,-4.2,-3.8);
	fit_func3->SetParLimits ( 3,6,10);
	fit_func3->SetParLimits ( 4,-3.1,-3);
	fit_func3->SetParLimits ( 5,10,15);
	fit_func3->SetParLimits ( 6,-2.8,-2.4);
	fit_func3->SetParLimits ( 2,0.2,0.25);
	
	hQ->Fit ( "fit_func3","","",-5,-1);
	hQ->GetXaxis()->SetRangeUser(-5,-1);*/
}