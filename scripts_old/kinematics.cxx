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



double kin(double *x, double *Q_C){ 
	
	float amu = 931.5; // atomic mass unit in MeV 
	float massEjec = 938.28;//mass of the proton in MeV/c2 
	float kBeam = 109.98; //put the correct value; beam energy
	float mbeam = 12 * amu;  //mass of the beam (12Be or 12C) in MeV
	float mrecoil = 13 * amu;  //mass of the recoil (13Be or 13C) in MeV
	float mejec = 1 * amu; //mass of the proton
	
	return pow(((TMath::Sqrt(mbeam*amu*kBeam)*TMath::Cos(TMath::DegToRad()*x[0])+TMath::Sqrt(mbeam*amu*kBeam*pow(TMath::Cos(TMath::DegToRad()*x[0]),2)+(mrecoil+amu)*(mrecoil*Q_C[0]+(mrecoil-mbeam)*kBeam)))/(mrecoil+amu)),2); 
} 

void kinematics(){
	double angle[17]={0};
	float Qgs13C = 2.722; // Ground state Q value of 13C
	float Exenergy = 3.853; //Third excited state energy
	float QTES13C = Qgs13C-Exenergy;
	
	float Qgs13Be = -2.448; // Ground state Q value of 13Be
	float Sn = -0.2283; //Neutron separation energy
	float Q13Be1 = Qgs13Be-(Sn+0.2283);
	float Q13Be2 = Qgs13Be-(Sn+0.85);
	float Q13Be3 = Qgs13Be-(Sn+1.9345);
	
	for (int i=0;i<17;i++){
		angle[i]=180-TMath::RadToDeg()*TMath::ATan((50+(15-i)*(4.94))/85);
		//cout << angle[i] << endl;
	}
	
	////////////////////////
	////	CARBON		////
	////////////////////////
	
	/*TFile *fCCalib_noP = new TFile("../../Analysis/C_calib/C_nopedestal_CCalib.root","READ");
	TH2D *hCCalib_noP = (TH2D*)fCCalib_noP->Get("hYuAnPID");
	
	TFile *fCCalib_P = new TFile("../../Analysis/C_calib/C_pedestal_CCalib.root","READ");
	TH2D *hCCalib_P = (TH2D*)fCCalib_P->Get("hYuAnPID");
	
	TFile *fCCalib_P_1_2peaks = new TFile("../../Analysis/C_calib/C_pedestal_1_2Peaks_CCalib.root","READ");
	TH2D *hCCalib_P_1_2peaks = (TH2D*)fCCalib_P_1_2peaks->Get("hYuAnPID");
	
	TFile *fCCalib_noP_1_2peaks = new TFile("../../Analysis/C_calib/C_nopedestal_1_2Peaks_CCalib.root","READ");
	TH2D *hCCalib_noP_1_2peaks = (TH2D*)fCCalib_noP_1_2peaks->Get("hYuAnPID");
	
	TFile *f_noP = new TFile("../../Analysis/no_pedestal/C_nopedestal.root","READ");
	TH2D *h_noP = (TH2D*)f_noP->Get("hYuAnPID");
	
	TFile *f_P = new TFile("../../Analysis/with_pedestal/C_pedestal_MicronDL_60sumT_85mm.root","READ");
	TH2D *h_P = (TH2D*)f_P->Get("hYuAnPID");
	
	TF1 *KinFcn = new TF1("KinFcn",kin,angle[0],angle[16],1);
	KinFcn->SetParameter(0,Qgs13C);
	TF1 *KinFcn2 = new TF1("KinFcn2",kin,angle[0],angle[16],1);
	KinFcn2->SetParameter(0,QTES13C);
	
	TCanvas *c1 = new TCanvas ( "c1" );
	c1->Divide ( 2,3 );
	
	c1->cd(1);
	hCCalib_noP->SetTitle("Carbon calibration, no pedestal");
	hCCalib_noP->Draw("colz");
	KinFcn->SetLineColor(kGreen);
	KinFcn->SetLineWidth(3);
	KinFcn->Draw("same");
	KinFcn2->SetLineColor(kRed);
	KinFcn2->SetLineWidth(3);
	KinFcn2->Draw("same");
	
	c1->cd(2);
	hCCalib_P->SetTitle("Carbon calibration, with pedestal");
	hCCalib_P->Draw("colz");
	KinFcn->SetLineColor(kGreen);
	KinFcn->SetLineWidth(3);
	KinFcn->Draw("same");
	KinFcn2->SetLineColor(kRed);
	KinFcn2->SetLineWidth(3);
	KinFcn2->Draw("same");
	
	c1->cd(3);
	hCCalib_noP_1_2peaks->SetTitle("Carbon calibration, first and alpha second peaks, no pedestal");
	hCCalib_noP_1_2peaks->Draw("colz");
	KinFcn->SetLineColor(kGreen);
	KinFcn->SetLineWidth(3);
	KinFcn->Draw("same");
	KinFcn2->SetLineColor(kRed);
	KinFcn2->SetLineWidth(3);
	KinFcn2->Draw("same");
	
	c1->cd(4);
	hCCalib_P_1_2peaks->SetTitle("Carbon calibration, first and second alpha peaks, with pedestal");
	hCCalib_P_1_2peaks->Draw("colz");
	KinFcn->SetLineColor(kGreen);
	KinFcn->SetLineWidth(3);
	KinFcn->Draw("same");
	KinFcn2->SetLineColor(kRed);
	KinFcn2->SetLineWidth(3);
	KinFcn2->Draw("same");
	
	c1->cd(5);
	h_noP->SetTitle("No pedestal");
	h_noP->Draw("colz");
	KinFcn->SetLineColor(kGreen);
	KinFcn->SetLineWidth(3);
	KinFcn->Draw("same");
	KinFcn2->SetLineColor(kRed);
	KinFcn2->SetLineWidth(3);
	KinFcn2->Draw("same");
	
	c1->cd(6);
	h_P->SetTitle("With pedestal");
	h_P->Draw("colz");
	KinFcn->SetLineColor(kGreen);
	KinFcn->SetLineWidth(3);
	KinFcn->Draw("same");
	KinFcn2->SetLineColor(kRed);
	KinFcn2->SetLineWidth(3);
	KinFcn2->Draw("same");*/
	
	////////////////////////
	////	BERILIUM	////
	////////////////////////
	
	TFile *f_P = new TFile("../../Analysis/Be/AlphaOnly/Be_pedestal.root","READ");
	TH2D *h_P = (TH2D*)f_P->Get("hYuAnPID");
	
	TFile *fCCalib_P = new TFile("../../Analysis/Be/CarbCalibAlpha/Be_pedestal_CCalib.root","READ");
	TH2D *hCCalib_P = (TH2D*)fCCalib_P->Get("hYuAnPID");
	
	TFile *fCCalib_noP_1_2peaks = new TFile("../../Analysis/Be/CarbCalibAlpha/Be_nopedestal_1_2Peaks_CCalib.root","READ");
	TH2D *hCCalib_noP_1_2peaks = (TH2D*)fCCalib_noP_1_2peaks->Get("hYuAnPID");
	
	TFile *fCCalib_P_1_2peaks = new TFile("../../Analysis/Be/CarbCalibAlpha/Be_pedestal_1_2Peaks_CCalib.root","READ");
	TH2D *hCCalib_P_1_2peaks = (TH2D*)fCCalib_P_1_2peaks->Get("hYuAnPID");
	
	TF1 *KinFcn = new TF1("KinFcn",kin,angle[0],angle[16],1);
	KinFcn->SetParameter(0,Q13Be1);
	TF1 *KinFcn2 = new TF1("KinFcn2",kin,angle[0],angle[16],1);
	KinFcn2->SetParameter(0,Q13Be2);
	TF1 *KinFcn3 = new TF1("KinFcn3",kin,angle[0],angle[16],1);
	KinFcn3->SetParameter(0,Q13Be3);
	
	TCanvas *c1 = new TCanvas ( "c1" );
	c1->Divide ( 2,2 );
	
	c1->cd(1);
	h_P->SetTitle("With pedestal");
	h_P->Draw("colz");
	h_P->Rebin2D(1,10);
	h_P->GetYaxis()->SetRangeUser(0.6,2.4);
	KinFcn->SetLineColor(kGreen);
	KinFcn->SetLineWidth(3);
	KinFcn->Draw("same");
	KinFcn2->SetLineColor(kRed);
	KinFcn2->SetLineWidth(3);
	KinFcn2->Draw("same");
	KinFcn3->SetLineColor(kBlack);
	KinFcn3->SetLineWidth(3);
	KinFcn3->Draw("same");
	
	c1->cd(2);
	hCCalib_P->SetTitle("Carbon calibration, with pedestal");
	hCCalib_P->Draw("colz");
	hCCalib_P->Rebin2D(1,10);
	hCCalib_P->GetYaxis()->SetRangeUser(0.6,2.4);
	KinFcn->SetLineColor(kGreen);
	KinFcn->SetLineWidth(3);
	KinFcn->Draw("same");
	KinFcn2->SetLineColor(kRed);
	KinFcn2->SetLineWidth(3);
	KinFcn2->Draw("same");
	KinFcn3->SetLineColor(kBlack);
	KinFcn3->SetLineWidth(3);
	KinFcn3->Draw("same");
	
	c1->cd(3);
	hCCalib_P_1_2peaks->SetTitle("Carbon calibration, first and second alpha peaks, with pedestal");
	hCCalib_P_1_2peaks->Draw("colz");
	hCCalib_P_1_2peaks->Rebin2D(1,10);
	hCCalib_P_1_2peaks->GetYaxis()->SetRangeUser(0.6,2.4);
	KinFcn->SetLineColor(kGreen);
	KinFcn->SetLineWidth(3);
	KinFcn->Draw("same");
	KinFcn2->SetLineColor(kRed);
	KinFcn2->SetLineWidth(3);
	KinFcn2->Draw("same");
	KinFcn3->SetLineColor(kBlack);
	KinFcn3->SetLineWidth(3);
	KinFcn3->Draw("same");
	
	c1->cd(4);
	hCCalib_noP_1_2peaks->SetTitle("Carbon calibration, first and second alpha peaks, without pedestal");
	hCCalib_noP_1_2peaks->Draw("colz");
	hCCalib_noP_1_2peaks->GetYaxis()->SetRangeUser(0.6,2.4);
	KinFcn->SetLineColor(kGreen);
	KinFcn->SetLineWidth(3);
	KinFcn->Draw("same");
	KinFcn2->SetLineColor(kRed);
	KinFcn2->SetLineWidth(3);
	KinFcn2->Draw("same");
	KinFcn3->SetLineColor(kBlack);
	KinFcn3->SetLineWidth(3);
	KinFcn3->Draw("same");
}