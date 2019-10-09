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
	float kBeam = 112.21; //put the correct value; beam energy
	float mbeam = 12 * amu;  //mass of the beam (12Be or 12C) in MeV
	float mrecoil = 13 * amu;  //mass of the recoil (13Be or 13C) in MeV
	float mejec = 1 * amu; //mass of the proton
	
	return pow(((TMath::Sqrt(mbeam*amu*kBeam)*TMath::Cos(TMath::DegToRad()*x[0])+TMath::Sqrt(mbeam*amu*kBeam*pow(TMath::Cos(TMath::DegToRad()*x[0]),2)+(mrecoil+amu)*(mrecoil*Q_C[0]+(mrecoil-mbeam)*kBeam)))/(mrecoil+amu)),2); 
} 

void kinematics(){
	double angle[17]={0};
    float TDistance = 80.88;
	float Qgs13C = 2.722; // Ground state Q value of 13C
	float QFES13C = -0.367;
    float QSES13C = -0.962;
	float Exenergy = 3.853; //Third excited state energy
	float QTES13C = Qgs13C-Exenergy;
	
	float Qgs13Be = -2.44; // Ground state Q value of 13Be
	float Sn = -0.2198; //Neutron separation energy
	float gs = 0.2;
    float fes = 0.8;
    float ses = 1.4;
	float Q13Be0 = Qgs13Be-Sn; // 0 MeV Sn
	//float Q13Be1 = Q13Be0-0.40;//-2.6202 (0.40)
	//float Q13Be2 = Qgs13Be;//-2.6802 (0.46)
	//float Q13Be3 = Qgs13Be-(fes-gs);//-3.0702 (0.85)
	//float Q13Be4 = Qgs13Be-(ses-gs);//-4.2202 (2);
	
	float Q13Be1 = -2.44;//Qgs13Be-(0.46-gs);//Relative to the threshold
	float Q13Be2 = -2.67;
	float Q13Be3 = -2.98;//Qgs13Be-(0.85-gs);
	float Q13Be4 = -4.16;//Qgs13Be-(2-gs);
	double Q2[2] = {Q13Be2,Q13Be2};
    double Q3[2] = {Q13Be3,Q13Be3};
    double Q4[2] = {Q13Be4,Q13Be4};
    
    double QC0[2] = {Qgs13C,Qgs13C};
    double QC1[2] = {QFES13C,QFES13C};
    double QC2[2] = {QSES13C,QSES13C};
    double QC3[2] = {QTES13C,QTES13C};
    double countRange[2] = {0,50};
    double angleRange[2] = {120,160};
	
	for (int i=0;i<17;i++){
		angle[i]=180-TMath::RadToDeg()*TMath::ATan((50+(16-i)*(4.94))/TDistance);
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
	TH2D *h_noP = (TH2D*)f_noP->Get("hYuAnPID");*/
	
	TFile *f_CMDL = new TFile("../../Analysis/with_pedestal/C_pedestal_MicronDL_60sumT_85mm.root","READ");
	TH2D *h_CMDL = (TH2D*)f_CMDL->Get("hYuAnPID");
	TH1D *h_CMDLQ = (TH1D*)f_CMDL->Get("hQval");
	
	TFile *f_CTDL = new TFile("../../Analysis/BeamOffset/C_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_Yu_TargetDistance80_88mm.root","READ");
	TH2D *h_CTDL = (TH2D*)f_CTDL->Get("hYuAnPID");
	TH1D *h_CTDLQ = (TH1D*)f_CTDL->Get("hQval");
    TH1D *h_CTDLQ0_7 = (TH1D*)f_CTDL->Get("hQvalRs0_7");
    TH1D *h_CTDLQ8_15 = (TH1D*)f_CTDL->Get("hQvalRs8_15");
    TH2D *h_CTDLQAn = (TH2D*)f_CTDL->Get("hQvalAn");
    TH1D *h_QvalR[16];
    for (int i=0;i<16;i++){
        h_QvalR[i] = (TH1D*)f_CTDL->Get(Form("hQvalR%d",i));;
    }
    TH1D *h_QvalS[8];
    for (int i=0;i<8;i++){
        h_QvalS[i] = (TH1D*)f_CTDL->Get(Form("hQvalS%d",i));;
    }
	
	TFile *f_BeMDL = new TFile("../../Analysis/Be/CarbonGain/Be_pedestal_Micron_DL_CarbonGain_Yu.root","READ");
	TH2D *h_BeMDL = (TH2D*)f_BeMDL->Get("hYuAnPID");
	TH1D *h_BeMDLQ = (TH1D*)f_BeMDL->Get("hQvalT");
	
	TFile *f_BeTDL = new TFile("../../Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod.root","READ");
	TH2D *h_BeTDL = (TH2D*)f_BeTDL->Get("hYuAnPID");
	TH1D *h_BeTDLQ = (TH1D*)f_BeTDL->Get("hQval");
    TH1D *h_BeTDLQ0_7 = (TH1D*)f_BeTDL->Get("hQval0_7");
    TH1D *h_BeTDLQ8_15 = (TH1D*)f_BeTDL->Get("hQval8_15");
    TH2D *h_BeTDLQAn = (TH2D*)f_BeTDL->Get("hQvalAn");
    TH1D *h_BeQvalR[16];
    for (int i=0;i<16;i++){
        h_BeQvalR[i] = (TH1D*)f_BeTDL->Get(Form("hQvalR%d",i));;
    }
    TH1D *h_BeQvalS[8];
    for (int i=0;i<8;i++){
        h_BeQvalS[i] = (TH1D*)f_BeTDL->Get(Form("hQvalS%d",i));;
    }
	
	TF1 *KinFcn = new TF1("KinFcn",kin,angle[0],angle[16],1);
	KinFcn->SetParameter(0,Qgs13C);
	TF1 *KinFcn1 = new TF1("KinFcn1",kin,angle[0],angle[16],1);
	KinFcn1->SetParameter(0,QFES13C);
    TF1 *KinFcn2 = new TF1("KinFcn2",kin,angle[0],angle[16],1);
	KinFcn2->SetParameter(0,QSES13C);
	TF1 *KinFcn33 = new TF1("KinFcn33",kin,angle[0],angle[16],1);
	KinFcn33->SetParameter(0,QTES13C);
	
	auto legend1 = new TLegend(0.7,0.7,0.9,0.9);
	legend1->AddEntry(KinFcn,"GS","l");
	legend1->AddEntry(KinFcn1,"1st 1/2+ (3.089 MeV)","l");
    legend1->AddEntry(KinFcn2,"2nd 3/2- (3.684 MeV)","l");
	legend1->AddEntry(KinFcn33,"3rd 5/2+ (3.853 MeV)","l");
	/*TCanvas *c1 = new TCanvas ( "c1" );
	//c1->Divide ( 2,3 );
	
	//c1->cd(1);
	//hCCalib_noP->SetTitle("Carbon calibration, no pedestal");
	//hCCalib_noP->Draw("colz");
	h_CMDL->SetTitle("12C(d,p)13C, Micron DL");
	h_CMDL->Draw("colz");
	KinFcn->SetLineColor(kGreen);
	KinFcn->SetLineWidth(3);
	KinFcn->Draw("same");
	KinFcn2->SetLineColor(kRed);
	KinFcn2->SetLineWidth(3);
	KinFcn2->Draw("same");
	legend1->Draw();*/
	
    
    
    //12C(d,p)13C energy vs. angle
	TCanvas *c2 = new TCanvas ( "c2" );
	
	//c1->cd(2);
	//hCCalib_P->SetTitle("Carbon calibration, with pedestal");
	//hCCalib_P->Draw("colz");
	h_CTDL->SetTitle("");//12C(d,p)13C, TRIUMF DL
	h_CTDL->Draw("colz");
    h_CTDL->GetXaxis()->SetTitle("Angle in Degrees"); h_CTDL->GetXaxis()->CenterTitle();
    h_CTDL->GetYaxis()->SetTitle("Proton Energy [MeV]"); h_CTDL->GetYaxis()->CenterTitle();
    h_CTDL->SetStats(0);
	KinFcn->SetLineColor(kGreen);
	KinFcn->SetLineWidth(3);
	KinFcn->Draw("same");
	KinFcn1->SetLineColor(kBlack);
	KinFcn1->SetLineWidth(3);
	KinFcn1->Draw("same");
	KinFcn2->SetLineColor(4);
	KinFcn2->SetLineWidth(3);
	KinFcn2->Draw("same");
    KinFcn33->SetLineColor(kRed);
	KinFcn33->SetLineWidth(3);
	KinFcn33->Draw("same");
	legend1->Draw();
	
    
    //12C(d,p)13C Q value spectrum compared with the actual values
    TCanvas *c3 = new TCanvas ( "c3" );
    h_CTDLQ->Draw();
    //h_CTDLQ->Rebin(2);
	h_CTDLQ->SetTitle("12C(d,p)13C Q value, TRIUMF DL");
    TGraph *gr2 = new TGraph(2,QC0,countRange);
    gr2->Draw("same");
    gr2->SetLineColor(kGreen);
    gr2->SetLineWidth(3);
    TGraph *gr3 = new TGraph(2,QC1,countRange);
    gr3->Draw("same");
    gr3->SetLineColor(kBlack);
    gr3->SetLineWidth(3);
    TGraph *gr22 = new TGraph(2,QC2,countRange);
    gr22->Draw("same");
    gr22->SetLineColor(4);
    gr22->SetLineWidth(3);
    TGraph *gr4 = new TGraph(2,QC3,countRange);
    gr4->Draw("same");
    gr4->SetLineColor(kRed);
    gr4->SetLineWidth(3);
    legend1->Draw();
    
    //Ring by ring
    TCanvas *c8 = new TCanvas ( "c8" );
    c8->Divide(4,4);
    for (int i=0;i<16;i++){
        c8->cd(i+1);
        h_QvalR[i]->Draw();
        gr2->Draw("same");
        gr2->SetLineColor(kGreen);
        gr2->SetLineWidth(3);
        gr3->Draw("same");
        gr3->SetLineColor(kBlack);
        gr3->SetLineWidth(3);
        gr4->Draw("same");
        gr4->SetLineColor(kRed);
        gr4->SetLineWidth(3);
        legend1->Draw();
    }
    
    //Rings 0 to 7 and 8 to 15
    TCanvas *c9 = new TCanvas ( "c9" );
	c9->Divide(1,2);
    c9->cd(1);
    h_CTDLQ0_7->Draw();
    gr2->Draw("same");
    gr2->SetLineColor(kGreen);
    gr2->SetLineWidth(3);
    gr3->Draw("same");
    gr3->SetLineColor(kBlack);
    gr3->SetLineWidth(3);
    gr4->Draw("same");
    gr4->SetLineColor(kRed);
    gr4->SetLineWidth(3);
    legend1->Draw();
    //legend2->Draw();
    c9->cd(2);
    h_CTDLQ8_15->Draw();
    gr2->Draw("same");
    gr2->SetLineColor(kGreen);
    gr2->SetLineWidth(3);
    gr3->Draw("same");
    gr3->SetLineColor(kBlack);
    gr3->SetLineWidth(3);
    gr4->Draw("same");
    gr4->SetLineColor(kRed);
    gr4->SetLineWidth(3);
    legend1->Draw();
    
    TCanvas *c10 = new TCanvas ( "c10" );
    c10->Divide(3,3);
    for (int i=0;i<8;i++){
        c10->cd(i+1);
        h_QvalS[i]->Draw();
        gr2->Draw("same");
        gr2->SetLineColor(kGreen);
        gr2->SetLineWidth(3);
        gr3->Draw("same");
        gr3->SetLineColor(kBlack);
        gr3->SetLineWidth(3);
        gr4->Draw("same");
        gr4->SetLineColor(kRed);
        gr4->SetLineWidth(3);
        legend1->Draw();
    }
    
    
	/*c1->cd(3);
	hCCalib_noP_1_2peaks->SetTitle("Carbon calibration, first and alpha second peaks, no pedestal");
	hCCalib_noP_1_2peaks->Draw("colz");
	KinFcn->SetLineColor(kGreen);
	KinFcn->SetLineWidth(3);
	KinFcn->Draw("same");
	KinFcn2->SetLineColor(kRed);
	KinFcn2->SetLineWidth(3);
	KinFcn2->Draw("same");
	
	
	
	//c1->cd(4);
	//hCCalib_P_1_2peaks->SetTitle("Carbon calibration, first and second alpha peaks, with pedestal");
	//hCCalib_P_1_2peaks->Draw("colz");
	KinFcn->SetLineColor(kGreen);
	KinFcn->SetLineWidth(3);
	KinFcn->Draw("same");
	KinFcn2->SetLineColor(kRed);
	KinFcn2->SetLineWidth(3);
	KinFcn2->Draw("same");*/
	
	/*c1->cd(5);
	h_noP->SetTitle("No pedestal");
	h_noP->Draw("colz");
	KinFcn->SetLineColor(kGreen);
	KinFcn->SetLineWidth(3);
	KinFcn->Draw("same");
	KinFcn2->SetLineColor(kRed);
	KinFcn2->SetLineWidth(3);
	KinFcn2->Draw("same");
	
	//c1->cd(6);
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
	
	/*TFile *f_P = new TFile("../../Analysis/Be/AlphaOnly/Be_pedestal.root","READ");
	TH2D *h_P = (TH2D*)f_P->Get("hYuAnPID1");
	
	TFile *fCCalib_P = new TFile("../../Analysis/Be/CarbCalibAlpha/Be_pedestal_CCalib.root","READ");
	TH2D *hCCalib_P = (TH2D*)fCCalib_P->Get("hYuAnPID1");
	
	TFile *fCCalib_noP_1_2peaks = new TFile("../../Analysis/Be/CarbCalibAlpha/Be_nopedestal_1_2Peaks_CCalib.root","READ");
	TH2D *hCCalib_noP_1_2peaks = (TH2D*)fCCalib_noP_1_2peaks->Get("hYuAnPID1");
	
	TFile *fCCalib_P_1_2peaks = new TFile("../../Analysis/Be/CarbCalibAlpha/Be_pedestal_1_2Peaks_CCalib.root","READ");
	TH2D *hCCalib_P_1_2peaks = (TH2D*)fCCalib_P_1_2peaks->Get("hYuAnPID1");*/
	
	TF1 *KinFcn0 = new TF1("KinFcn0",kin,angle[0],angle[16],1);
	KinFcn0->SetParameter(0,Q13Be0);
	TF1 *KinFcn3 = new TF1("KinFcn3",kin,angle[0],angle[16],1);
	KinFcn3->SetParameter(0,Q13Be1);
	TF1 *KinFcn4 = new TF1("KinFcn4",kin,angle[0],angle[16],1);
	KinFcn4->SetParameter(0,Q13Be2);
	TF1 *KinFcn5 = new TF1("KinFcn5",kin,angle[0],angle[16],1);
	KinFcn5->SetParameter(0,Q13Be3);
	TF1 *KinFcn6 = new TF1("KinFcn6",kin,angle[0],angle[16],1);
	KinFcn6->SetParameter(0,Q13Be4);
	
	auto legend2 = new TLegend(0.7,0.7,0.9,0.9);
	legend2->AddEntry(KinFcn0,"0.17 MeV","l");
	//legend2->AddEntry(KinFcn3,"0.40 MeV","l");
	legend2->AddEntry(KinFcn4,Form("%.2f MeV",0.42),"l");//Q13Be2
	legend2->AddEntry(KinFcn5,Form("%.2f MeV",0.76),"l");//Q13Be3
	legend2->AddEntry(KinFcn6,Form("%.2f MeV",1.91),"l");//Q13Be4
	
	//TCanvas *c3 = new TCanvas ( "c3" );
	//c1->Divide ( 2,2 );
	
	/*c1->cd(1);
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
	KinFcn3->SetLineWidth(3);
	KinFcn3->Draw("same");*/
	KinFcn3->SetLineColor(kBlack);
	
	
	
	//c1->cd(3);
	/*hCCalib_P_1_2peaks->SetTitle("Carbon calibration, first and second alpha peaks, with pedestal");
	hCCalib_P_1_2peaks->Draw("colz");
	hCCalib_P_1_2peaks->Rebin2D(1,10);
	hCCalib_P_1_2peaks->GetYaxis()->SetRangeUser(0.6,2.4);*/
	/*h_BeMDL->SetTitle("12Be(d,p)13Be, Micron DL");
	h_BeMDL->Draw("colz");
	//h_BeMDL->Rebin2D(1,10);
	//h_BeMDL->GetYaxis()->SetRangeUser(0.6,2.4);
	KinFcn0->SetLineColor(kYellow);
	KinFcn0->SetLineWidth(3);
	KinFcn0->Draw("same");
	KinFcn3->SetLineColor(kGreen);
	KinFcn3->SetLineWidth(3);
	KinFcn3->Draw("same");
	KinFcn4->SetLineColor(kRed);
	KinFcn4->SetLineWidth(3);
	KinFcn4->Draw("same");
	KinFcn5->SetLineColor(kBlack);
	KinFcn5->SetLineWidth(3);
	KinFcn5->Draw("same");*/
	/*KinFcn6->SetLineColor(2);
	KinFcn6->SetLineWidth(3);
	KinFcn6->Draw("same");*/
	//legend2->Draw();
	
	TCanvas *c4 = new TCanvas ( "c4" );
	
	//c1->cd(4);
	/*hCCalib_noP_1_2peaks->SetTitle("Carbon calibration, first and second alpha peaks, without pedestal");
	hCCalib_noP_1_2peaks->Draw("colz");
	hCCalib_noP_1_2peaks->GetYaxis()->SetRangeUser(0.6,2.4);*/
	h_BeTDL->SetTitle("");//12Be(d,p)13Be, TRIUMF DL
	h_BeTDL->Draw("colz");
    h_BeTDL->GetXaxis()->SetTitle("Angle in Degrees"); h_BeTDL->GetXaxis()->CenterTitle();
    h_BeTDL->GetYaxis()->SetTitle("Proton Energy [MeV]"); h_BeTDL->GetYaxis()->CenterTitle();
    h_BeTDL->SetStats(0);
	//h_BeTDL->Rebin2D(1,2);
	h_BeTDL->GetYaxis()->SetRangeUser(0.6,2.4);
	KinFcn0->SetLineColor(kYellow);
	KinFcn0->SetLineWidth(3);
	KinFcn0->Draw("same");
	KinFcn3->SetLineColor(kBlack);
	KinFcn3->SetLineWidth(3);
	//KinFcn3->Draw("same");
	KinFcn4->SetLineColor(kRed);
	KinFcn4->SetLineWidth(3);
	KinFcn4->Draw("same");
	KinFcn5->SetLineColor(kGreen);
	KinFcn5->SetLineWidth(3);
	KinFcn5->Draw("same");
	KinFcn6->SetLineColor(6);
	KinFcn6->SetLineWidth(3);
	KinFcn6->Draw("same");
	
	legend2->Draw();
	
	TCanvas *c5 = new TCanvas ( "c5" );
	/*c5->Divide(1,2);
	
	c5->cd(1);
	h_CMDLQ->Draw();
	h_CMDLQ->SetTitle("12C(d,p)13C Q value, Micron DL");
	
	c5->cd(2);
	h_CTDLQ->Draw();
	h_CTDLQ->SetTitle("12C(d,p)13C Q value, TRIUMF DL");
	
	TCanvas *c6 = new TCanvas ( "c6" );
	c6->Divide(1,2);
	
	c6->cd(1);
	h_BeMDLQ->Draw();
	
	h_BeMDLQ->GetXaxis()->SetRangeUser(-4.4,-1);
	h_BeMDLQ->SetTitle("12Be(d,p)13Be Q value, Micron DL");*/
	
    
    //Q values for 12Be(d,p)13Be, commented out on 31/07/2019 from line 448 to 508
	//c6->cd(2);
	h_BeTDLQ->Draw();
    h_BeTDLQ->Rebin(4);
	h_BeTDLQ->GetXaxis()->SetRangeUser(-5,-2);
	h_BeTDLQ->SetTitle("12Be(d,p)13Be Q value, TRIUMF DL");
    TGraph *gr5 = new TGraph(2,Q2,countRange);
    gr5->Draw("same");
    gr5->SetLineColor(kRed);
    gr5->SetLineWidth(3);
    TGraph *gr6 = new TGraph(2,Q3,countRange);
    gr6->Draw("same");
    gr6->SetLineColor(kGreen);
    gr6->SetLineWidth(3);
    TGraph *gr7 = new TGraph(2,Q4,countRange);
    gr7->Draw("same");
    gr7->SetLineColor(6);
    gr7->SetLineWidth(3);
    legend2->Draw();
    
    TCanvas *c6 = new TCanvas ( "c6" );
	h_BeTDLQAn->Draw("colz");
    TGraph *grAn2 = new TGraph(2,angleRange,Q2);
    grAn2->Draw("same");
    grAn2->SetLineColor(kRed);
    grAn2->SetLineWidth(3);
    TGraph *grAn3 = new TGraph(2,angleRange,Q3);
    grAn3->Draw("same");
    grAn3->SetLineColor(kGreen);
    grAn3->SetLineWidth(3);
    TGraph *grAn4 = new TGraph(2,angleRange,Q4);
    grAn4->Draw("same");
    grAn4->SetLineColor(6);
    grAn4->SetLineWidth(3);
    legend2->Draw();
    
    TCanvas *c7 = new TCanvas ( "c7" );
	c7->Divide(1,2);
    c7->cd(1);
    h_BeTDLQ0_7->Draw();
    h_BeTDLQ0_7->Rebin(4);
    gr5->Draw("same");
    gr5->SetLineColor(kRed);
    gr5->SetLineWidth(3);
    gr6->Draw("same");
    gr6->SetLineColor(kGreen);
    gr6->SetLineWidth(3);
    gr7->Draw("same");
    gr7->SetLineColor(6);
    gr7->SetLineWidth(3);
    //legend2->Draw();
    c7->cd(2);
    h_BeTDLQ8_15->Draw();
    h_BeTDLQ8_15->Rebin(4);
    gr5->Draw("same");
    gr5->SetLineColor(kRed);
    gr5->SetLineWidth(3);
    gr6->Draw("same");
    gr6->SetLineColor(kGreen);
    gr6->SetLineWidth(3);
    gr7->Draw("same");
    gr7->SetLineColor(6);
    gr7->SetLineWidth(3);
    
    //Ring by ring
    TCanvas *c11 = new TCanvas ( "c11" );
    c11->Divide(4,4);
    for (int i=0;i<16;i++){
        c11->cd(i+1);
        h_BeQvalR[i]->Draw();
        h_BeQvalR[i]->Rebin(4);
        h_BeQvalR[i]->GetXaxis()->SetRangeUser(-5,-2);
        gr5->Draw("same");
        gr5->SetLineColor(kRed);
        gr5->SetLineWidth(3);
        gr6->Draw("same");
        gr6->SetLineColor(kGreen);
        gr6->SetLineWidth(3);
        gr7->Draw("same");
        gr7->SetLineColor(6);
        gr7->SetLineWidth(3);
        legend2->Draw();
    }
    

}
