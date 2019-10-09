//Calculates and plots the residual vs. the angle for different beam energies.
//Also plots the fit parameters of the residuals vs the beam energy and finds the minima for both the slope and the intercept.

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



double kin(double *x, double *kBeam){ 
	
	float amu = 931.5; // atomic mass unit in MeV 
	float massEjec = 938.28;//mass of the proton in MeV/c2 
	//float kBeam = 114; //put the correct value; beam energy
	float mbeam = 12 * amu;  //mass of the beam (12Be or 12C) in MeV
	float mrecoil = 13 * amu;  //mass of the recoil (13Be or 13C) in MeV
	float mejec = 1 * amu; //mass of the proton
	float Qgs13C = 2.722; // Ground state Q value of 13C
	//float Exenergy = 0;//GS
	float Exenergy = 3.853; //Third excited state energy=3.853 MeV
	float Q_C = Qgs13C-Exenergy;
	
	return pow(((TMath::Sqrt(mbeam*amu*kBeam[0])*TMath::Cos(TMath::DegToRad()*x[0])+TMath::Sqrt(mbeam*amu*kBeam[0]*pow(TMath::Cos(TMath::DegToRad()*x[0]),2)+(mrecoil+amu)*(mrecoil*Q_C+(mrecoil-mbeam)*kBeam[0])))/(mrecoil+amu)),2); 
} 

void projectionMulti2() {
	
	double angle[17]={0}, anglelow[16]={0}, anglehigh[16]={0}, mid_ang[16]={0},  xaxis[2]={0}, errorbar[16]={0};
	int BinLow[16]={0}, BinHigh[16]={0};
	double ExcitedStateCutEn[16][2]={{2,3},{2,3},{2,3},{1.8,2.8},{1.6,2.6},{1.6,2.6},{1.6,2.6},{1.6,2.6},{1.6,2.6},{1.4,2.4},{1.4,2.4},{1.2,2.4},{1.2,2.2},{1.3,2.1},{1.2,2.2},{1.2,2}};//for third excited state
	double ExcitedStateCutBin[16][2];
	double Values[16][40];
	int en_length=20;
	float KEBeam[en_length];
	//float KEBeamInitial=102;//for ground state
	float KEBeamInitial=94;//for third excited state
	for (int i=0;i<en_length;i++){
		KEBeam[i]=KEBeamInitial;
		KEBeamInitial=KEBeamInitial+2;
	}
	// Bin numbers 571 to 937 (energies from 2.85 to 4.68) were selected after inspecting all the 16 histograms for the ground state
	int MinBin = 571; //with pedestal
	//int MinBin = 583; //without pedestal (583 for alpha calibration and 578 for Carbon secondary calibration)
	int MaxBin = 937;
	TH1D *h1[16];
	TH1D *h2[16];
	TGraphErrors *gr1[en_length];
	
	TCanvas *c1 = new TCanvas ( "c1" ); //create a canvas
	//c1->Divide (en_length/2,2);
	c1->Divide (5,4);
	
	TCanvas *c2 = new TCanvas ( "c2" ); //create a canvas
	c2->Divide (1,2);
	
	TF1 *fit_lin = new TF1 ("fit_lin","pol1",0,200);
	TF1 *fit_quad = new TF1 ("fit_quad","[0]+[1]*x+[2]*pow(x,2)",0,200);
	double x[2]={120,160};
	TGraph *gr2 = new TGraph(2,x,xaxis);
	TGraphErrors *gr3 = new TGraphErrors();
	TGraphErrors *gr4 = new TGraphErrors();
	
	TFile *f = new TFile("../Analysis/C_calib/C_nopedestal_CCalib.root","READ");//with pedestal
	//TFile *f = new TFile("../Analysis/no_pedestal/C_nopedestal.root","READ");//without pedestal
	TH2D *h_Eloss = (TH2D*)f->Get("hYuAnPID");
	TH2D *h_noEloss = (TH2D*)f->Get("hYuAnPID1");
	
	for (int i=0;i<17;i++){
		angle[i]=180-TMath::RadToDeg()*TMath::ATan((50+(15-i)*(4.94))/85);        
	}
	
	double EnergyValue=0;
	
	for (int i=0;i<16;i++){
		anglelow[i]=angle[i]+0.01;
		anglehigh[i]=angle[i+1]-0.01;
		BinLow[i] = h_Eloss->GetXaxis()->FindBin(anglelow[i]);
		BinHigh[i] = h_Eloss->GetXaxis()->FindBin(anglehigh[i]);
		mid_ang[i]=(angle[i]+angle[i+1])/2;
		errorbar[i]=(angle[i+1]-angle[i])/2;
		h1[i] = h_Eloss->ProjectionY(Form("h1_%d",i),BinLow[i],BinHigh[i],"");
		h2[i] = h_noEloss->ProjectionY(Form("h2_%d",i),BinLow[i],BinHigh[i],"");
		ExcitedStateCutBin[i][0]=h1[i]->GetXaxis()->FindBin(ExcitedStateCutEn[i][0]);
		ExcitedStateCutBin[i][1]=h1[i]->GetXaxis()->FindBin(ExcitedStateCutEn[i][1]);
		//ExcitedStateCutBin[i][1]=h1[i]->GetXaxis()->FindBin((ExcitedStateCutEn[i][0]+ExcitedStateCutEn[i][1])/2);
	}
	
	ofstream parameters;
	parameters.open("/home/jerome/12Be_exp/Analysis/with_pedestal/fom1.txt"); //open a .txt file to store the results of the fit; change the path and name acordingly
	parameters << "p0\tp0_error\tp1\tp1_error \n";
	
	for (int KE=0;KE<en_length;KE++){
		gr1[KE] = new TGraphErrors();
		double kin_value[16]={0}, WeightedSum[16]={0}, TotalCounts[16]={0}, Residual[16]={0};
		TF1 *KinFcn = new TF1("KinFcn",kin,0,180,1);
		KinFcn->SetParameter(0,KEBeam[KE]);
		TF1 *KinFcn2 = new TF1("KinFcn2",kin,angle[0],angle[16],1);
		KinFcn2->SetParameter(0,KEBeam[KE]);
		int NumberPoints = 0;
		//cout << "The value at 125 degrees is : "<< KinFcn->Eval(125.0,0,1) << endl;
		for (int i=0;i<16;i++){
			kin_value[i]=KinFcn->Eval(mid_ang[i],0,1);
			int Counts = 0;
			int NumCounts = 0;
			
			//for (int j=MinBin;j<MaxBin;j++){// GS
			for (int j=ExcitedStateCutBin[i][0];j<ExcitedStateCutBin[i][1];j++){// TES
				EnergyValue=h1[i]->GetXaxis()->GetBinCenter(j);
				Counts=h1[i]->GetBinContent(j);
				WeightedSum[i]=EnergyValue*Counts+WeightedSum[i];
				TotalCounts[i]=Counts+TotalCounts[i];
				for(int c=0; c<Counts; c++){
					Values[i][NumCounts]=EnergyValue-kin_value[i];
					NumCounts++;
				}
			}
			if (TotalCounts[i]==0){ //don't divide by zero
				Residual[i]=0;
			}else{
				Residual[i]=(WeightedSum[i]/TotalCounts[i])-kin_value[i];
			}
			for(int k=0; k<TotalCounts[i]; k++){
					gr1[KE]->SetPoint(NumberPoints,mid_ang[i],Values[i][k]);
					NumberPoints++;
			}
		}
		gr1[KE]->Fit(fit_lin,"QE","",angle[0],angle[16]);
		parameters << fit_lin->GetParameter(0) <<"\t"<< fit_lin->GetParError(0)<<"\t"<< fit_lin->GetParameter(1) <<"\t"<< fit_lin->GetParError(1) <<endl;
		c1->cd(KE+1);
		gr1[KE]->SetLineColor(kBlack);
		gr1[KE]->GetXaxis()->SetLimits(angle[0],angle[16]);
		gr1[KE]->SetMarkerStyle(21);
		string name2 = Form ( "Residual vs. Angle_%d",KE );
		gr1[KE]->SetTitle(name2.c_str());
		//gr1[KE]->SetTitle("Residual vs. Angle");
		gr1[KE]->Draw("AP");
		gr2->SetLineColor(kBlue);
		gr2->Draw("same");
		
		c2->cd(1);
		gr3->SetPoint(KE,KEBeam[KE],abs(fit_lin->GetParameter(0)));
		gr3->SetPointError(KE,0,fit_lin->GetParError(0));
		
		
		c2->cd(2);
		gr4->SetPoint(KE,KEBeam[KE],abs(fit_lin->GetParameter(1)));
		gr4->SetPointError(KE,0,fit_lin->GetParError(1));
		
	}
	c2->SetTitle("Third Excited State");
	c2->cd(1);
	gr3->SetLineColor(kBlack);
	gr3->SetTitle("Intercept");
	gr3->Draw("AP*");
	gr3->GetXaxis()->SetTitle("Energy (MeV)");
	gr3->GetYaxis()->SetTitle("Intercept");
		
	c2->cd(2);
	gr4->SetLineColor(kBlack);
	gr4->SetTitle("Slope");
	gr4->Draw("AP*");
	gr4->GetXaxis()->SetTitle("Energy (MeV)");
	gr4->GetYaxis()->SetTitle("Slope");
	
	/*fit_quad->SetParLimits(0,0,0.002);
	fit_quad->SetParLimits(1,-3.15e-05,-3.0e-05);
	fit_quad->SetParLimits(2,1.2e-07,1.22e-07);*/
	gr3->Fit(fit_quad,"QE","",KEBeam[0],KEBeam[en_length-1]);
	cout << "Minimum in Slope @ " << -fit_quad->GetParameter(1)/(2*fit_quad->GetParameter(2)) << " MeV" << endl;
	gr4->Fit(fit_quad,"QE","",KEBeam[0],KEBeam[en_length-1]);
	cout << "Minimum in Intercept @ " << -fit_quad->GetParameter(1)/(2*fit_quad->GetParameter(2)) << " MeV" << endl;

		
		/*TCanvas *c2 = new TCanvas ( "c2" ); //create a canvas
		c2->Divide ( 1,2 );

		c1->cd(1);
		gr1->SetLineColor(kBlack);
		gr1->GetXaxis()->SetLimits(angle[0],angle[16]);
		gr1->SetMarkerStyle(21);
		gr1->SetTitle("Residual vs. Angle");
		gr1->Draw("AP");*/
		//gr1->Fit(fit_lin,"E","",angle[0],angle[16]);
		//gr1->Fit("pol1");
		
		/*gr2->SetLineColor(kRed);
		gr2->Draw("same");
		c2->Update();
		
		c2->cd(2);
		h_Eloss->Draw("colz");
		KinFcn2->SetLineColor(kRed);
		KinFcn2->SetLineWidth(3);
		KinFcn2->Draw("same");*/
		
		//cout << fit_lin->GetParameter(0) <<"\t"<< fit_lin->GetParError(0)<<"\t"<< fit_lin->GetParameter(1) <<"\t"<< fit_lin->GetParError(1) <<endl;
}

