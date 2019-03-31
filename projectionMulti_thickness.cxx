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



double kin(double *x, double *par){ 
	
	float amu = 931.5; // atomic mass unit in MeV 
	float massEjec = 938.28;//mass of the proton in MeV/c2 
	float kBeam = 110; //put the correct value; beam energy
	float mbeam = 12 * amu;  //mass of the beam (12Be or 12C) in MeV
	float mrecoil = 13 * amu;  //mass of the recoil (13Be or 13C) in MeV
	float mejec = 1 * amu; //mass of the proton
	float Qgs13C = 2.722; // Ground state Q value of 13C
	float Exenergy = 3.853; //Third excited state energy
	float Q_C = Qgs13C-Exenergy;
	
	return pow(((TMath::Sqrt(mbeam*amu*kBeam)*TMath::Cos(TMath::DegToRad()*x[0])+TMath::Sqrt(mbeam*amu*kBeam*pow(TMath::Cos(TMath::DegToRad()*x[0]),2)+(mrecoil+amu)*(mrecoil*Q_C+(mrecoil-mbeam)*kBeam)))/(mrecoil+amu)),2); 
} 

void projectionMulti_thickness() {
	
	double ExcitedStateCutEn[16][2]={{2,3},{2,3},{2,3},{1.8,2.8},{1.6,2.6},{1.6,2.6},{1.6,2.6},{1.6,2.6},{1.6,2.6},{1.4,2.4},{1.4,2.4},{1.2,2.4},{1.2,2.2},{1.3,2.1},{1.2,2.2},{1.2,2}};//TES
	double ExcitedStateCutBin[16][2], xaxis[2]={0};
	
	// Bin numbers 571 to 937 (energy from 2.85 to 4.68) were selected after inspecting all the 16 histograms for the GS
	int MinBin = 571;
	int MaxBin = 937;
	TH1D *h1[16];
	TH1D *h2[16];
	TCanvas *c1 = new TCanvas ( "c1" ); //create a canvas
	c1->Divide ( 4,4 );
	TCanvas *c2 = new TCanvas ( "c2" ); //create a canvas
	c2->Divide ( 1,2 );
	double x[2]={120,160};
	TGraph *gr2 = new TGraph(2,x,xaxis);
	TGraphErrors *gr3 = new TGraphErrors();
	TGraphErrors *gr4 = new TGraphErrors();

	
	TF1 *KinFcn = new TF1("KinFcn",kin,0,180,0);
	TF1 *fit_lin = new TF1 ( "fit_lin","pol1",0,200 );
	TF1 *fit_quad = new TF1 ("fit_quad","[0]+[1]*x+[2]*pow(x,2)",0,200);
	int th=67;
	int iter=0;
	
	while (th<=97){
		double angle[17]={0}, anglelow[16]={0}, anglehigh[16]={0}, mid_ang[16]={0}, kin_value[16]={0}, WeightedSum[16]={0}, TotalCounts[16]={0}, Residual[16]={0}, errorbar[16]={0};
		int BinLow[16]={0}, BinHigh[16]={0};
		double Values[16][40]={0};
		string f_name=Form("../Analysis/with_pedestal/C_pedestal_MicronDL_60sumT_%imm.root",th);
		TFile *f = new TFile(f_name.c_str(),"READ");
		TH2D *h_Eloss = (TH2D*)f->Get("hYuAnPID");
		TH2D *h_noEloss = (TH2D*)f->Get("hYuAnPID1");
		for (int i=0;i<17;i++){
			angle[i]=180-TMath::RadToDeg()*TMath::ATan((50+(15-i)*(4.94))/th);        
			//cout << angle[i] << endl;  
		}
		
		double EnergyValue=0;
		
		for (int i=0;i<16;i++){
			anglelow[i]=angle[i]+0.01;
			anglehigh[i]=angle[i+1]-0.01;
			BinLow[i] = h_Eloss->GetXaxis()->FindBin(anglelow[i]);
			BinHigh[i] = h_Eloss->GetXaxis()->FindBin(anglehigh[i]);
			mid_ang[i]=(angle[i]+angle[i+1])/2;
			kin_value[i]=KinFcn->Eval(mid_ang[i],0,1);
			errorbar[i]=(angle[i+1]-angle[i])/2;
			h1[i] = h_Eloss->ProjectionY(Form("h1_%d",i),BinLow[i],BinHigh[i],"");
			h2[i] = h_noEloss->ProjectionY(Form("h2_%d",i),BinLow[i],BinHigh[i],"");
			ExcitedStateCutBin[i][0]=h1[i]->GetXaxis()->FindBin(ExcitedStateCutEn[i][0]);
			//ExcitedStateCutBin[i][1]=h1[i]->GetXaxis()->FindBin(ExcitedStateCutEn[i][1]);
			ExcitedStateCutBin[i][1]=h1[i]->GetXaxis()->FindBin((ExcitedStateCutEn[i][0]+ExcitedStateCutEn[i][1])/2);
			
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
		}
		
		int NumberPoints = 0;
		TGraphErrors *gr1 = new TGraphErrors();
		for (int i=0;i<16;i++){
	/*
			if (TotalCounts[i]==0){
				continue;
			}else{
				cout << "Setting point " << NumberPoints << "\tValue : " << Residual[i] << "\tAngle : " << mid_ang[i] << " +/- " << errorbar[i] << endl;
				gr1->SetPoint(NumberPoints,mid_ang[i],Residual[i]);
				//gr1->SetPointError(NumberPoints,errorbar[i],0.1);
				NumberPoints++;
			}
	*/
			for(int k=0; k<TotalCounts[i]; k++){
				gr1->SetPoint(NumberPoints,mid_ang[i],Values[i][k]);
				NumberPoints++;	
			}	
			
		}
		
		/*TGraph *gr1 = new TGraph();
		for (int i=0;i<16;i++){
			if (TotalCounts[i]==0){
				continue;
			}else{
				gr1->SetPoint(i,mid_ang[i],Residual[i]);
			}
		}*/
		
		c1->cd(iter+1);
		gr1->SetLineColor(kBlack);
		gr1->GetXaxis()->SetLimits(angle[0],angle[16]);
		gr1->SetMarkerStyle(21);
		gr1->SetTitle("Residual vs. Angle");
		gr1->Draw("AP");
		gr1->Fit(fit_lin,"QE","",angle[0],angle[16]);
		gr2->SetLineColor(kRed);
		gr2->Draw("same");
		
		c2->cd(1);
		gr3->SetPoint(iter,th,abs(fit_lin->GetParameter(0)));
		gr3->SetPointError(iter,0,fit_lin->GetParError(0));
		
		
		c2->cd(2);
		gr4->SetPoint(iter,th,abs(fit_lin->GetParameter(1)));
		gr4->SetPointError(iter,0,fit_lin->GetParError(1));
		th=th+2;
		iter++;
		//cout << "Next" << endl;
	}
	
	c2->SetTitle("Third Excited State");
	c2->cd(1);
	gr3->SetLineColor(kBlack);
	gr3->SetTitle("Intercept");
	gr3->Draw("AP*");
	gr3->GetXaxis()->SetTitle("Target Distance (mm)");
	gr3->GetYaxis()->SetTitle("Intercept");
		
	c2->cd(2);
	gr4->SetLineColor(kBlack);
	gr4->SetTitle("Slope");
	gr4->Draw("AP*");
	gr4->GetXaxis()->SetTitle("Target Distance (mm)");
	gr4->GetYaxis()->SetTitle("Slope");
	
	/*fit_quad->SetParLimits(0,0,0.002);
	fit_quad->SetParLimits(1,-3.15e-05,-3.0e-05);
	fit_quad->SetParLimits(2,1.2e-07,1.22e-07);*/
	gr3->Fit(fit_quad,"QE","",67,97);
	cout << "Minimum in Slope @ " << -fit_quad->GetParameter(1)/(2*fit_quad->GetParameter(2)) << " mm" << endl;
	gr4->Fit(fit_quad,"QE","",67,97);
	cout << "Minimum in Intercept @ " << -fit_quad->GetParameter(1)/(2*fit_quad->GetParameter(2)) << " mm" << endl;
	
}

