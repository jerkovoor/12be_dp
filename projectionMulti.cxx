//Calculates and plots the residual vs. the angle.
//Change the script accordingly to look at the ground state and the third excited state

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
	float Exenergy = 0;
	//float Exenergy = 3.853; //Third excited state energy
	float Q_C = Qgs13C-Exenergy;
	
	return pow(((TMath::Sqrt(mbeam*amu*kBeam)*TMath::Cos(TMath::DegToRad()*x[0])+TMath::Sqrt(mbeam*amu*kBeam*pow(TMath::Cos(TMath::DegToRad()*x[0]),2)+(mrecoil+amu)*(mrecoil*Q_C+(mrecoil-mbeam)*kBeam)))/(mrecoil+amu)),2); 
} 

void projectionMulti() {
	
	double angle[17]={0}, anglelow[16]={0}, anglehigh[16]={0}, mid_ang[16]={0}, kin_value[16]={0}, WeightedSum[16]={0}, TotalCounts[16]={0}, Residual[16]={0}, xaxis[2]={0}, errorbar[16]={0};
	int BinLow[16]={0}, BinHigh[16]={0};
	double ExcitedStateCutEn[16][2]={{2,3},{2,3},{2,3},{1.8,2.8},{1.6,2.6},{1.6,2.6},{1.6,2.6},{1.6,2.6},{1.6,2.6},{1.4,2.4},{1.4,2.4},{1.2,2.4},{1.2,2.2},{1.3,2.1},{1.2,2.2},{1.2,2}};//TES
	double ExcitedStateCutBin[16][2];
	double Values[16][40];
	
	// Bin numbers 571 to 937 (energy from 2.85 to 4.68) were selected after inspecting all the 16 histograms for the GS
	int MinBin = 578;
	int MaxBin = 937;
	
	/*TCanvas *c1 = new TCanvas ( "c1" ); //create a canvas
	c1->Divide ( 4,4 );*/
	

	
		TF1 *KinFcn = new TF1("KinFcn",kin,0,180,0);


		//cout << "The value at 125 degrees is : "<< KinFcn->Eval(125.0,0,0) << endl;
		
		TFile *f = new TFile("../Analysis/C_calib/C_nopedestal_CCalib.root","READ");
		TH2D *h_Eloss = (TH2D*)f->Get("hYuAnPID");
		TH2D *h_noEloss = (TH2D*)f->Get("hYuAnPID1");

		//h->Draw("colz");
		for (int i=0;i<17;i++){
			angle[i]=180-TMath::RadToDeg()*TMath::ATan((50+(15-i)*(4.94))/85);        
			//cout << angle[i] << endl;  
		}
		
		/*TCanvas *c1 = new TCanvas ( "c1" ); //create a canvas
		c1->Divide ( 4,4 );*/ //divide the Canvas in 16
		
		TH1D *h1[16];
		TH1D *h2[16];
		
		double EnergyValue=0;
		
		for (int i=0;i<16;i++){
			anglelow[i]=angle[i]+0.01;
			anglehigh[i]=angle[i+1]-0.01;
			BinLow[i] = h_Eloss->GetXaxis()->FindBin(anglelow[i]);
			BinHigh[i] = h_Eloss->GetXaxis()->FindBin(anglehigh[i]);
			mid_ang[i]=(angle[i]+angle[i+1])/2;
			kin_value[i]=KinFcn->Eval(mid_ang[i],0,1);
			errorbar[i]=(angle[i+1]-angle[i])/2;
			
			
			/*cout << mid_ang[i] << endl;
			cout << kin_value[i] << endl;
			
			cout << "Angle = " << anglelow[i] << "    Bin = " << BinLow[i] << endl;
			cout << "Angle = " << anglehigh[i] << "    Bin = " << BinHigh[i] << endl;*/
			
			
			h1[i] = h_Eloss->ProjectionY(Form("h1_%d",i),BinLow[i],BinHigh[i],"");
			h2[i] = h_noEloss->ProjectionY(Form("h2_%d",i),BinLow[i],BinHigh[i],"");
			ExcitedStateCutBin[i][0]=h1[i]->GetXaxis()->FindBin(ExcitedStateCutEn[i][0]);
			ExcitedStateCutBin[i][1]=h1[i]->GetXaxis()->FindBin(ExcitedStateCutEn[i][1]);
			
			int Counts = 0;
			int NumCounts = 0;
			for (int j=MinBin;j<MaxBin;j++){// GS
			//for (int j=ExcitedStateCutBin[i][0];j<ExcitedStateCutBin[i][1];j++){// TES
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
			//cout << Residual[i] << endl;
			
			/*c1->cd(i+1);
			
			h1[i]->SetLineColor(kBlue);
			h1[i]->Draw();
			h2[i]->SetLineColor(kRed);
			h2[i]->Draw("same");*/	
			
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
		
		
		double x[2]={120,160};
		
		TF1 *fit_lin = new TF1 ( "fit_lin","pol1",0,200 );
		//fit_lin->SetParameters(-0.1,0);
		
		TGraph *gr2 = new TGraph(2,x,xaxis);
		
		TF1 *KinFcn2 = new TF1("KinFcn2",kin,angle[0],angle[16],0);
		
		
		TCanvas *c2 = new TCanvas ( "c2" ); //create a canvas
		c2->Divide ( 1,2 );

		c2->cd(1);
		gr1->SetLineColor(kBlack);
		gr1->GetXaxis()->SetLimits(angle[0],angle[16]);
		gr1->SetMarkerStyle(21);
		gr1->SetTitle("Residual vs. Angle");
		gr1->Draw("AP");
		gr1->Fit(fit_lin,"E","",angle[0],angle[16]);
		gr1->Fit("pol1");
		
		gr2->SetLineColor(kRed);
		gr2->Draw("same");
		c2->Update();
		
		c2->cd(2);
		h_Eloss->Draw("colz");
		KinFcn2->SetLineColor(kRed);
		KinFcn2->SetLineWidth(3);
		KinFcn2->Draw("same");
		
		cout << fit_lin->GetParameter(0) <<"\t"<< fit_lin->GetParError(0)<<"\t"<< fit_lin->GetParameter(1) <<"\t"<< fit_lin->GetParError(1) <<endl;
	
}

