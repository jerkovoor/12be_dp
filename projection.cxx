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
	float kBeam = 114; //put the correct value; beam energy
	float mbeam = 12 * amu;  //mass of the beam (12Be or 12C) in MeV
	float mrecoil = 13 * amu;  //mass of the recoil (13Be or 13C) in MeV
	float mejec = 1 * amu; //mass of the proton
	float Qgs13C = 2.722; // Ground state Q value of 13C
	float Exenergy = 3.853;
	//float Exenergy = 0;
	float Q_C = Qgs13C-Exenergy;
	
		return pow(((TMath::Sqrt(mbeam*amu*kBeam)*TMath::Cos(TMath::DegToRad()*x[0])+TMath::Sqrt(mbeam*amu*kBeam*pow(TMath::Cos(TMath::DegToRad()*x[0]),2)+(mrecoil+amu)*(mrecoil*Q_C+(mrecoil-mbeam)*kBeam)))/(mrecoil+amu)),2); 
} 

void projection() {
	double angle[17]={0}, anglelow[16]={0}, anglehigh[16]={0}, mid_ang[16]={0}, kin_value[16]={0}, xaxis[2]={0}, errorbar[16]={0};
	int BinLow[16]={0}, BinHigh[16]={0};
	
	// Bin numbers 571 to 937 (energy from 2.85 to 4.68)were selected after inspecting all the 16 histograms
	int MinBin = 571;
	int MaxBin = 937;
	double GroundStateCutEn[16][2]={{2500,2800},{2500,2800},{2300,2600},{2200,2500},{2000,2300},{1900,2300},{1900,2300},{1900,2200},{1900,2200},{2100,2200},{1800,2200},{1900,2100},{1900,2100},{1600,1700},{1600,1800},{1500,1700}};
	double GroundStateCutBin[16][2];
	double ExcitedStateCutEn[16][2]={{1100,1600},{1100,1650},{1000,1500},{900,1500},{680,1400},{780,1400},{950,1500},{800,1400},{800,1400},{700,1300},{700,1300},{650,1300},{650,1200},{700,1200},{600,1100},{600,1100}};
	double ExcitedStateCutBin[16][2];
	
	TF1 *KinFcn = new TF1("KinFcn",kin,0,180,0);

	//cout << "The value at 125 degrees is : "<< KinFcn->Eval(125.0,0,0) << endl;
	
	TFile *f = new TFile("../Analysis/C_calib/C_pedestal_MicronDL_60sumT_calib_S04.root","READ");
	
	ofstream carb_cal;
	carb_cal.open("/home/jerome/12Be_exp/Analysis/C_calib/carbon_cal_TES_channel_SD.txt"); //open a .txt file to store the results of the fit; change the path and name acordingly
	carb_cal << "Strip\t"<< "ChannelNum\t"<< "SD\t" << "Energy" << endl;
	
	for (int i=0;i<17;i++){
			angle[i]=180-TMath::RadToDeg()*TMath::ATan((50+(15-i)*(4.94))/85);        
			//cout << angle[i] << endl;  
		}
	
	for(int sec=0;sec<8;sec++){
		TH2D *h_Eloss = (TH2D*)f->Get(Form("hYuAnPID%d",sec));
		//TH2D *h_noEloss = (TH2D*)f->Get("hYuAnPID1");

		
		TH1D *h1[16];
		TH1D *h2[16];
		
		double EnergyValue=0, energyValue=0, WeightedSum[16]={0}, weightedSum[16]={0}, TotalCounts[16]={0}, totalCounts[16]={0}, Average[16]={0}, average[16]={0}, SD[16]={0}, Residual[16]={0};
		
		for (int i=0;i<16;i++){
			anglelow[i]=angle[i]+0.01;
			anglehigh[i]=angle[i+1]-0.01;
			BinLow[i] = h_Eloss->GetXaxis()->FindBin(anglelow[i]);
			BinHigh[i] = h_Eloss->GetXaxis()->FindBin(anglehigh[i]);
			mid_ang[i]=(angle[i]+angle[i+1])/2;
			kin_value[i]=KinFcn->Eval(mid_ang[i],0,0);
			errorbar[i]=(angle[i+1]-angle[i])/2;
			h1[i] = h_Eloss->ProjectionY(Form("h1_%d",i),BinLow[i],BinHigh[i],"");
			//h2[i] = h_noEloss->ProjectionY(Form("h2_%d",i),BinLow[i],BinHigh[i],"");
			GroundStateCutBin[i][0]=h1[i]->GetXaxis()->FindBin(GroundStateCutEn[i][0]);
			GroundStateCutBin[i][1]=h1[i]->GetXaxis()->FindBin(GroundStateCutEn[i][1]);
			ExcitedStateCutBin[i][0]=h1[i]->GetXaxis()->FindBin(ExcitedStateCutEn[i][0]);
			ExcitedStateCutBin[i][1]=h1[i]->GetXaxis()->FindBin(ExcitedStateCutEn[i][1]);
			
			double Counts =0;
			for (int j=ExcitedStateCutBin[i][0];j<ExcitedStateCutBin[i][1];j++){
				EnergyValue=h1[i]->GetXaxis()->GetBinCenter(j);
				Counts=h1[i]->GetBinContent(j);
				WeightedSum[i]=EnergyValue*Counts+WeightedSum[i];
				TotalCounts[i]=Counts+TotalCounts[i];
			}
			
			double counts=0;
			double SD1[16]={0};
			for (int j=ExcitedStateCutBin[i][0];j<ExcitedStateCutBin[i][1];j++){
				energyValue=h1[i]->GetXaxis()->GetBinCenter(j);
				counts=h1[i]->GetBinContent(j);
				weightedSum[i]=energyValue*counts+weightedSum[i];
				totalCounts[i]=counts+totalCounts[i];
			}
			
			for (int j=ExcitedStateCutBin[i][0];j<ExcitedStateCutBin[i][1];j++){
				energyValue=h1[i]->GetXaxis()->GetBinCenter(j);
				counts=h1[i]->GetBinContent(j);
				if (totalCounts[i]==0){ //don't divide by zero
					SD1[i]=0;
				}else{
					average[i]=weightedSum[i]/totalCounts[i];
					SD1[i]=(pow(((average[i]-energyValue)*counts),2)/totalCounts[i])+SD1[i];
				}
			}
			
			SD[i]=sqrt(SD1[i]);
			
			if (TotalCounts[i]==0){ //don't divide by zero
				Residual[i]=0;
				carb_cal <<16*sec+i << "\t\t" << Average[i] << "\t\t\t" << 0 << "\t"<< kin_value[i] <<endl;
			}else{
				Average[i]=WeightedSum[i]/TotalCounts[i];
				Residual[i]=Average[i]-kin_value[i];
				carb_cal <<16*sec+i << "\t\t" << Average[i] << "\t\t\t" << SD[i]/sqrt(TotalCounts[i]) << "\t"<< kin_value[i] <<endl;
			}
			
		}
		
		int NumberPoints = 0;
		TGraphErrors *gr1 = new TGraphErrors();
		for (int i=0;i<16;i++){
			if (TotalCounts[i]==0){
				continue;
			}else{
				//cout << "Setting point " << NumberPoints << "\tValue : " << Residual[i] << "\tAngle : " << mid_ang[i] << " +/- " << errorbar[i] << endl;
				gr1->SetPoint(NumberPoints,mid_ang[i],Residual[i]);
				//gr1->SetPointError(NumberPoints,errorbar[i],0.1);
				NumberPoints++;
			}
		}
	}
	/*double x[2]={120,160};
	
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
	KinFcn2->Draw("same");*/
	
}

