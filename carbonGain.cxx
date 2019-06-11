using namespace std;

#include <TFile.h> 
#include <TMath.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>



double kin(double *x, double *Q_C){ 
	
	float amu = 931.5; // atomic mass unit in MeV 
	float massEjec = 938.28;//mass of the proton in MeV/c2 
	float kBeam = 110; //put the correct value; beam energy
	float mbeam = 12 * amu;  //mass of the beam (12Be or 12C) in MeV
	float mrecoil = 13 * amu;  //mass of the recoil (13Be or 13C) in MeV
	float mejec = 1 * amu; //mass of the proton
	
		return pow(((TMath::Sqrt(mbeam*amu*kBeam)*TMath::Cos(TMath::DegToRad()*x[0])+TMath::Sqrt(mbeam*amu*kBeam*pow(TMath::Cos(TMath::DegToRad()*x[0]),2)+(mrecoil+amu)*(mrecoil*Q_C[0]+(mrecoil-mbeam)*kBeam)))/(mrecoil+amu)),2); 
} 

void carbonGain() {
	
	double angle[17]={0}, anglelow[16]={0}, anglehigh[16]={0}, mid_ang[16]={0}, kin_valueGS[16]={0}, kin_valueTES[16]={0}, xaxis[2]={0}, errorbar[16]={0};
	int BinLow[16]={0}, BinHigh[16]={0};
	// Bin numbers 571 to 937 (energies from 2.85 to 4.68) were selected after inspecting all the 16 histograms for the ground state
	double GroundStateCutEn1 = 2.85, GroundStateCutEn2 = 4.68;//with pedestal
	//int GroundStateCutBin = 583; //without pedestal (583 for alpha calibration and 578 for Carbon secondary calibration)
	int GroundStateCutBin1, GroundStateCutBin2;
	double ExcitedStateCutEn[16][2]={{2,3},{2,3},{2,3},{1.8,2.8},{1.6,2.6},{1.6,2.6},{1.6,2.6},{1.6,2.6},{1.6,2.6},{1.4,2.4},{1.4,2.4},{1.2,2.4},{1.2,2.2},{1.3,2.1},{1.2,2.2},{1.2,2}};//for third excited state
	double ExcitedStateCutBin[16][2];
	float Qgs13C = 2.722; // Ground state Q value of 13C
	float Exenergy = 3.853;//TES
	float QTES13C = Qgs13C-Exenergy;
	
	TF1 *KinFcn = new TF1("KinFcn",kin,0,180,1);
	KinFcn->SetParameter(0,Qgs13C);
	TF1 *KinFcn2 = new TF1("KinFcn2",kin,0,180,1);
	KinFcn2->SetParameter(0,QTES13C);

	//cout << "The value at 125 degrees is : "<< KinFcn->Eval(125.0,0,0) << endl;
	
	//TFile *f = new TFile("../Analysis/CarbonGain/C_pedestal_Micron_DL_Yu.root","READ");
	TFile *f = new TFile("../Analysis/CarbonGain/C_pedestal_Micron_DL_Yu_wholeShifted.root","READ");
	TH2D *h_YuAnPID = (TH2D*)f->Get("hYuAnPID");
	
	/*ofstream carb_gainGS;
	carb_gainGS.open("/home/jerome/12Be_exp/Analysis/CarbonGain/carbon_gain_GS_shift.txt");//GS
	carb_gainGS << "Strip\t"<< "DiffSq\t" << "Diff" << endl;*/
	
	/*ofstream carb_gainTES;
	carb_gainTES.open("/home/jerome/12Be_exp/Analysis/CarbonGain/carbon_gain_TES_shift.txt"); //TES
	carb_gainTES << "Strip\t"<< "DiffSq\t" << "Diff" << endl;*/
	
	for (int i=0;i<17;i++){
			angle[i]=180-TMath::RadToDeg()*TMath::ATan((50+(15-i)*(4.94))/85);        
			//cout << angle[i] << endl;  
	}
	
	
	TCanvas *c1 = new TCanvas ( "c1" );
	c1->Divide(2,4);
	TCanvas *c2 = new TCanvas ( "c2" );
	c2->Divide(2,4);
	
	for(int sec=0;sec<8;sec++){
		TH2D *h_Eloss = (TH2D*)f->Get(Form("hYuAnPIDSec_%d",sec));
		
		TH1D *h1[16];

		double EnergyValueGS=0, EnergyValueTES=0, WeightedSumGS[16]={0}, WeightedSumTES[16]={0}, TotalCountsGS[16]={0}, TotalCountsTES[16]={0}, AverageGS[16]={0}, AverageTES[16]={0}, DiffGS[16]={0}, DiffTES[16]={0}, DiffSqGS[16]={0}, DiffSqTES[16]={0}, shiftGS[16]={0}, shiftTES[16]={0};
		
		for (int i=0;i<16;i++){
			anglelow[i]=angle[i]+0.01;
			anglehigh[i]=angle[i+1]-0.01;
			BinLow[i] = h_Eloss->GetXaxis()->FindBin(anglelow[i]);
			BinHigh[i] = h_Eloss->GetXaxis()->FindBin(anglehigh[i]);
			mid_ang[i]=(angle[i]+angle[i+1])/2;
			kin_valueGS[i]=KinFcn->Eval(mid_ang[i],0,0);//GS, TES
			kin_valueTES[i]=KinFcn2->Eval(mid_ang[i],0,0);
			errorbar[i]=(angle[i+1]-angle[i])/2;
			h1[i] = h_Eloss->ProjectionY(Form("h1_%d",i),BinLow[i],BinHigh[i],"");
			GroundStateCutBin1=h1[i]->GetXaxis()->FindBin(GroundStateCutEn1);
			GroundStateCutBin2=h1[i]->GetXaxis()->FindBin(GroundStateCutEn2);
			ExcitedStateCutBin[i][0]=h1[i]->GetXaxis()->FindBin(ExcitedStateCutEn[i][0]);
			ExcitedStateCutBin[i][1]=h1[i]->GetXaxis()->FindBin(ExcitedStateCutEn[i][1]);
			
			double CountsGS =0;
			for (int j=GroundStateCutBin1;j<GroundStateCutBin2;j++){//GS
				EnergyValueGS=h1[i]->GetXaxis()->GetBinCenter(j);
				CountsGS=h1[i]->GetBinContent(j);
				DiffSqGS[i]=pow((EnergyValueGS-kin_valueGS[i]),2)*CountsGS+DiffSqGS[i];
				WeightedSumGS[i]=EnergyValueGS*CountsGS+WeightedSumGS[i];
				TotalCountsGS[i]=CountsGS+TotalCountsGS[i];
			}
			
			double CountsTES =0;
			for (int j=ExcitedStateCutBin[i][0];j<ExcitedStateCutBin[i][1];j++){//TES
				EnergyValueTES=h1[i]->GetXaxis()->GetBinCenter(j);
				CountsTES=h1[i]->GetBinContent(j);
				DiffSqTES[i]=pow((EnergyValueTES-kin_valueTES[i]),2)*CountsTES+DiffSqTES[i];
				WeightedSumTES[i]=EnergyValueTES*CountsTES+WeightedSumTES[i];
				TotalCountsTES[i]=CountsTES+TotalCountsTES[i];
			}
			
			AverageGS[i]=WeightedSumGS[i]/TotalCountsGS[i];
			AverageTES[i]=WeightedSumTES[i]/TotalCountsTES[i];
			DiffGS[i]=AverageGS[i]-kin_valueGS[i];
			DiffTES[i]=AverageTES[i]-kin_valueTES[i];
			shiftGS[i]=DiffSqGS[i]/TotalCountsGS[i];
			shiftTES[i]=DiffSqTES[i]/TotalCountsTES[i];

			/*if (TotalCountsGS[i]==0){ //don't divide by zero
				carb_gainGS <<16*sec+i << "\t\t" << 0 <<endl;
			}else{
				//AverageGS[i]=WeightedSumGS[i]/TotalCountsGS[i];
				//ResidualGS[i]=AverageGS[i]-kin_valueGS[i];
				carb_gainGS <<16*sec+i << "\t\t" << shiftGS[i] << "\t\t" << DiffGS[i] <<endl;
			}
			
			if (TotalCountsTES[i]==0){ //don't divide by zero
				carb_gainTES <<16*sec+i << "\t\t" << 0 <<endl;
			}else{
				//AverageTES[i]=WeightedSumTES[i]/TotalCountsTES[i];
				//ResidualTES[i]=AverageTES[i]-kin_valueTES[i];
				carb_gainTES <<16*sec+i << "\t\t" << shiftTES[i] << "\t\t" << DiffTES[i] <<endl;
			}*/
		}
		
		int NumberPointsGS = 0;
		TGraph* gr1 = new TGraph();
		c1->cd(sec+1);
		for (int i=0;i<16;i++){
			if (TotalCountsGS[i]==0){
				continue;
			}else{
				gr1->SetPoint(NumberPointsGS,mid_ang[i],AverageGS[i]);
				NumberPointsGS++;
			}
		}
		gr1->Draw ( "AP*" );
		gr1->GetYaxis()->SetRangeUser(2,5);//GS
		gr1->GetXaxis()->SetRangeUser(120,160);
		KinFcn->Draw("same");//GS
		
		int NumberPointsTES = 0;
		TGraph* gr2 = new TGraph();
		c2->cd(sec+1);
		for (int i=0;i<16;i++){
			if (TotalCountsTES[i]==0){
				continue;
			}else{
				gr2->SetPoint(NumberPointsTES,mid_ang[i],AverageTES[i]);
				NumberPointsTES++;
			}
		}
		gr2->Draw ( "AP*" );
		gr2->GetYaxis()->SetRangeUser(1,3);//TES
		gr2->GetXaxis()->SetRangeUser(120,160);
		KinFcn2->Draw("same");//TES
	}
	TCanvas *c3 = new TCanvas ( "c3" );
	TGraph* gr3 = new TGraph();
	h_YuAnPID->Draw("colz");
	KinFcn->Draw("same");
	KinFcn->SetLineColor(kGreen);
	KinFcn->SetLineWidth(3);
	KinFcn2->Draw("same");
	KinFcn2->SetLineColor(kRed);
	KinFcn2->SetLineWidth(3);
}