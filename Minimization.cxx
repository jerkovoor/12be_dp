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

void Minimization() {
	
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
	double gainInitial=0.95;
	int gainLength=11;
	double Gain[gainLength];
	for (int i=0;i<gainLength;i++){
		Gain[i]=gainInitial+0.01*i;
	}
	int offsetLength=41;
	double Offset[offsetLength];
	
	double OffsetInitial=-0.2;
	for (int i=0;i<offsetLength;i++){
		Offset[i]=OffsetInitial+0.01*i;
	}
	
	TF1 *KinFcn = new TF1("KinFcn",kin,0,180,1);
	KinFcn->SetParameter(0,Qgs13C);
	TF1 *KinFcn2 = new TF1("KinFcn2",kin,0,180,1);
	KinFcn2->SetParameter(0,QTES13C);
	
	ofstream ChiSq;
	ChiSq.open("/home/jerome/12Be_exp/Analysis/CarbonGain/ChiSquare.txt");//GS
	ChiSq << "Gain\t"<< "Offset\t" << "Chi Square" << endl;

	//cout << "The value at 125 degrees is : "<< KinFcn->Eval(125.0,0,0) << endl;
	
	TFile *f = new TFile("../Analysis/CarbonGain/C_pedestal_TRIUMF_DL_CarbonGain_Yu_binTest.root","READ");
	//TFile *f = new TFile("../Analysis/CarbonGain/C_pedestal_Micron_DL_Yu_shifted.root","READ");
	TH2D *h_YuAnPID = (TH2D*)f->Get("hYuAnPID_20");
	for (int i=0;i<17;i++){
			angle[i]=180-TMath::RadToDeg()*TMath::ATan((50+(16-i)*(4.94))/85);        
			//cout << angle[i] << endl;  
	}
	
	
	TGraph* gr1 = new TGraph();
	TGraph* gr2 = new TGraph();
	TH1D *h1[16];
	TF1* f1 = new TF1 ( "f1","[0]+[1]*x^(-1)+[2]*x^(-2)",0,180); 
	
	
	
	for (int GainNumber=0;GainNumber<gainLength;GainNumber++){
		for (int OffsetNumber=0;OffsetNumber<offsetLength;OffsetNumber++){
			int NumberPointsGS = 0, NumberPointsTES = 0;
			double EnergyValueGS=0, EnergyValueTES=0, ModEnergyValueTES, ModEnergyValueGS, WeightedSumGS[16]={0}, WeightedSumTES[16]={0}, TotalCountsGS[16]={0}, TotalCountsTES[16]={0}, AverageGS[16]={0}, AverageTES[16]={0}, DiffGS[16]={0}, DiffTES[16]={0}, DiffSqGS[16]={0}, DiffSqTES[16]={0}, shiftGS[16]={0}, shiftTES[16]={0}, ChiSquare[16]={0};
			double TotChiSquare=0;
			for (int i=0;i<16;i++){
				anglelow[i]=angle[i]+0.01;
				anglehigh[i]=angle[i+1]-0.01;
				BinLow[i] = h_YuAnPID->GetXaxis()->FindBin(anglelow[i]);
				BinHigh[i] = h_YuAnPID->GetXaxis()->FindBin(anglehigh[i]);
				mid_ang[i]=(angle[i]+angle[i+1])/2;
				kin_valueGS[i]=KinFcn->Eval(mid_ang[i],0,0);//GS
				kin_valueTES[i]=KinFcn2->Eval(mid_ang[i],0,0);// TES
				errorbar[i]=(angle[i+1]-angle[i])/2;
				h1[i] = h_YuAnPID->ProjectionY(Form("h1_%d",i),BinLow[i],BinHigh[i],"");
				GroundStateCutBin1=h1[i]->GetXaxis()->FindBin(GroundStateCutEn1);
				GroundStateCutBin2=h1[i]->GetXaxis()->FindBin(GroundStateCutEn2);
				ExcitedStateCutBin[i][0]=h1[i]->GetXaxis()->FindBin(ExcitedStateCutEn[i][0]);
				ExcitedStateCutBin[i][1]=h1[i]->GetXaxis()->FindBin(ExcitedStateCutEn[i][1]);
				
				double CountsGS =0;
				for (int j=GroundStateCutBin1;j<GroundStateCutBin2;j++){//GS
					EnergyValueGS=h1[i]->GetXaxis()->GetBinCenter(j);
					ModEnergyValueGS=Gain[GainNumber]*EnergyValueGS+Offset[OffsetNumber];
					CountsGS=h1[i]->GetBinContent(j);
					DiffSqGS[i]=pow((ModEnergyValueGS-kin_valueGS[i]),2)*CountsGS+DiffSqGS[i];
					WeightedSumGS[i]=ModEnergyValueGS*CountsGS+WeightedSumGS[i];
					TotalCountsGS[i]=CountsGS+TotalCountsGS[i];
				}
				
				double CountsTES =0;
				for (int j=ExcitedStateCutBin[i][0];j<ExcitedStateCutBin[i][1];j++){//TES
					EnergyValueTES=h1[i]->GetXaxis()->GetBinCenter(j);
					ModEnergyValueTES=Gain[GainNumber]*EnergyValueTES+Offset[OffsetNumber];
					CountsTES=h1[i]->GetBinContent(j);
					DiffSqTES[i]=pow((ModEnergyValueTES-kin_valueTES[i]),2)*CountsTES+DiffSqTES[i];
					WeightedSumTES[i]=ModEnergyValueTES*CountsTES+WeightedSumTES[i];
					TotalCountsTES[i]=CountsTES+TotalCountsTES[i];
				}
				
				AverageGS[i]=WeightedSumGS[i]/TotalCountsGS[i];
				AverageTES[i]=WeightedSumTES[i]/TotalCountsTES[i];
				DiffGS[i]=AverageGS[i]-kin_valueGS[i];
				DiffTES[i]=AverageTES[i]-kin_valueTES[i];
				if (TotalCountsTES[i]==0){
					shiftTES[i]=0;
				}else{
					shiftTES[i]=DiffSqTES[i]/TotalCountsTES[i];
				}
				if (TotalCountsGS[i]==0){
					shiftGS[i]=0;
				}else{
					shiftGS[i]=DiffSqGS[i]/TotalCountsGS[i];
				}
				ChiSquare[i]=shiftGS[i]+shiftTES[i];
				TotChiSquare=ChiSquare[i]+TotChiSquare;
				//cout << shiftTES[i] << endl;
			}
			ChiSq <<Gain[GainNumber] << "\t\t" << Offset[OffsetNumber] << "\t\t" << TotChiSquare <<endl;
		}
	}
}