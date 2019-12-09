using namespace std;


#include <TMath.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include "EnergyLoss.h"
#include <TSpectrum.h>
#include <string>


double gauss1(double *x, double *p) {
    return p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2]));
}

void targetThicknessRunByRunBe(){
    
    
    TFile *g = new TFile("/home/jerome/12Be_exp/Analysis/TargetThickness/FromBeBeam/Be_target/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2_target_thickness_Sd1rAlphaCal_Sd2rNewInBeamCal_SingleRun.root","READ");
    //TFile *g1 = new TFile("Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2_target_thickness.root
    TH1D *hCSd1r_0_1_2_En[140];// number of runs = 140
    TH1D *hCSd2r_0_1_2_En[140];
    
    TFile *g1 = new TFile("/home/jerome/12Be_exp/Analysis/TargetThickness/FromBeBeam/Be_target/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2_target_thickness_Sd1rAlphaCal_Sd2rNewInBeamCal_SingleRun_withFit.root","RECREATE");
    
    TF1* fit_gauss1 = new TF1("fit_gauss1",gauss1,0,200,3);
    
    
    fit_gauss1->SetParLimits (2,0.1,20);// Standard deviation

    EnergyLoss* Be_Al = new EnergyLoss("Be_in_Al.dat");
    EnergyLoss* Be_SiO2 = new EnergyLoss("Be_in_SiO2.dat");
	EnergyLoss* Be_B = new EnergyLoss("Be_in_B.dat");
    EnergyLoss* Be_Si = new EnergyLoss("Be_in_Si.dat");
    EnergyLoss* Be_P = new EnergyLoss("Be_in_P.dat");
    EnergyLoss* Be_Ag = new EnergyLoss("Be_in_Ag.dat");
    EnergyLoss* Be_D2 = new EnergyLoss("Be_in_D2_density0_201.dat");
    
    /*Be_Al->AddBackHigh(15.);
    Be_SiO2->AddBackHigh(15.);
    Be_B->AddBackHigh(15.);
    */
    
    double dAg=4.6405042/1000.;
    
    //Deadlayer D1
    double dSdAl1=1.5/1000.;
    double dSdSiO2=3.5/1000.;
    double dSdAl2=0.3/1000.;
    double dSdB=0.5/1000.;
    
    //Active Si layer
    double dSd1Si=61./1000.;
    double dSd2Si=493./1000.;
    
    //Deadlayer D2
    double dSdP=0.5/1000.;
    double dSdAl3=0.3/1000.;
    

    
    double E_Be_after_Ag=108.57;
    
    int length = 140; //Number of runs
    
    double ESd1Si_Target[length], ESd2Si_Target[length], E_Be_Sd2_P_Target[length], E_Be_Sd2_Al3_Target[length], E_Be_Sd1_Al3_Target[length], E_Be_Sd1_P_Target[length], Sd1rE_Target[length], E_Be_B_Target[length], E_Be_Al2_Target[length], E_Be_SiO2_Target[length], E_Be_Al1_Target[length], E_Be_before_Ag_Target[length], runNumber[length], errorSd1r[length], errorSd2r[length], errorTotal[length], thickness[length]; 
    
    //Without Target
      
    double energy = 77.44;//77.05;
        
    double ESd1Si_noTarget=8.987;//24.54;
    double ESd2Si_noTarget=101.79;//77.05;
        
    //cout << "sd2 Energy: " << energy << endl;
        
    double E_Be_Sd2_P_noTarget=Be_P->AddBack(ESd2Si_noTarget,dSdP);
    energy = Be_P->AddBack(energy, dSdP);
        
    //cout << "energy before P dead layer (sd2): " << energy << endl;
        
    double E_Be_Sd2_Al3_noTarget=Be_Al->AddBack(E_Be_Sd2_P_noTarget,dSdAl3);
    energy = Be_Al->AddBack(energy, dSdAl3);
        
    //cout << "energy before Al dead layer (sd2): " << energy << endl;
        
    double E_Be_Sd1_Al3_noTarget=Be_Al->AddBack(E_Be_Sd2_Al3_noTarget,dSdAl3);
    energy = Be_Al->AddBack(energy, dSdAl3);
        
    //cout << "energy before Al dead layer (sd1): " << energy << endl;
        
    double E_Be_Sd1_P_noTarget=Be_P->AddBack(E_Be_Sd1_Al3_noTarget,dSdP);
    energy = Be_P->AddBack(energy, dSdP);
        
    //cout << "energy before P dead layer (sd1): " << energy << endl;
        
    double Sd1rE_noTarget=ESd1Si_noTarget+E_Be_Sd1_P_noTarget;
    energy = ESd1Si_noTarget + energy;
        
    //cout << "sd1 energy: " << ESd1Si_noTarget << endl;
    //cout << "sd1 energy + sd2 energy with dead layers: " << energy << endl;
        
    //cout << "sd1 deposited from SRIM: " << 101.8 - Be_Si->CalcRemainder(101.8, dSd1Si) << endl;
        
    double E_Be_B_noTarget = Be_B->AddBack(Sd1rE_noTarget,dSdB);
    double E_Be_Al2_noTarget = Be_Al->AddBack(E_Be_B_noTarget,dSdAl2);
    double E_Be_SiO2_noTarget = Be_SiO2->AddBack(E_Be_Al2_noTarget,dSdSiO2);
    double E_Be_Al1_noTarget = Be_Al->AddBack(E_Be_SiO2_noTarget,dSdAl1);
    double E_Be_before_Ag_noTarget = Be_Ag->AddBack(E_Be_Al1_noTarget,dAg); //Initial Energy
        
    //cout << E_Be_Sd2_Al3_noTarget << endl;
        
    cout << E_Be_before_Ag_noTarget << endl;
    
    //cout << ESi_noTarget << endl;
    //cout << E_Be_Si << endl;
    TGraphErrors *gr=new TGraphErrors();
    TGraphErrors *gr1=new TGraphErrors();
    TGraphErrors *gr2=new TGraphErrors();
    TGraphErrors *gr3=new TGraphErrors();
    g1->cd();
    
    int npeaks=1;
    TSpectrum* s = new TSpectrum ( 2*npeaks );
    int nfound=0;
    double* txp = s->GetPositionX();
    
    double Sd2rxp[3];
    double Sd2ra[length];
    
    for ( int j=0; j<length; j++ ){
        hCSd1r_0_1_2_En[j] = (TH1D*)g->Get(Form("hCSd1r_0_1_2_En_%i",j));
        hCSd2r_0_1_2_En[j] = (TH1D*)g->Get(Form("hCSd2r_0_1_2_En_%i",j));
        
        nfound = s->Search ( hCSd2r_0_1_2_En[j],2,"",0.1 );
       
        for ( int p=0; p<nfound; p++ ){Sd2rxp[p] = txp[p];}
       
        Sd2ra[j] = Sd2rxp[0];
    
    }
    
    
    for (int i=0;i<length;i++){
        
        std::string sd1Title = hCSd1r_0_1_2_En[i]->GetTitle();
        sd1Title.erase(sd1Title.begin(), sd1Title.end() - 4);
        runNumber[i] = std::stoi(sd1Title);
        // std::cout << sd1Title << '\t' << runNumber << std::endl;
        
        //With Target
        
        fit_gauss1->SetParLimits (1,8,10);
        hCSd1r_0_1_2_En[i]->Fit ("fit_gauss1","Q","",7,11);
        ESd1Si_Target[i]=fit_gauss1->GetParameter(1);//24.70;
        errorSd1r[i]=fit_gauss1->GetParError(1);
        
        fit_gauss1->SetParLimits (1,Sd2ra[i]-1,Sd2ra[i]+1);
        hCSd2r_0_1_2_En[i]->Fit ("fit_gauss1","Q","",Sd2ra[i]-3,Sd2ra[i]+1);
        ESd2Si_Target[i]=fit_gauss1->GetParameter(1);//76.88;
        errorSd2r[i]=fit_gauss1->GetParError(1);
        
        E_Be_Sd2_P_Target[i]=Be_P->AddBack(ESd2Si_Target[i],dSdP);
        E_Be_Sd2_Al3_Target[i]=Be_Al->AddBack(E_Be_Sd2_P_Target[i],dSdAl3);
        
        E_Be_Sd1_Al3_Target[i]=Be_Al->AddBack(E_Be_Sd2_Al3_Target[i],dSdAl3);
        E_Be_Sd1_P_Target[i]=Be_P->AddBack(E_Be_Sd1_Al3_Target[i],dSdP);
        Sd1rE_Target[i]=ESd1Si_Target[i]+E_Be_Sd1_P_Target[i];
        
        errorTotal[i]=TMath::Sqrt(errorSd1r[i]*errorSd1r[i]+errorSd2r[i]*errorSd2r[i]);
        
        E_Be_B_Target[i] = Be_B->AddBack(Sd1rE_Target[i],dSdB);
        E_Be_Al2_Target[i] = Be_Al->AddBack(E_Be_B_Target[i],dSdAl2);
        E_Be_SiO2_Target[i] = Be_SiO2->AddBack(E_Be_Al2_Target[i],dSdSiO2);
        E_Be_Al1_Target[i] = Be_Al->AddBack(E_Be_SiO2_Target[i],dSdAl1);
        E_Be_before_Ag_Target[i] = Be_Ag->AddBack(E_Be_Al1_Target[i],dAg); //Final Energy
        
        thickness[i] = 1000*Be_D2->CalcRange(E_Be_before_Ag_noTarget,E_Be_before_Ag_Target[i]);// thickness in um.
        
        //cout << E_Be_before_Ag_Target[i] << endl;
        
        
        //cout << thickness << endl;
        hCSd1r_0_1_2_En[i]->GetXaxis()->SetRangeUser(6,12);
        hCSd1r_0_1_2_En[i]->Write();
        hCSd2r_0_1_2_En[i]->GetXaxis()->SetRangeUser(96,106);
        hCSd2r_0_1_2_En[i]->Write();
    }
    
    
    /*
    TCanvas *c4 = new TCanvas ( "c4" );
    c4->Divide(9,9);
    
    TCanvas *c5 = new TCanvas ( "c5" );
    c5->Divide(9,9);
    
    for (int i=0;i<140;i++){
        c4->cd(i+1);
        hCSd1r_0_1_2_En[i]->Draw();
        c5->cd(i+1);
        hCSd2r_0_1_2_En[i]->Draw();
    }
    */
    g1->Close();
    
    int index =0;
    for (int i=0;i<length;i++){
        if (runNumber[i]==5100){
            continue;
        }else{
            gr->SetPoint(index,runNumber[i],E_Be_before_Ag_Target[i]);
            gr->SetPointError(index,0,errorTotal[i]);
            
            gr1->SetPoint(index,runNumber[i],thickness[i]);
            gr1->SetPointError(index,0,0);
            
            gr2->SetPoint(index,runNumber[i],ESd1Si_Target[i]);
            gr2->SetPointError(index,0,0);
            
            gr3->SetPoint(index,runNumber[i],ESd2Si_Target[i]);
            gr3->SetPointError(index,0,0);
            
            //cout << runNumber[i] << endl;
            index = index+1;
        }
    }
    TCanvas *c6 = new TCanvas ( "c6" );
    gr->Draw("AP*");
    TCanvas *c7 = new TCanvas ( "c7" );
    gr1->Draw("AP*");
    TCanvas *c8 = new TCanvas ( "c8" );
    gr2->Draw("AP*");
    TCanvas *c9 = new TCanvas ( "c9" );
    gr3->Draw("AP*");
    
    
    
}
