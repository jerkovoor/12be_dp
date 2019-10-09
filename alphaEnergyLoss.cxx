using namespace std;


#include <TMath.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include "EnergyLoss.h"


double gauss1(double *x, double *p) {
    return p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2]));
}

void alphaEnergyLoss(){
    
    
    TFile *g1 = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2.root","READ");
    //TFile *g1 = new TFile("Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2_target_thickness.root
    TH1D *hBeSd1r = (TH1D*)g1->Get("hCSd1rEn");
    
    TF1* fit_gauss1 = new TF1("fit_gauss1",gauss1,-10,10,3);
    
    //hBeSd1r->Fit ( "fit_gauss1","","",9,11);

    EnergyLoss* He_Al = new EnergyLoss("He_in_Al.dat");
    EnergyLoss* He_SiO2 = new EnergyLoss("He_in_SiO2.dat");
	EnergyLoss* He_B = new EnergyLoss("He_in_B.dat");
    EnergyLoss* He_Si = new EnergyLoss("He_in_Si.dat");
    EnergyLoss* He_P = new EnergyLoss("He_in_P.dat");
    
    
    /*He_Al->AddBackHigh(15.);
    He_SiO2->AddBackHigh(15.);
    He_B->AddBackHigh(15.);
    */
    
    double Sd1rAngles[24];
    
    double SD1_DISTANCE = 600.;
    double SD_INNER_RADIUS = 11.;
    double SD_OUTER_RADIUS = 35.;
    
    double width_Sd1r=1;
    
    for (int i=0;i<24;i++){
        Sd1rAngles[i]=TMath::RadToDeg()*TMath::ATan((SD_INNER_RADIUS+(i+0.5)*width_Sd1r)/(SD1_DISTANCE));
        //cout << Sd1rAngles[i] << endl;
    }
    
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
    

    
    double alpha1=5.156;
    double alpha2=5.486;
    double alpha3=5.805;
    
    double alpha1_Al1[24], alpha1_SiO2[24], alpha1_Al2[24], alpha1_B[24];
    double alpha2_Al1[24], alpha2_SiO2[24], alpha2_Al2[24], alpha2_B[24];
    double alpha3_Al1[24], alpha3_SiO2[24], alpha3_Al2[24], alpha3_B[24];
    
    
    for (int i=0;i<24;i++){
    
        alpha1_Al1[i] = He_Al->CalcRemainder(alpha1,dSdAl1/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        alpha1_SiO2[i] = He_SiO2->CalcRemainder(alpha1_Al1[i],dSdSiO2/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        alpha1_Al2[i] = He_Al->CalcRemainder(alpha1_SiO2[i],dSdAl2/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        alpha1_B[i] = He_B->CalcRemainder(alpha1_Al2[i],dSdB/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        
        alpha2_Al1[i] = He_Al->CalcRemainder(alpha2,dSdAl1/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        alpha2_SiO2[i] = He_SiO2->CalcRemainder(alpha2_Al1[i],dSdSiO2/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        alpha2_Al2[i] = He_Al->CalcRemainder(alpha2_SiO2[i],dSdAl2/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        alpha2_B[i] = He_B->CalcRemainder(alpha2_Al2[i],dSdB/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        
        alpha3_Al1[i] = He_Al->CalcRemainder(alpha3,dSdAl1/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        alpha3_SiO2[i] = He_SiO2->CalcRemainder(alpha3_Al1[i],dSdSiO2/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        alpha3_Al2[i] = He_Al->CalcRemainder(alpha3_SiO2[i],dSdAl2/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        alpha3_B[i] = He_B->CalcRemainder(alpha3_Al2[i],dSdB/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        
        //cout << alpha1_B[i] << endl;
    }
    
    //4.9098 4.3241 4.2696 4.1663
    //cout << alpha1_Al1[1] << "  " << alpha1_SiO2[1] << "  " << alpha1_Al2[1] << "  " << alpha1_B[1] << endl;
    
    TF1 *fit_lin = new TF1 ( "fit_lin","[0]*x+[1]",0,4000 );
    
    
    /////////////////////////////Sd1r Calibration////////////////////////////////
    
    double stripSd1r[24], peakX1Sd1r[24], peakX2Sd1r[24], peakX3Sd1r[24];
    
    TGraph* grSd1r[24];
    
    ifstream Alpha_peaks_Sd1r;
    Alpha_peaks_Sd1r.open("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1r_AlphaPeaks.txt");
    
    
    if(Alpha_peaks_Sd1r.is_open()){
        Alpha_peaks_Sd1r.ignore(256,'\n');
        for(int i=0;i<24;i++){
            Alpha_peaks_Sd1r >> stripSd1r[i] >> peakX1Sd1r[i] >> peakX2Sd1r[i] >> peakX3Sd1r[i];
            //Alpha_peaks_Sd1r >> stripSd1r[i] >> peakX1Sd1r[i] >> peakX2Sd1r[i] >> peakX3Sd1r[i];
            
            double EexpSd1r[3]={peakX1Sd1r[i],peakX2Sd1r[i],peakX3Sd1r[i]};
            double EtheoSd1r[3]={1000*alpha1_B[i],1000*alpha2_B[i],1000*alpha3_B[i]}; //Energy in keV
            
            grSd1r[i] = new TGraph(3,EexpSd1r,EtheoSd1r);
        }
    }
    
    else{cout << "No file found " << endl;}
    Alpha_peaks_Sd1r.close();
    
    
    ofstream Sd1r_alphaCal;
    Sd1r_alphaCal.open("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1rAlphaCalibration.txt");
    Sd1r_alphaCal << " stripSd1r / gain / offset  \n";


    TCanvas *c2 = new TCanvas ( "c2" );
    c2->Divide ( 6,4 );
        
        
    for ( int i=0; i<24; i++ ){
        c2->cd ( i+1 );

        grSd1r[i]->Draw ( "AP*" );
        grSd1r[i]->Fit ( "fit_lin","","",0,4000 );
        Sd1r_alphaCal << i << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;
        
    }
    
    TFile *fSd1r_out = new TFile("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1rAlphaCalibration.root","RECREATE");
    
    fSd1r_out->cd();
    
    for ( int i=0; i<24; i++ ){
        string name2 = Form ( "3points_%i",i );
        grSd1r[i]->Write(name2.c_str());
    }
      
    fSd1r_out->Close();
    
    
    /////////////////////////////Sd1s Calibration////////////////////////////////
    
    
    double stripSd1s[32], peakX1Sd1s[32], peakX2Sd1s[32], peakX3Sd1s[32];
    
    TGraph* grSd1s[32];
    
    ifstream Alpha_peaks_Sd1s;
    Alpha_peaks_Sd1s.open("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1s_AlphaPeaks.txt");
    
    
    if(Alpha_peaks_Sd1s.is_open()){
        Alpha_peaks_Sd1s.ignore(256,'\n');
        for(int i=0;i<32;i++){
            Alpha_peaks_Sd1s >> stripSd1s[i] >> peakX1Sd1s[i] >> peakX2Sd1s[i] >> peakX3Sd1s[i];
            //Alpha_peaks_Sd1s >> stripSd1s[i] >> peakX1Sd1s[i] >> peakX2Sd1s[i] >> peakX3Sd1s[i];
            
            double EexpSd1s[3]={peakX1Sd1s[i],peakX2Sd1s[i],peakX3Sd1s[i]};
            double EtheoSd1s[3]={1000*alpha1_B[1],1000*alpha2_B[1],1000*alpha3_B[1]}; //Energy in keV
            
            grSd1s[i] = new TGraph(3,EexpSd1s,EtheoSd1s);
        }
    }
    
    else{cout << "No file found " << endl;}
    Alpha_peaks_Sd1s.close();
    
    
    ofstream Sd1s_alphaCal;
    Sd1s_alphaCal.open("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1sAlphaCalibration.txt");
    Sd1s_alphaCal << " stripSd1s / gain / offset  \n";


    TCanvas *c3 = new TCanvas ( "c3" );
    c3->Divide ( 8,4 );
        
        
    for ( int i=0; i<32; i++ ){
        c3->cd ( i+1 );

        grSd1s[i]->Draw ( "AP*" );
        grSd1s[i]->Fit ( "fit_lin","","",0,4000 );
        Sd1s_alphaCal << i << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;
        
    }
    
    TFile *fSd1s_out = new TFile("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1sAlphaCalibration.root","RECREATE");
    
    fSd1s_out->cd();
    
    for ( int i=0; i<32; i++ ){
        string name2 = Form ( "3points_%i",i );
        grSd1s[i]->Write(name2.c_str());
    }
      
    fSd1s_out->Close();
        
    
}
