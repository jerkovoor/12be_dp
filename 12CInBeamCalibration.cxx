//In-beam calibration of Sd1 and Sd2 using 12C beam and 4He particles

using namespace std;


#include <TMath.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include "EnergyLoss.h"


double gauss1(double *x, double *p) {
    return p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2]));
}

void CInBeamCalibration(){
    
    
    TFile *g1 = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2.root","READ");
    //TFile *g1 = new TFile("Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2_target_thickness.root
    TH1D *hBeSd1r = (TH1D*)g1->Get("hCSd1rEn");
    
    TF1* fit_gauss1 = new TF1("fit_gauss1",gauss1,-10,10,3);
    
    //hBeSd1r->Fit ( "fit_gauss1","","",9,11);

    EnergyLoss* C_Al = new EnergyLoss("C_in_Al.dat");
    EnergyLoss* C_SiO2 = new EnergyLoss("C_in_SiO2.dat");
	EnergyLoss* C_B = new EnergyLoss("C_in_B.dat");
    EnergyLoss* C_Si = new EnergyLoss("C_in_Si.dat");
    EnergyLoss* C_P = new EnergyLoss("C_in_P.dat");
    
    
    /*He_Al->AddBackHigh(15.);
    He_SiO2->AddBackHigh(15.);
    He_B->AddBackHigh(15.);
    */
    
    double Sd1rAngles[24];
    double Sd2rAngles[24];
    
    double SD1_DISTANCE = 600.;
    double SD2_DISTANCE = 690.;
    double SD_INNER_RADIUS = 11.;
    double SD_OUTER_RADIUS = 35.;
    
    double width_Sdr=1;
    
    for (int i=0;i<24;i++){
        Sd1rAngles[i]=TMath::RadToDeg()*TMath::ATan((SD_INNER_RADIUS+(i+0.5)*width_Sdr)/(SD1_DISTANCE));
        Sd2rAngles[i]=TMath::RadToDeg()*TMath::ATan((SD_INNER_RADIUS+(i+0.5)*width_Sdr)/(SD2_DISTANCE));
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
    

    
    double EnHeSd1=2.2944;
    double EnHeSd2=32.745;
    
    double EnC=107.01;  //Energy of Carbon beam after passing through the silver backing
    
    

    
    double EnHe_Al1[24], EnHe_SiO2[24], EnHe_Al2[24], EnHe_B[24];
    double EnC_Al1[24], EnC_SiO2[24], EnC_Al2[24], EnC_B[24], EnC_Si[24], EnC_P[24], EnC_Al3[24];
    double EnC_Al3_Sd2[24], EnC_P_Sd2[24];
    double EnCSd1[24], EnCSd2[24];

    
    for (int i=0;i<24;i++){
    /*
        EnHe_Al1[i] = He_Al->CalcRemainder(EnHe,dSdAl1/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        EnHe_SiO2[i] = He_SiO2->CalcRemainder(EnHe_Al1[i],dSdSiO2/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        EnHe_Al2[i] = He_Al->CalcRemainder(EnHe_SiO2[i],dSdAl2/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        EnHe_B[i] = He_B->CalcRemainder(EnHe_Al2[i],dSdB/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
     */ 
    //Decided not to consider angle effects in the energy loss of 4He particles.
    
        ///////////////////Sd1////////////////////
        
        EnC_Al1[i] = C_Al->CalcRemainder(EnC,dSdAl1/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        EnC_SiO2[i] = C_SiO2->CalcRemainder(EnC_Al1[i],dSdSiO2/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        EnC_Al2[i] = C_Al->CalcRemainder(EnC_SiO2[i],dSdAl2/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        EnC_B[i] = C_B->CalcRemainder(EnC_Al2[i],dSdB/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        EnC_Si[i] = C_Si->CalcRemainder(EnC_B[i],dSd1Si/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        
        EnCSd1[i]=EnC_B[i]-EnC_Si[i]; //Energy recorded in Sd1
        
        EnC_P[i] = C_P->CalcRemainder(EnC_Si[i],dSdP/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        EnC_Al3[i] = C_Al->CalcRemainder(EnC_P[i],dSdAl3/fabs(cos(Sd1rAngles[i]*M_PI/180.)));
        
        ///////////////////Sd2////////////////////
        EnC_Al3_Sd2[i] = C_Al->CalcRemainder(EnC_Al3[i],dSdAl3/fabs(cos(Sd2rAngles[i]*M_PI/180.)));
        EnC_P_Sd2[i] = C_P->CalcRemainder(EnC_Al3_Sd2[i],dSdP/fabs(cos(Sd2rAngles[i]*M_PI/180.)));
        
        EnCSd2[i] = EnC_P_Sd2[i];  //Full energy is deposited in the Sd2
        
        
        //cout << EnHe_B[i] << endl;
    }
    
    
    
    //4.9098 4.3241 4.2696 4.1663
    //cout << EnHe_Al1[1] << "  " << EnHe_SiO2[1] << "  " << EnHe_Al2[1] << "  " << EnHe_B[1] << endl;
    
    TF1 *fit_lin = new TF1 ( "fit_lin","[0]*x+[1]",0,4000 );
    
    
    
    /////////////////////////////Sd1r Calibration////////////////////////////////
    
    double stripSd1r[24], peakX1Sd1r[24], peakX2Sd1r[24];
    
    ////////////He peaks///////////////
    
    ifstream He_peaks_Sd1r;
    He_peaks_Sd1r.open("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1r_4He_Peaks.txt");
    
    
    if(He_peaks_Sd1r.is_open()){
        He_peaks_Sd1r.ignore(256,'\n');
        for(int i=0;i<24;i++){
            He_peaks_Sd1r >> stripSd1r[i] >> peakX1Sd1r[i];
        }
    }else{
        cout << "No file found " << endl;
    }
    He_peaks_Sd1r.close();
    
     ////////////C peaks///////////////
    
    ifstream C_peaks_Sd1r;
    C_peaks_Sd1r.open("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1r_12C_Peaks.txt");
    
    
    if(C_peaks_Sd1r.is_open()){
        C_peaks_Sd1r.ignore(256,'\n');
        for(int i=0;i<24;i++){
            C_peaks_Sd1r >> stripSd1r[i] >> peakX2Sd1r[i];
        }
    }else{
        cout << "No file found " << endl;
    }
    C_peaks_Sd1r.close();
    
    
    ///////////Fitting////////////////
    
    TGraph* grSd1r[24];
    
    for(int i=0;i<24;i++){
        double EexpSd1r[2]={peakX1Sd1r[i],peakX2Sd1r[i]};
        double EtheoSd1r[2]={1000*EnHeSd1,1000*EnCSd1[i]}; //Energy in keV
            
        grSd1r[i] = new TGraph(2,EexpSd1r,EtheoSd1r);
    }
    
    
    ofstream Sd1r_InBeamCal;
    Sd1r_InBeamCal.open("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1rInBeamCalibration.txt");
    Sd1r_InBeamCal << " stripSd1r / gain / offset  \n";


    TCanvas *c2 = new TCanvas ( "c2" );
    c2->Divide ( 6,4 );
        
        
    for ( int i=0; i<24; i++ ){
        c2->cd ( i+1 );

        grSd1r[i]->Draw ( "AP*" );
        grSd1r[i]->Fit ( "fit_lin","Q","",0,4000 );
        Sd1r_InBeamCal << i << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;
        
    }
    
    TFile *fSd1r_out = new TFile("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1rInBeamCalibration.root","RECREATE");
    
    fSd1r_out->cd();
    
    for ( int i=0; i<24; i++ ){
        string name2 = Form ( "3points_%i",i );
        grSd1r[i]->Write(name2.c_str());
    }
      
    fSd1r_out->Close();
    
    
    
    
    
    
    
    
    
    /////////////////////////////Sd1s Calibration////////////////////////////////
    
    
    double stripSd1s[32], peakX1Sd1s[32], peakX2Sd1s[32];
    
    ////////////He peaks///////////////
    
    
    ifstream He_peaks_Sd1s;
    He_peaks_Sd1s.open("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1s_4He_Peaks.txt");
    
    
    if(He_peaks_Sd1s.is_open()){
        He_peaks_Sd1s.ignore(256,'\n');
        for(int i=0;i<32;i++){
            He_peaks_Sd1s >> stripSd1s[i] >> peakX1Sd1s[i];
        }
    }else{
        cout << "No file found " << endl;
    }
    He_peaks_Sd1s.close();
    
    ////////////C peaks///////////////
    
    ifstream C_peaks_Sd1s;
    C_peaks_Sd1s.open("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1s_12C_Peaks.txt");
    
    
    if(C_peaks_Sd1s.is_open()){
        C_peaks_Sd1s.ignore(256,'\n');
        for(int i=0;i<24;i++){
            C_peaks_Sd1s >> stripSd1s[i] >> peakX2Sd1s[i];
        }
    }else{
        cout << "No file found " << endl;
    }
    C_peaks_Sd1s.close();
    
    ///////////Fitting/////////////
    
    TGraph* grSd1s[32];
    
    for(int i=0;i<32;i++){
        double EexpSd1s[2]={peakX1Sd1s[i],peakX2Sd1s[i]};
        double EtheoSd1s[2]={1000*EnHeSd1,1000*EnCSd1[i]}; //Energy in keV
            
        grSd1s[i] = new TGraph(2,EexpSd1s,EtheoSd1s);
    }
    
    
    ofstream Sd1s_alphaCal;
    Sd1s_alphaCal.open("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1sInBeamCalibration.txt");
    Sd1s_alphaCal << " stripSd1s / gain / offset  \n";


    TCanvas *c3 = new TCanvas ( "c3" );
    c3->Divide ( 8,4 );
        
        
    for ( int i=0; i<32; i++ ){
        c3->cd ( i+1 );

        grSd1s[i]->Draw ( "AP*" );
        grSd1s[i]->Fit ( "fit_lin","Q","",0,4000 );
        Sd1s_alphaCal << i << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;
        
    }
    
    TFile *fSd1s_out = new TFile("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1sInBeamCalibration.root","RECREATE");
    
    fSd1s_out->cd();
    
    for ( int i=0; i<32; i++ ){
        string name2 = Form ( "3points_%i",i );
        grSd1s[i]->Write(name2.c_str());
    }
      
    fSd1s_out->Close();
    
    
    
    
    
    
    
    /////////////////////////////Sd2r Calibration////////////////////////////////
    
    double stripSd2r[24], peakX1Sd2r[24], peakX2Sd2r[24];
    
    ////////////He peaks///////////////
    
    ifstream He_peaks_Sd2r;
    He_peaks_Sd2r.open("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd2r_4He_Peaks.txt");
    
    
    if(He_peaks_Sd2r.is_open()){
        He_peaks_Sd2r.ignore(256,'\n');
        for(int i=0;i<24;i++){
            He_peaks_Sd2r >> stripSd2r[i] >> peakX1Sd2r[i];
        }
    }else{
        cout << "No file found " << endl;
    }
    He_peaks_Sd2r.close();
    
     ////////////C peaks///////////////
    
    ifstream C_peaks_Sd2r;
    C_peaks_Sd2r.open("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd2r_12C_Peaks.txt");
    
    
    if(C_peaks_Sd2r.is_open()){
        C_peaks_Sd2r.ignore(256,'\n');
        for(int i=0;i<24;i++){
            C_peaks_Sd2r >> stripSd2r[i] >> peakX2Sd2r[i];
        }
    }else{
        cout << "No file found " << endl;
    }
    C_peaks_Sd2r.close();
    
    
    ///////////Fitting////////////////
    
    TGraph* grSd2r[24];
    
    for(int i=0;i<24;i++){
        double EexpSd2r[2]={peakX1Sd2r[i],peakX2Sd2r[i]};
        double EtheoSd2r[2]={1000*EnHeSd2,1000*EnCSd2[i]}; //Energy in keV
            
        grSd2r[i] = new TGraph(2,EexpSd2r,EtheoSd2r);
    }
    
    
    ofstream Sd2r_InBeamCal;
    Sd2r_InBeamCal.open("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd2rInBeamCalibration.txt");
    Sd2r_InBeamCal << " stripSd2r / gain / offset  \n";


    TCanvas *c4 = new TCanvas ( "c4" );
    c4->Divide ( 6,4 );
        
        
    for ( int i=0; i<24; i++ ){
        c4->cd ( i+1 );

        grSd2r[i]->Draw ( "AP*" );
        grSd2r[i]->Fit ( "fit_lin","Q","",0,4000 );
        Sd2r_InBeamCal << i << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;
        
    }
    
    TFile *fSd2r_out = new TFile("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd2rInBeamCalibration.root","RECREATE");
    
    fSd2r_out->cd();
    
    for ( int i=0; i<24; i++ ){
        string name2 = Form ( "3points_%i",i );
        grSd2r[i]->Write(name2.c_str());
    }
      
    fSd2r_out->Close();
    
    
    
    
    
    
    
    
    
    /////////////////////////////Sd2s Calibration////////////////////////////////
    
    
    double stripSd2s[32], peakX1Sd2s[32], peakX2Sd2s[32];
    
    ////////////He peaks///////////////
    
    
    ifstream He_peaks_Sd2s;
    He_peaks_Sd2s.open("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd2s_4He_Peaks.txt");
    
    
    if(He_peaks_Sd2s.is_open()){
        He_peaks_Sd2s.ignore(256,'\n');
        for(int i=0;i<32;i++){
            He_peaks_Sd2s >> stripSd2s[i] >> peakX1Sd2s[i];
        }
    }else{
        cout << "No file found " << endl;
    }
    He_peaks_Sd2s.close();
    
    ////////////C peaks///////////////
    
    ifstream C_peaks_Sd2s;
    C_peaks_Sd2s.open("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd2s_12C_Peaks.txt");
    
    
    if(C_peaks_Sd2s.is_open()){
        C_peaks_Sd2s.ignore(256,'\n');
        for(int i=0;i<24;i++){
            C_peaks_Sd2s >> stripSd2s[i] >> peakX2Sd2s[i];
        }
    }else{
        cout << "No file found " << endl;
    }
    C_peaks_Sd2s.close();
    
    ///////////Fitting/////////////
    
    TGraph* grSd2s[32];
    
    for(int i=0;i<32;i++){
        double EexpSd2s[2]={peakX1Sd2s[i],peakX2Sd2s[i]};
        double EtheoSd2s[2]={1000*EnHeSd2,1000*EnCSd2[i]}; //Energy in keV
            
        grSd2s[i] = new TGraph(2,EexpSd2s,EtheoSd2s);
    }
    
    
    ofstream Sd2s_alphaCal;
    Sd2s_alphaCal.open("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd2sInBeamCalibration.txt");
    Sd2s_alphaCal << " stripSd2s / gain / offset  \n";


    TCanvas *c5 = new TCanvas ( "c5" );
    c5->Divide ( 8,4 );
        
        
    for ( int i=0; i<32; i++ ){
        c5->cd ( i+1 );

        grSd2s[i]->Draw ( "AP*" );
        grSd2s[i]->Fit ( "fit_lin","Q","",0,4000 );
        Sd2s_alphaCal << i << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;
        
    }
    
    TFile *fSd2s_out = new TFile("/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd2sInBeamCalibration.root","RECREATE");
    
    fSd2s_out->cd();
    
    for ( int i=0; i<32; i++ ){
        string name2 = Form ( "3points_%i",i );
        grSd2s[i]->Write(name2.c_str());
    }
      
    fSd2s_out->Close();
        
    
}
