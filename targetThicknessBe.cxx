using namespace std;


#include <TMath.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include "EnergyLoss.h"


double gauss1(double *x, double *p) {
    return p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2]));
}

void targetThicknessBe(){
    
    
    TFile *g1 = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2.root","READ");
    //TFile *g1 = new TFile("Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2_target_thickness.root
    TH1D *hBeSd1r = (TH1D*)g1->Get("hCSd1rEn");
    
    TF1* fit_gauss1 = new TF1("fit_gauss1",gauss1,-10,10,3);
    
    hBeSd1r->Fit ( "fit_gauss1","","",9,11);

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
    
    
    //cout << ESi_noTarget << endl;
    //cout << E_Be_Si << endl;
    
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
    
    //With Target
    
    double ESd1Si_Target=9.067;//24.70;
    double ESd2Si_Target=100.87;//76.88;
    
    double E_Be_Sd2_P_Target=Be_P->AddBack(ESd2Si_Target,dSdP);
    double E_Be_Sd2_Al3_Target=Be_Al->AddBack(E_Be_Sd2_P_Target,dSdAl3);
    
    double E_Be_Sd1_Al3_Target=Be_Al->AddBack(E_Be_Sd2_Al3_Target,dSdAl3);
    double E_Be_Sd1_P_Target=Be_P->AddBack(E_Be_Sd1_Al3_Target,dSdP);
    double Sd1rE_Target=ESd1Si_Target+E_Be_Sd1_P_Target;
    double E_Be_B_Target = Be_B->AddBack(Sd1rE_Target,dSdB);
    double E_Be_Al2_Target = Be_Al->AddBack(E_Be_B_Target,dSdAl2);
    double E_Be_SiO2_Target = Be_SiO2->AddBack(E_Be_Al2_Target,dSdSiO2);
    double E_Be_Al1_Target = Be_Al->AddBack(E_Be_SiO2_Target,dSdAl1);
    double E_Be_before_Ag_Target = Be_Ag->AddBack(E_Be_Al1_Target,dAg); //Final Energy
    
    double thickness = Be_D2->CalcRange(E_Be_before_Ag_noTarget,E_Be_before_Ag_Target);
    
    cout << E_Be_before_Ag_Target << endl;
    cout << thickness << endl;
    
    
    
    
}
