using namespace std;


#include <TMath.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include "EnergyLoss.h"


double gauss1(double *x, double *p) {
    return p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2]));
}

void targetThicknessC(){
    
    
    TFile *g1 = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2.root","READ");
    //TFile *g1 = new TFile("Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2_target_thickness.root
    TH1D *hCSd1r = (TH1D*)g1->Get("hCSd1rEn");
    
    TF1* fit_gauss1 = new TF1("fit_gauss1",gauss1,-10,10,3);
    
    hCSd1r->Fit ( "fit_gauss1","","",9,11);

    EnergyLoss* C_Al = new EnergyLoss("C_in_Al.dat");
    EnergyLoss* C_SiO2 = new EnergyLoss("C_in_SiO2.dat");
	EnergyLoss* C_B = new EnergyLoss("C_in_B.dat");
    EnergyLoss* C_Si = new EnergyLoss("C_in_Si.dat");
    EnergyLoss* C_P = new EnergyLoss("C_in_P.dat");
    EnergyLoss* C_Ag = new EnergyLoss("C_in_Ag.dat");
    EnergyLoss* C_D2 = new EnergyLoss("C_in_D2_density0_201.dat");
    
    /*C_Al->AddBackHigh(15.);
    C_SiO2->AddBackHigh(15.);
    C_B->AddBackHigh(15.);
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
    

    
    double E_C_after_Ag=108.57;
    
    
    //cout << ESi_noTarget << endl;
    //cout << E_C_Si << endl;
    
    //Without Target
    
    double energy = 77.44;//77.05;
    
    double ESd1Si_noTarget=21.54;//24.54;
    double ESd2Si_noTarget=83.19;//77.05;
    
    //cout << "sd2 Energy: " << energy << endl;
    
    double E_C_Sd2_P_noTarget=C_P->AddBack(ESd2Si_noTarget,dSdP);
    energy = C_P->AddBack(energy, dSdP);
    
    //cout << "energy before P dead layer (sd2): " << energy << endl;
    
    double E_C_Sd2_Al3_noTarget=C_Al->AddBack(E_C_Sd2_P_noTarget,dSdAl3);
    energy = C_Al->AddBack(energy, dSdAl3);
    
    //cout << "energy before Al dead layer (sd2): " << energy << endl;
    
    double E_C_Sd1_Al3_noTarget=C_Al->AddBack(E_C_Sd2_Al3_noTarget,dSdAl3);
    energy = C_Al->AddBack(energy, dSdAl3);
    
    //cout << "energy before Al dead layer (sd1): " << energy << endl;
    
    double E_C_Sd1_P_noTarget=C_P->AddBack(E_C_Sd1_Al3_noTarget,dSdP);
    energy = C_P->AddBack(energy, dSdP);
    
    //cout << "energy before P dead layer (sd1): " << energy << endl;
    
    double Sd1rE_noTarget=ESd1Si_noTarget+E_C_Sd1_P_noTarget;
    energy = ESd1Si_noTarget + energy;
    
    //cout << "sd1 energy: " << ESd1Si_noTarget << endl;
    //cout << "sd1 energy + sd2 energy with dead layers: " << energy << endl;
    
    //cout << "sd1 deposited from SRIM: " << 101.8 - C_Si->CalcRemainder(101.8, dSd1Si) << endl;
    
    double E_C_B_noTarget = C_B->AddBack(Sd1rE_noTarget,dSdB);
    double E_C_Al2_noTarget = C_Al->AddBack(E_C_B_noTarget,dSdAl2);
    double E_C_SiO2_noTarget = C_SiO2->AddBack(E_C_Al2_noTarget,dSdSiO2);
    double E_C_Al1_noTarget = C_Al->AddBack(E_C_SiO2_noTarget,dSdAl1);
    double E_C_before_Ag_noTarget = C_Ag->AddBack(E_C_Al1_noTarget,dAg); //Initial Energy
    
    //cout << E_C_Sd2_Al3_noTarget << endl;
    
    cout << E_C_before_Ag_noTarget << endl;
    
    //With Target
    
    double ESd1Si_Target=22.00;//24.70;
    double ESd2Si_Target=80.44;//76.88;
    
    double E_C_Sd2_P_Target=C_P->AddBack(ESd2Si_Target,dSdP);
    double E_C_Sd2_Al3_Target=C_Al->AddBack(E_C_Sd2_P_Target,dSdAl3);
    
    double E_C_Sd1_Al3_Target=C_Al->AddBack(E_C_Sd2_Al3_Target,dSdAl3);
    double E_C_Sd1_P_Target=C_P->AddBack(E_C_Sd1_Al3_Target,dSdP);
    double Sd1rE_Target=ESd1Si_Target+E_C_Sd1_P_Target;
    double E_C_B_Target = C_B->AddBack(Sd1rE_Target,dSdB);
    double E_C_Al2_Target = C_Al->AddBack(E_C_B_Target,dSdAl2);
    double E_C_SiO2_Target = C_SiO2->AddBack(E_C_Al2_Target,dSdSiO2);
    double E_C_Al1_Target = C_Al->AddBack(E_C_SiO2_Target,dSdAl1);
    double E_C_before_Ag_Target = C_Ag->AddBack(E_C_Al1_Target,dAg); //Final Energy
    
    double thickness = C_D2->CalcRange(E_C_before_Ag_noTarget,E_C_before_Ag_Target);
    
    cout << E_C_before_Ag_Target << endl;
    cout << thickness << endl;
    
    
    
    
}
