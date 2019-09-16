using namespace std;

#include "TFile.h" 
#include <cmath>
#include "iostream"
#include "fstream"
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TMath.h>

#include <Fit/Fitter.h>


double kin(double *x, double *Q_C){ 
	
	float amu = 931.5; // atomic mass unit in MeV 
	float massEjec = 938.28;//mass of the proton in MeV/c2 
	float kBeam = 112.21; //put the correct value; beam energy
	float mbeam = 12 * amu;  //mass of the beam (12Be or 12C) in MeV
	float mrecoil = 13 * amu;  //mass of the recoil (13Be or 13C) in MeV
	float mejec = 1 * amu; //mass of the proton
	
	return pow(((TMath::Sqrt(mbeam*amu*kBeam)*TMath::Cos(TMath::DegToRad()*x[0])+TMath::Sqrt(mbeam*amu*kBeam*pow(TMath::Cos(TMath::DegToRad()*x[0]),2)+(mrecoil+amu)*(mrecoil*Q_C[0]+(mrecoil-mbeam)*kBeam)))/(mrecoil+amu)),2); 
} 


double shaper(double *x, double *p) {
    return p[0]*exp(-(x[0] - p[1])*(x[0] - p[1])/(2.*p[2]*p[2])) + p[3]*exp(-(x[0] - p[4])*(x[0] - p[4])/(2.*p[2]*p[2])) + p[5]*exp(-(x[0] - p[6])*(x[0] - p[6])/(2.*p[2]*p[2]));
}

double shaper2(double *x, double *p) {
    return p[0]*exp(-(x[0] - p[1])*(x[0] - p[1])/(2.*p[2]*p[2])) + p[3]*exp(-(x[0] - p[4])*(x[0] - p[4])/(2.*p[2]*p[2])) + p[5]*exp(-(x[0] - p[6])*(x[0] - p[6])/(2.*p[2]*p[2])) + 
        + p[7]*exp(-(x[0] - p[8])*(x[0] - p[8])/(2.*p[2]*p[2]));
}

double exp_background(double *x, double *p) {
    return p[0]*TMath::Exp(p[1]*x[0]);
}

double poly6_background(double *x, double *p) {
    return p[0]*pow(x[0],6);
}


double gauss1(double *x, double *p) {
    return p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2]));
}

/*double voig1(double *x, double *p) {
    return p[5]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[5])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)));//+p[6]*TMath::Exp(p[7]*x[0]);
}*/

double voig1(double *x, double *p) {
    return p[5]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[5])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)));
}

/*double totalFit(double *x, double *p){
    return voig1(x,p)+voig1(x,p[8])+voig1(x,&p[16])+voig1(x,&p[24]);
}*/

double voig2(double *x, double *p) {
    return p[9]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[9])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[9]*(p[5]*exp(-pow((x[0]-p[6]),2)/(2*p[2]*p[2])))+(1-p[9])*(p[7]/(pow((x[0]-p[6]),2)+pow((p[8]/2),2)))+p[10]*TMath::Exp(p[11]*x[0]);
}

/*double voig3(double *x, double *p) {
    return p[13]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[13])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[13]*(p[5]*exp(-pow((x[0]-p[6]),2)/(2*p[2]*p[2])))+(1-p[13])*(p[7]/(pow((x[0]-p[6]),2)+pow((p[8]/2),2)))+p[13]*(p[9]*exp(-pow((x[0]-p[10]),2)/(2*p[2]*p[2])))+(1-p[13])*(p[11]/(pow((x[0]-p[10]),2)+pow((p[12]/2),2)))+p[14]*TMath::Exp(p[15]*x[0]);
}*/

double voig3(double *x, double *p) {
    return p[5]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[5])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[11]*(p[6]*exp(-pow((x[0]-p[7]),2)/(2*p[8]*p[9])))+(1-p[11])*(p[9]/(pow((x[0]-p[7]),2)+pow((p[10]/2),2)))+p[17]*(p[12]*exp(-pow((x[0]-p[13]),2)/(2*p[14]*p[14])))+(1-p[17])*(p[15]/(pow((x[0]-p[13]),2)+pow((p[16]/2),2)))+p[18]*TMath::Exp(p[19]*x[0]);
}

double voig4(double *x, double *p) {
    return p[17]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[17])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[18]*(p[5]*exp(-pow((x[0]-p[6]),2)/(2*p[2]*p[2])))+(1-p[18])*(p[7]/(pow((x[0]-p[6]),2)+pow((p[8]/2),2)))+p[19]*(p[9]*exp(-pow((x[0]-p[10]),2)/(2*p[2]*p[2])))+(1-p[19])*(p[11]/(pow((x[0]-p[10]),2)+pow((p[12]/2),2)))+p[20]*(p[13]*exp(-pow((x[0]-p[14]),2)/(2*p[21]*p[21])))+(1-p[20])*(p[15]/(pow((x[0]-p[14]),2)+pow((p[16]/2),2)))+p[22]*pow(x[0],6);
}


double lorentz3gauss4(double *x, double *p) {
    return p[16]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[16])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[16]*(p[5]*exp(-pow((x[0]-p[6]),2)/(2*p[2]*p[2])))+p[16]*(p[7]*exp(-pow((x[0]-p[8]),2)/(2*p[2]*p[2])))+(1-p[16])*(p[9]/(pow((x[0]-p[8]),2)+pow((p[10]/2),2)))+p[16]*(p[11]*exp(-pow((x[0]-p[12]),2)/(2*p[13]*p[13])))+(1-p[16])*(p[14]/(pow((x[0]-p[12]),2)+pow((p[15]/2),2)))+p[17]*TMath::Exp(p[18]*x[0]);
}
/*double voig4(double *x, double *p) {
    return p[5]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[5])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[11]*(p[6]*exp(-pow((x[0]-p[7]),2)/(2*p[8]*p[8])))+(1-p[11])*(p[9]/(pow((x[0]-p[7]),2)+pow((p[10]/2),2)))+p[17]*(p[12]*exp(-pow((x[0]-p[13]),2)/(2*p[14]*p[14])))+(1-p[17])*(p[15]/(pow((x[0]-p[13]),2)+pow((p[16]/2),2)))+p[23]*(p[18]*exp(-pow((x[0]-p[19]),2)/(2*p[20]*p[20])))+(1-p[23])*(p[21]/(pow((x[0]-p[19]),2)+pow((p[22]/2),2)))+p[24]*TMath::Exp(p[25]*x[0]);
}*/


void BeQval(){
    float Qgs13C = 2.722; // Ground state Q value of 13C
	float QFES13C = -0.367;
    float QSES13C = -0.962;
	float Exenergy = 3.853; //Third excited state energy
	float QTES13C = Qgs13C-Exenergy;
    double QC0[2] = {Qgs13C,Qgs13C};
    double QC1[2] = {QFES13C,QFES13C};
    double QC2[2] = {QSES13C,QSES13C};
    double QC3[2] = {QTES13C,QTES13C};
    double countRange[2] = {0,50};
    double angle[17]={0};
    float TDistance = 80.88;
    
    TFile *f = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/C_pedestal_TRIUMF_DL_NoOffset_Shift0_0883_Yu_TargetDistance80_87mm.root","READ");
	//TFile *f = new TFile("/home/jerome/12Be_exp/Analysis/CarbonGain/C_pedestal_TRIUMF_DL_QvalMinGS_Yu.root","READ");
    TFile *fg = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod.root","READ");
	TFile *g = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n2.root","READ");
    TFile *g1 = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2.root","READ");
    TFile *g4 = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1.root","READ");
    TFile *gh1 = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2_FirstHalf.root","READ");
    TFile *gh2 = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2_SecondHalf.root","READ");


	TF1* fit_func2 = new TF1 ( "fit_func2","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp([6]*x)",-10,10);
    TF1* fit_func3 = new TF1 ( "fit_func3","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp([8]*x)",-10,10);
	TF1* fit_func4 = new TF1 ( "fit_func4","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp([10]*x)",-10,10);
    TF1* lor2 = new TF1 ("lor2","([0]*[2]*[2])/(pow((x-[1]),2)+pow(([2]/2),2))+([3]*[2]*[2])/(pow((x-[4]),2)+pow(([5]/2),2))+[6]*TMath::Exp([7]*x)",-10,10);
    TF1* lor3 = new TF1 ("lor3","([0]*[2]*[2])/(pow((x-[1]),2)+pow(([2]/2),2))+([3]*[2]*[2])/(pow((x-[4]),2)+pow(([2]/2),2))+([5]*[2]*[2])/(pow((x-[6]),2)+pow(([2]/2),2))+[7]*TMath::Exp([8]*x)",-10,10);
    
    TF1* voigt1 = new TF1("voigt1",voig1,-10,10,6);
    
    TF1* voigt2 = new TF1("voigt2",voig2,-10,10,12);
    TF1* voigt3 = new TF1("voigt3",voig3,-10,10,20);
    TF1* voigt4 = new TF1("voigt4",voig4,-10,10,23);
    TF1* voigt4I = new TF1("voigt4I",voig4,-10,10,23);
    TF1* voigt3_gaussian1 = new TF1("voigt3_gaussian1",lorentz3gauss4,-10,10,19);
    
    
    
    TF1* Sfunc = new TF1("shape", shaper, -10, 10, 7);
    TF1* Sfunc2 = new TF1("shape2", shaper2, -10, 10, 9);
	TH1D *hCQ = (TH1D*)f->Get("hQval");
    TH1D *hBeQ = (TH1D*)fg->Get("hQval");
    TH1D *hBeQ1 = (TH1D*)g1->Get("hQval");
    TH1D *hBeQ2 = (TH1D*)g->Get("hQval");
    TH1D *hBeQ4 = (TH1D*)g4->Get("hQval");
    TH1D *hBeQh1 = (TH1D*)gh1->Get("hQval");
    TH1D *hBeQh2 = (TH1D*)gh2->Get("hQval");
    
    TH1D *h_BeQvalR[16];
    for (int i=0;i<16;i++){
        h_BeQvalR[i] = (TH1D*)g->Get(Form("hQvalR%d",i));
    }
    
    TH1D *h_BeQval2R[8];
    for (int i=0;i<8;i++){
        h_BeQval2R[i] = (TH1D*)g->Get(Form("hQval2R_%d_%d",2*i,2*i+1));
    }
    
    TH1D *h_BeQval3R[5];
    for (int i=0;i<5;i++){
        h_BeQval3R[i] = (TH1D*)g->Get(Form("hQval3R_%d_%d_%d",3*i,3*i+1,3*i+2));
    }
    
    TH1D *h_BeQval4R[4];
    for (int i=0;i<4;i++){
        h_BeQval4R[i] = (TH1D*)g->Get(Form("hQval4R_%d_%d_%d_%d",4*i,4*i+1,4*i+2,4*i+3));
    }
    
    TH1D *h_BeQval2S[4];
    for (int i=0;i<4;i++){
        h_BeQval2S[i] = (TH1D*)fg->Get(Form("hQval%d_%d",90*i,90*(i+1)));
    }
    

    TH2D *h_BeTDL = (TH2D*)g->Get("hYuAnPID");
    
    for (int i=0;i<17;i++){
		angle[i]=180-TMath::RadToDeg()*TMath::ATan((50+(16-i)*(4.94))/TDistance);
		//cout << angle[i] << endl;
	}
    
    

////////////////////////
	////	BERILIUM	////
	////////////////////////
	
	TF1* fit_func1a = new TF1 ( "fit_func1a","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",-4.5,-3.8);
	TF1* fit_func1b = new TF1 ( "fit_func1b","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",-3.35,-2.8);
	TF1* fit_func1c = new TF1 ( "fit_func1c","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",-3,-2);
    
    //Double_t p[32];
    /*TF1* voigt1a = new TF1("voigt1a",voig1,-4.35,-3.95,6);
    voigt1a->SetParameter(0,Qgs13C);
    TF1* voigt1b = new TF1("voigt1b",voig1,-3.4,-2.8,6);
    TF1* voigt1c = new TF1("voigt1c",voig1,-2.9,-2.55,6);
    TF1* voigt1d = new TF1("voigt1d",voig1,-2.45,-2,6);
    */
    //TF1* total = new TF1("total","voigt1(0)+voigt1(8)+voigt1(16)+voigt1(24)",-4.35,-2);
    
    hBeQ->Draw();
	hBeQ->Rebin(8);
    hBeQ->SetLineWidth(2);
    
	// 2-Gaussian Fit for Beryllium
	fit_func2->SetParLimits ( 0,10,50);
	fit_func2->SetParLimits ( 1,-4.3,-4.05);
	fit_func2->SetParLimits ( 3,8,15);
	fit_func2->SetParLimits ( 4,-2.6,-2.4);
	fit_func2->SetParLimits ( 2,0.18,0.2);
    fit_func2->SetParLimits ( 5,0.18,0.2);
    fit_func2->SetParLimits ( 6,-1.1,-0.9);
	
	
	
	// 3-Gaussian Fit for Beryllium
	fit_func3->SetParLimits (0,15,30);
	fit_func3->SetParLimits (1,-4.18,-4.05); //( 1,-4.18,-4.05);
	fit_func3->SetParLimits (3,3,20);
	fit_func3->SetParLimits (4,-3.5,-2.8);
	fit_func3->SetParLimits (5,10,15);
	fit_func3->SetParLimits (6,-2.7,-2.54);//( 6,-2.6,-2.5);
	fit_func3->SetParLimits (2,0.05,0.2);
    
    // 2-Lorentzian Fit for Beryllium
    //lor2->SetParLimits (0,0.1,2);
	lor2->SetParLimits (1,-4.18,-4.05); //( 1,-4.18,-4.05);
	lor2->SetParLimits (2,0.01,0.2);
	//lor2->SetParLimits (3,0.01,1.5);
	lor2->SetParLimits (4,-2.8,-2.4);
    lor2->SetParLimits (5,0.01,0.2);
	
    
    
    // 3-Lorentzian Fit for Beryllium
    //lor3->SetParLimits (0,15,30);
	lor3->SetParLimits (1,-4.18,-4.05); //( 1,-4.18,-4.05);
	//lor3->SetParLimits (3,3,20);
	lor3->SetParLimits (4,-3.4,-2.8);
	lor3->SetParLimits (5,0.01,1.5);
	lor3->SetParLimits (6,-2.58,-2.54);//( 6,-2.6,-2.5);
	lor3->SetParLimits (2,0.05,0.2);
    //lor3->SetParLimits (9,0.05,0.2);
    
    // 2-Voigt Fit for Beryllium
    voigt2->SetParLimits (0,0,10);//Amplitude of first gaussian
    voigt2->SetParLimits (1,-4.18,-4.05); //( 1,-4.18,-4.05);//centroid of first gaussian and lorentzian
	voigt2->SetParLimits (2,0.01,0.2);//SD of all the gaussians
	voigt2->SetParLimits (3,0.01,1.5);//Amplitude of first lorentzian
    voigt2->SetParLimits (4,0.01,0.2);//width of first lorentzian
    voigt2->SetParLimits (5,0.01,0.2);//Amplitude of second gaussian
	voigt2->SetParLimits (6,-2.75,-2.6);//centroid of second gaussian and lorentzian
    //voigt2->SetParLimits (7,0.001,0.05);//Amplitude of second lorentzian
    voigt2->SetParLimits (8,0.01,0.5);//width of second lorentzian
	voigt2->SetParLimits (9,0.01,0.6);//voigt fraction
    
    
    
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
    ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(1000000);
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.001);
    
    fit_func1a->SetParameters(20,-3.92,0.22);

	//fit_func1b->SetParameters(12,-0.36,0.22);
    fit_func1b->SetParLimits (1,-0.6,0.4);
    fit_func1b->SetParameter (2,-3.02);//(2,0.2,0.28);
	
	fit_func1c->SetParameters(10,-2.62,0.22);
    
    
    
    
    /*hBeQ->Fit ( fit_func1a, "R");
	hBeQ->Fit ( fit_func1b, "R+");
	hBeQ->Fit ( fit_func1c, "R+");*/
    
    // 3-Voigt Fit for Beryllium
    voigt3->SetParLimits (0,0.01,15);//Amplitude of first gaussian
    voigt3->SetParLimits (1,-4.25,-4.1); //( 1,-4.18,-4.05);//centroid of first gaussian and lorentzian
	voigt3->SetParLimits (2,0.01,0.26);//SD of all the gaussians
	voigt3->SetParLimits (3,0.01,1);//Amplitude of first lorentzian
    voigt3->SetParLimits (4,0.01,0.3);//width of first lorentzian
    voigt3->SetParLimits (5,0.01,10);//Amplitude of second gaussian
	voigt3->SetParLimits (6,-2.9,-2.55);//centroid of second gaussian and lorentzian
    //voigt3->SetParLimits (7,0.001,0.05);//Amplitude of second lorentzian
    voigt3->SetParLimits (8,0.1,0.5);//width of second lorentzian
    //voigt3->SetParLimits (9,0.01,10);//Amplitude of third gaussian
    voigt3->SetParLimits (10,-2.5,-2);//centroid of third gaussian and lorentzian
    //voigt3->SetParLimits (11,0.01,0.05);//Amplitude of third lorentzian
    voigt3->SetParLimits (12,0.01,0.6);//width of third lorentzian
	voigt3->SetParLimits (13,0.01,0.6);//voigt fraction
    
    // Individual Voigt Fits for Beryllium
    
    /*voigt1a->SetParLimits (1,-4.25,-4.1);
    voigt1a->SetParLimits (2,0.01,0.4);//width of the gaussian
    voigt1a->SetParLimits (4,0.05,0.5);//width of the lorentzian
    voigt1a->SetParLimits (5,0.01,1);//voigt fraction
    
    voigt1b->SetParLimits (1,-3.0,-2.8);
    voigt1b->SetParLimits (2,0.01,0.4);//width of the gaussian
    voigt1b->SetParLimits (4,0.05,0.5);//width of the lorentzian
    voigt1b->SetParLimits (5,0.01,1);//voigt fraction
    
    voigt1c->SetParLimits (1,-2.8,-2.6);
    voigt1c->SetParLimits (2,0.01,0.4);//width of the gaussian
    voigt1c->SetParLimits (4,0.01,0.5);//width of the lorentzian
    voigt1c->SetParLimits (5,0.01,1);//voigt fraction
    
    voigt1d->SetParLimits (1,-2.4,-2.25);
    voigt1d->SetParLimits (2,0.01,0.3);//width of the gaussian
    voigt1d->SetParLimits (4,0.05,0.5);//width of the lorentzian
    voigt1d->SetParLimits (5,0.01,1);//voigt fraction
    
    hBeQ->Fit ( voigt1a, "R");
	hBeQ->Fit ( voigt1b, "R+");
	hBeQ->Fit ( voigt1c, "R+");
    hBeQ->Fit ( voigt1d, "R+");*/
    
    /*Double_t p[26];
    voigt1a->GetParameters(&p[0]);
    voigt1b->GetParameters(&p[6]);
    voigt1c->GetParameters(&p[12]);
    voigt1d->GetParameters(&p[18]);
    
    voigt4->SetParameters(p);*/
    
    // 4-Voigt Fit for Beryllium
    voigt4->SetParLimits (0,0.01,15);//Amplitude of first gaussian
    voigt4->SetParLimits (1,-4.3,-4.13); //( 1,-4.18,-4.05);//centroid of first gaussian and lorentzian
	voigt4->SetParLimits (2,0.01,0.26);//SD of all the gaussians
	//voigt4->SetParLimits (3,0.01,1);//Amplitude of first lorentzian
    voigt4->SetParLimits (4,0.01,0.3);//width of first lorentzian
    voigt4->SetParLimits (5,0.01,10);//Amplitude of second gaussian
	voigt4->SetParLimits (6,-3.2,-2.95);//centroid of second gaussian and lorentzian
    voigt4->SetParLimits (7,0,0.05);//Amplitude of second lorentzian
    voigt4->SetParLimits (8,0.1,0.5);//width of second lorentzian
    voigt4->SetParLimits (9,0,50);//Amplitude of third gaussian
    voigt4->SetParLimits (10,-2.7,-2.6);//centroid of third gaussian and lorentzian
    voigt4->SetParLimits (11,0,0.1);//Amplitude of third lorentzian
    voigt4->SetParLimits (12,0.01,1);//width of third lorentzian
	//voigt4->SetParLimits (13,0,0.05);//Amplitude of fourth gaussian
    voigt4->SetParLimits (14,-2.4,-2.2);//centroid of fourth gaussian and lorentzian
    voigt4->SetParLimits (15,0,0.5);//Amplitude of fourth lorentzian
    voigt4->SetParLimits (16,0.01,0.7);//width of fourth lorentzian
    voigt4->SetParLimits (17,0.01,0.6);//voigt1 fraction
    voigt4->SetParLimits (18,0.01,0.9);//voigt2 fraction
    voigt4->SetParLimits (19,0.01,0.9);//voigt3 fraction
    voigt4->SetParLimits (20,0.01,0.9);//voigt4 fraction
	voigt4->SetParLimits (21,0.01,0.2);
    
    
    
    
	//hBeQ->Fit ( "total","R+");
    hBeQ->Fit ( "voigt4","+","",-4.35,-2.15);
    hBeQ->GetFunction("voigt4")->SetLineColor(kBlack);
    hBeQ->GetFunction("voigt4")->SetLineWidth(3);
    
    
    
    
    TF1* voigt1a = new TF1("voigt1a",voig1,-4.35,-3.9,6);
    voigt1a->SetParameters(voigt4->GetParameter(0),voigt4->GetParameter(1),voigt4->GetParameter(2),voigt4->GetParameter(3),voigt4->GetParameter(4),voigt4->GetParameter(17));
    TF1* voigt1b = new TF1("voigt1b",voig1,-3.2,-2.8,6);
    voigt1b->SetParameters(voigt4->GetParameter(5),voigt4->GetParameter(6),voigt4->GetParameter(2),voigt4->GetParameter(7),voigt4->GetParameter(8),voigt4->GetParameter(18));
    TF1* voigt1c = new TF1("voigt1c",voig1,-2.9,-2.55,6);
    voigt1c->SetParameters(voigt4->GetParameter(9),voigt4->GetParameter(10),voigt4->GetParameter(2),voigt4->GetParameter(11),voigt4->GetParameter(12),voigt4->GetParameter(19));
    TF1* voigt1d = new TF1("voigt1d",voig1,-2.5,-2.15,6);
    voigt1d->SetParameters(voigt4->GetParameter(13),voigt4->GetParameter(14),voigt4->GetParameter(21),voigt4->GetParameter(15),voigt4->GetParameter(16),voigt4->GetParameter(20));
    TF1* background = new TF1("background",poly6_background,-4.35,-2.15,2);
    background->SetParameter(0,voigt4->GetParameter(22));
    
    voigt1a->Draw("same");
    voigt1a->SetLineWidth(3);
    voigt1b->Draw("same");
    voigt1b->SetLineWidth(3);
    voigt1c->Draw("same");
    voigt1c->SetLineWidth(3);
    voigt1d->Draw("same");
    voigt1d->SetLineWidth(3);
    background->Draw("same");
    background->SetLineWidth(3);
    background->SetLineColor(3);
    
    /*voigt3_gaussian1->SetParLimits (0,0.01,15);//Amplitude of first gaussian
    voigt3_gaussian1->SetParLimits (1,-4.3,-4.13); //( 1,-4.18,-4.05);//centroid of first gaussian and lorentzian
	voigt3_gaussian1->SetParLimits (2,0.01,0.3);//SD of all the gaussians
	//voigt3_gaussian1->SetParLimits (3,0.01,1);//Amplitude of first lorentzian
    voigt3_gaussian1->SetParLimits (4,0.01,0.4);//width of first lorentzian
    voigt3_gaussian1->SetParLimits (5,0.01,20);//Amplitude of second gaussian
	voigt3_gaussian1->SetParLimits (6,-3.3,-2.9);//centroid of second gaussian
    voigt3_gaussian1->SetParLimits (7,0,50);//Amplitude of third gaussian
    voigt3_gaussian1->SetParLimits (8,-2.75,-2.6);//centroid of third gaussian and lorentzian
    voigt3_gaussian1->SetParLimits (9,0,0.7);//Amplitude of third lorentzian
    voigt3_gaussian1->SetParLimits (10,0.01,5);//width of third lorentzian
	voigt3_gaussian1->SetParLimits (11,0,0.05);//Amplitude of fourth gaussian
    voigt3_gaussian1->SetParLimits (12,-2.4,-2.2);//centroid of fourth gaussian and lorentzian
    voigt3_gaussian1->SetParLimits (13,0.01,0.2);//SD of fourth gaussian
    //voigt3_gaussian1->SetParLimits (14,0,0.5);//Amplitude of fourth lorentzian
    voigt3_gaussian1->SetParLimits (15,0.01,0.6);//width of fourth lorentzian
    voigt3_gaussian1->SetParLimits (16,0.01,0.6);//voigt fraction
    
    hBeQ->Fit ( "voigt3_gaussian1","+","",-4.35,-2.15);
    hBeQ->GetFunction("voigt3_gaussian1")->SetLineColor(kBlack);
    hBeQ->GetFunction("voigt3_gaussian1")->SetLineWidth(3);
    
    TF1* voigt1a = new TF1("voigt1a",voig1,-4.35,-3.9,6);
    voigt1a->SetParameters(voigt3_gaussian1->GetParameter(0),voigt3_gaussian1->GetParameter(1),voigt3_gaussian1->GetParameter(2),voigt3_gaussian1->GetParameter(3),voigt3_gaussian1->GetParameter(4),voigt3_gaussian1->GetParameter(16));
    TF1* gauss1b = new TF1("gauss1b",gauss1,-3.2,-2.8,6);
    gauss1b->SetParameters(voigt3_gaussian1->GetParameter(5),voigt3_gaussian1->GetParameter(6),voigt3_gaussian1->GetParameter(2));
    TF1* voigt1c = new TF1("voigt1c",voig1,-2.9,-2.55,6);
    voigt1c->SetParameters(voigt3_gaussian1->GetParameter(7),voigt3_gaussian1->GetParameter(8),voigt3_gaussian1->GetParameter(2),voigt3_gaussian1->GetParameter(9),voigt3_gaussian1->GetParameter(10),voigt3_gaussian1->GetParameter(16));
    TF1* voigt1d = new TF1("voigt1d",voig1,-2.5,-2.15,6);
    voigt1d->SetParameters(voigt3_gaussian1->GetParameter(11),voigt3_gaussian1->GetParameter(12),voigt3_gaussian1->GetParameter(13),voigt3_gaussian1->GetParameter(14),voigt3_gaussian1->GetParameter(15),voigt3_gaussian1->GetParameter(16));
    
    voigt1a->Draw("same");
    voigt1a->SetLineWidth(3);
    gauss1b->Draw("same");
    gauss1b->SetLineWidth(3);
    voigt1c->Draw("same");
    voigt1c->SetLineWidth(3);
    voigt1d->Draw("same");
    voigt1d->SetLineWidth(3);
    */
    
	hBeQ->GetXaxis()->SetRangeUser(-5,-2);
    hBeQ->SetTitle("");
    hBeQ->GetXaxis()->SetTitle("Q Value [MeV]");
    hBeQ->GetYaxis()->SetTitle("Counts");
    hBeQ->GetXaxis()->SetLabelSize(0.05);
    hBeQ->GetXaxis()->SetTitleSize(0.05);
    hBeQ->GetXaxis()->SetTitleOffset(0.9);
    hBeQ->GetYaxis()->SetLabelSize(0.05);
    hBeQ->GetYaxis()->SetTitleSize(0.05);
    hBeQ->GetYaxis()->SetTitleOffset(0.9);
    
    //TCanvas *c4 = new TCanvas ( "c4" );
    h_BeTDL->SetTitle("");
	//h_BeTDL->Draw("colz");
    h_BeTDL->GetYaxis()->SetRangeUser(0.6,2.4);
    h_BeTDL->GetXaxis()->SetTitle("Angle [degrees]");
    h_BeTDL->GetYaxis()->SetTitle("Energy of the protons [MeV]");
    h_BeTDL->GetXaxis()->SetLabelSize(0.05);
    h_BeTDL->GetXaxis()->SetTitleSize(0.05);
    h_BeTDL->GetXaxis()->SetTitleOffset(0.9);
    h_BeTDL->GetYaxis()->SetLabelSize(0.05);
    h_BeTDL->GetYaxis()->SetTitleSize(0.05);
    h_BeTDL->GetYaxis()->SetTitleOffset(0.9);
    TF1 *KinFcn = new TF1("KinFcn",kin,angle[0],angle[16],1);
    float Q13Be0 = -2.2202;
	KinFcn->SetParameter(0,Q13Be0);
    KinFcn->SetLineColor(kRed);
	KinFcn->SetLineWidth(3);
	//KinFcn->Draw("same");
    
    Double_t p[23];
    voigt4->GetParameters(&p[0]);
    //voigt4I->SetParameters(p);
    
    // 4-Voigt Fit for Beryllium
    voigt4I->SetParLimits (0,0.01,15);//Amplitude of first gaussian
    voigt4I->FixParameter (1,p[1]); //( 1,-4.18,-4.05);//centroid of first gaussian and lorentzian
	voigt4I->FixParameter (2,p[2]);//SD of all the gaussians
	//voigt4I->SetParLimits (3,0.01,1);//Amplitude of first lorentzian
    voigt4I->FixParameter (4,p[4]);//width of first lorentzian
    //voigt4I->SetParLimits (5,0.01,10);//Amplitude of second gaussian
	voigt4I->FixParameter (6,p[6]);//centroid of second gaussian and lorentzian
    //voigt4I->SetParLimits (7,0,0.05);//Amplitude of second lorentzian
    voigt4I->FixParameter (8,p[8]);//width of second lorentzian
    //voigt4I->SetParLimits (9,0,50);//Amplitude of third gaussian
    voigt4I->FixParameter (10,p[10]);//centroid of third gaussian and lorentzian
    //voigt4I->SetParLimits (11,0,0.1);//Amplitude of third lorentzian
    voigt4I->FixParameter (12,p[12]);//width of third lorentzian
	//voigt4I->SetParLimits (13,0,10);//Amplitude of fourth gaussian
    voigt4I->FixParameter (14,p[14]);//centroid of fourth gaussian and lorentzian
    //voigt4I->SetParLimits (15,0,0.5);//Amplitude of fourth lorentzian
    voigt4I->FixParameter (16,p[16]);//width of fourth lorentzian
    voigt4I->FixParameter (17,p[17]);//voigt1 fraction
    voigt4I->FixParameter (18,p[18]);//voigt2 fraction
    voigt4I->FixParameter (19,p[19]);//voigt3 fraction
    voigt4I->FixParameter (20,p[20]);//voigt4 fraction
    voigt4I->SetParLimits (22,-1.5,0); //Background amplitude
	voigt4I->FixParameter (21,p[21]); //SD of fourth gaussian
    
    
    
   /* TCanvas *c2 = new TCanvas ( "c2" );
    c2->Divide(4,4);
    for (int i=0;i<16;i++){
        c2->cd(i+1);
        h_BeQvalR[i]->Draw();
        h_BeQvalR[i]->Rebin(4);
        h_BeQvalR[i]->Fit ( "voigt4I","+","",-4.35,-2.15);
        h_BeQvalR[i]->GetXaxis()->SetRangeUser(-5,-2);
    }
    
    TCanvas *c3 = new TCanvas ( "c3" );
    c3->Divide(3,3);
    for (int i=0;i<8;i++){
        c3->cd(i+1);
        h_BeQval2R[i]->Draw();
        h_BeQval2R[i]->Rebin(4);
        h_BeQval2R[i]->Fit ( "voigt4I","+","",-4.35,-2.15);
        h_BeQval2R[i]->GetXaxis()->SetRangeUser(-5,-2);
    }
    
    TCanvas *c4 = new TCanvas ( "c4" );
    c4->Divide(2,3);
    for (int i=0;i<5;i++){
        c4->cd(i+1);
        h_BeQval3R[i]->Draw();
        h_BeQval3R[i]->Rebin(4);
        h_BeQval3R[i]->Fit ( "voigt4I","+","",-4.35,-2.15);
        h_BeQval3R[i]->GetXaxis()->SetRangeUser(-5,-2);
    }
    */
   
   ///////////       CALCULATION OF SOLID ANGLES START          //////////////
   
    //read in the geometry file to calculate the angles
	ifstream geometry;
	geometry.open ("/home/jerome/12Be_exp/scripts/geometry_s1506.txt"); //open the geometry file; change the path and name accordingly
	if(!geometry.is_open()){
		cout << " No Geometry file found " << endl;
		return -1;
	}

	string read_geometry;
	istringstream iss;
	string name, dummy;
	float Ydr0, Ydr1, Yuz; // Ydz distance from the target, Ydr0 inner radius, Ydz outer radius, same for the Sd

	while ( getline ( geometry,read_geometry ) ){ //in the calibration file named geometry start reading the lines

		if ( read_geometry.find ( "YD_INNER_RADIUS",0 ) !=string::npos ){
			iss.clear();
			iss.str ( read_geometry );
			iss >> name >> dummy >> Ydr0;
		}

		if ( read_geometry.find ( "YD_OUTER_RADIUS",0 ) !=string::npos ){
			iss.clear();
			iss.str ( read_geometry );
			iss >> name >> dummy >> Ydr1;
		}

		if ( read_geometry.find ( "YU_DISTANCE",0 ) !=string::npos ){
			iss.clear();
			iss.str ( read_geometry );
			iss >> name >> dummy >> Yuz;
		}
	}
	//end of the geometry file
	
	double width=(Ydr1-Ydr0)/16;
    
    int n=4;
    int distLength=(16/n)+1, omegaLength=(16/n);
    double dist[distLength], omega[omegaLength];
    double angle1[omegaLength], solidAngleBins[distLength];
    

    
    for (int i=0;i<distLength;i++){
        dist[i]=Ydr0+n*i*width;
        solidAngleBins[i]=180-TMath::RadToDeg()*TMath::ATan((Ydr0+(distLength-1-i)*n*width)/(0-Yuz));
        //cout << solidAngleBins[i] << endl;
    }
    
    for (int i=0;i<omegaLength;i++){
        omega[i]=42*(M_PI/180)*(0-Yuz)*((1/TMath::Sqrt(pow(dist[i],2)+Yuz*Yuz))-(1/TMath::Sqrt(pow(dist[i+1],2)+Yuz*Yuz)));
        angle1[i]=(solidAngleBins[i]+solidAngleBins[i+1])/2;
        //cout << omega[i] << endl;
    }
    
    ///////////////          CALCULATION OF SOLID ANGLES End            ////////////////////
    
    int counts[4][4] = {0}, counts_noB[4][4] = {0}, totalCounts[4][4] = {0}, totalCounts_noB[4][4] = {0}; //first index for angles and second index for the states
    double eValue=0;
    double energyCut[4][2]={{-4.3,-4},{-3.1,-2.85},{-2.8,-2.45},{-2.48,-2.2}};
    double energyCutBin[4][2]={0}, error[4][4];

    
    for (int i=0;i<4;i++){
        h_BeQval4R[i]->Rebin(6);
    }
    for (int i=0;i<4;i++){
        energyCutBin[i][0]=h_BeQval4R[0]->GetXaxis()->FindBin(energyCut[i][0]);
        energyCutBin[i][1]=h_BeQval4R[0]->GetXaxis()->FindBin(energyCut[i][1]);
    }
    
    
    TCanvas *c5 = new TCanvas ( "c5" );
    c5->Divide(2,2);
    for (int i=0;i<4;i++){ // 4 angles
        c5->cd(i+1);
        h_BeQval4R[i]->Draw();
        h_BeQval4R[i]->Fit ( "voigt4I","+MQ","",-4.35,-2.15);
        h_BeQval4R[i]->GetXaxis()->SetRangeUser(-5,-2);
        
        /*background->SetParameters(voigt4I->GetParameter(21),voigt4I->GetParameter(22));
        background->Draw("same");
        background->SetLineWidth(3);
        background->SetLineColor(3);
        */
        
        //cout << voigt4I->GetParameter(21) << "  " << background->GetParameter(0) << "  " << voigt4I->GetParameter(22) << " " << background->GetParameter(1) << endl;
        
        for (int k=0;k<4;k++){ // 4 states
            for (int j=energyCutBin[k][0];j<energyCutBin[k][1];j++){ //looping through the bins of each state
                eValue=h_BeQval4R[i]->GetXaxis()->GetBinCenter(j);
                double a=voigt4I->GetParameter(22);
                //double b=voigt4I->GetParameter(22);
                
                counts[i][k]=(h_BeQval4R[i]->GetBinContent(j));
                counts_noB[i][k]=(h_BeQval4R[i]->GetBinContent(j))-(a*pow(eValue,6)); //subtracting the bakground
                
                totalCounts[i][k]=totalCounts[i][k]+counts[i][k];
                totalCounts_noB[i][k]=totalCounts_noB[i][k]+counts_noB[i][k]; //No background
            }
            error[i][k]=sqrt(totalCounts[i][k]);
            //cout << totalCounts[i][k] << endl;
            //cout << totalCounts_noB[i][k] << endl;
        }
    }
    
    TGraphErrors *gr1[4];
    TCanvas *c6 = new TCanvas ( "c6" );
    c6->Divide(2,2);
    for (int i=0;i<4;i++){ // 4 states
        gr1[i] = new TGraphErrors();
        for (int k=0;k<4;k++){ // 4 angles
            gr1[i]->SetPoint(k,angle1[k],totalCounts[k][i]/omega[k]);
            gr1[i]->SetPointError(k,0,error[k][i]/omega[k]);
        }
        c6->cd(i+1);
        gr1[i]->Draw("AP*");
        gr1[i]->SetTitle(Form("State %d",i+1));
    }
    
    TCanvas *c7 = new TCanvas ( "c7" );
    c7->Divide(2,2);
    for (int i=0;i<4;i++){ // 4 angle ranges in phi
        c7->cd(i+1);
        h_BeQval2S[i]->Rebin(6);
        h_BeQval2S[i]->Draw();
        h_BeQval2S[i]->GetXaxis()->SetRangeUser(-5,-2);
    }
    
    TCanvas *c8 = new TCanvas ( "c8" );
    c8->Divide(2,2);
    c8->cd(1);
    //hBeQ->Rebin(6);
    hBeQ->Draw();
    hBeQ->GetXaxis()->SetRangeUser(-5,-2);
        
    c8->cd(2);
    hBeQ1->Rebin(8);
    hBeQ1->Draw();
    hBeQ1->GetXaxis()->SetRangeUser(-5,-2);
     
    c8->cd(3);
    hBeQ2->Rebin(8);
    hBeQ2->Draw();
    hBeQ2->GetXaxis()->SetRangeUser(-5,-2);
        
    c8->cd(4);
    hBeQ4->Rebin(8);
    hBeQ4->Draw();
    hBeQ4->GetXaxis()->SetRangeUser(-5,-2);
        
    TCanvas *c9 = new TCanvas ( "c9" );
    c9->Divide(2,1);
    c9->cd(1);
    hBeQh1->Rebin(8);
    hBeQh1->Draw();
    hBeQh1->GetXaxis()->SetRangeUser(-5,-2);
    
    c9->cd(2);
    hBeQh2->Rebin(8);
    hBeQh2->Draw();
    hBeQh2->GetXaxis()->SetRangeUser(-5,-2);
    
}
