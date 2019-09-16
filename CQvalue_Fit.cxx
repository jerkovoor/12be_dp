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
    return p[17]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[17])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[18]*(p[5]*exp(-pow((x[0]-p[6]),2)/(2*p[2]*p[2])))+(1-p[18])*(p[7]/(pow((x[0]-p[6]),2)+pow((p[8]/2),2)))+p[19]*(p[9]*exp(-pow((x[0]-p[10]),2)/(2*p[2]*p[2])))+(1-p[19])*(p[11]/(pow((x[0]-p[10]),2)+pow((p[12]/2),2)))+p[20]*(p[13]*exp(-pow((x[0]-p[14]),2)/(2*p[23]*p[23])))+(1-p[20])*(p[15]/(pow((x[0]-p[14]),2)+pow((p[16]/2),2)))+p[21]*TMath::Exp(p[22]*x[0]);
}


double lorentz3gauss4(double *x, double *p) {
    return p[16]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[16])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[16]*(p[5]*exp(-pow((x[0]-p[6]),2)/(2*p[2]*p[2])))+p[16]*(p[7]*exp(-pow((x[0]-p[8]),2)/(2*p[2]*p[2])))+(1-p[16])*(p[9]/(pow((x[0]-p[8]),2)+pow((p[10]/2),2)))+p[16]*(p[11]*exp(-pow((x[0]-p[12]),2)/(2*p[13]*p[13])))+(1-p[16])*(p[14]/(pow((x[0]-p[12]),2)+pow((p[15]/2),2)))+p[17]*TMath::Exp(p[18]*x[0]);
}
/*double voig4(double *x, double *p) {
    return p[5]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[5])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[11]*(p[6]*exp(-pow((x[0]-p[7]),2)/(2*p[8]*p[8])))+(1-p[11])*(p[9]/(pow((x[0]-p[7]),2)+pow((p[10]/2),2)))+p[17]*(p[12]*exp(-pow((x[0]-p[13]),2)/(2*p[14]*p[14])))+(1-p[17])*(p[15]/(pow((x[0]-p[13]),2)+pow((p[16]/2),2)))+p[23]*(p[18]*exp(-pow((x[0]-p[19]),2)/(2*p[20]*p[20])))+(1-p[23])*(p[21]/(pow((x[0]-p[19]),2)+pow((p[22]/2),2)))+p[24]*TMath::Exp(p[25]*x[0]);
}*/


void CQvalue_Fit(){
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
    TF1* voigt4 = new TF1("voigt4",voig4,-10,10,24);
    TF1* voigt4I = new TF1("voigt4I",voig4,-10,10,24);
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
	////	CARBON		////
	////////////////////////
    
    TF1* fit_func1d = new TF1 ( "fit_func1d","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",-1.65,-1.2);
    TF1* fit_func1a = new TF1 ( "fit_func1a","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",-1.2,-0.5);
	TF1* fit_func1b = new TF1 ( "fit_func1b","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",-0.6,0.4);
	TF1* fit_func1c = new TF1 ( "fit_func1c","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",2,3.5);
    
	hCQ->Draw ();
	//hCQ->Rebin(4);
	
    

	// 3-Gaussian Fit for Carbon
	fit_func3->SetParLimits ( 0,25,200);
	fit_func3->SetParLimits ( 1,-1.5,-0.8);
	fit_func3->SetParLimits ( 3,5,40);
	fit_func3->SetParameter (4,-0.36);//( 4,-0.6,0.4);
	fit_func3->SetParLimits ( 5,5,50);
	fit_func3->SetParLimits ( 6,2.1,3.2);
	fit_func3->SetParLimits ( 2,0.2,0.28);
    
    fit_func4->SetParLimits ( 0,18,30);
    fit_func4->SetParLimits ( 1,-1.5,-1);
	//fit_func4->FixParameter ( 1,-1.3);
    fit_func4->SetParLimits ( 3,25,35);
	fit_func4->SetParLimits ( 4,-1.5,-0.8);
    //fit_func4->FixParameter ( 4,-0.98);
	fit_func4->SetParLimits ( 5,5,40);
	fit_func4->SetParameter (6,-0.36);//( 4,-0.6,0.4);
	fit_func4->SetParLimits ( 7,5,50);
	fit_func4->SetParLimits ( 8,2.1,3.2);
	fit_func4->SetParLimits ( 2,0.2,0.28);
	
	
    
    Sfunc->SetParameters(1, -1.1, 0.5, 1, -0.3, 1, 2.72);
    Sfunc->SetParLimits(1, -1.3, -1.0);
    Sfunc->SetParLimits(4, -1.3, -0.9);
    Sfunc->SetParLimits(6, 2.6, 2.8);
    
    Sfunc2->SetParameters(1, -1.1, 0.5, 0.2, -0.9, 1, -0.3, 1, 2.72);
    Sfunc2->SetParLimits(1, -1.4, -1.0);
    Sfunc2->SetParLimits(4, -1.2, -0.8);
    Sfunc2->SetParLimits(6, -0.5, 0.);
    Sfunc2->SetParLimits(8, 2.6, 2.8);
    
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
    ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(1000000);
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.001);
    
    hCQ->Fit ( "fit_func4","R","",-3,4);
    hCQ->GetFunction("fit_func4")->SetLineColor(kRed);
    Double_t p[11];
    fit_func4->GetParameters(&p[0]);
    
    fit_func1a->FixParameter(0,p[3]);
    fit_func1a->FixParameter(1,p[4]);
    fit_func1a->FixParameter(2,p[2]);
    
    fit_func1b->FixParameter(0,p[5]);
    fit_func1b->FixParameter(1,p[6]);
    fit_func1b->FixParameter(2,p[2]);
    
    fit_func1c->FixParameter(0,p[7]);
    fit_func1c->FixParameter(1,p[8]);
    fit_func1c->FixParameter(2,p[2]);
    
    fit_func1d->FixParameter(0,p[0]);
    fit_func1d->FixParameter(1,p[1]);
    fit_func1d->FixParameter(2,p[2]);

		
	hCQ->Fit ( fit_func1a, "R+");
	hCQ->Fit ( fit_func1b, "R+");
	hCQ->Fit ( fit_func1c, "R+");
    hCQ->Fit ( fit_func1d, "R+");
	
    hCQ->GetFunction("fit_func1a")->SetLineColor(kBlack);
    hCQ->GetFunction("fit_func1b")->SetLineColor(kBlack);
    hCQ->GetFunction("fit_func1c")->SetLineColor(kBlack);
    hCQ->GetFunction("fit_func1d")->SetLineColor(kBlack);
    
    hCQ->SetTitle("12C(d,p)13C Q value, TRIUMF DL");
    /*TGraph *gr2 = new TGraph(2,QC0,countRange);
    gr2->Draw("same");
    gr2->SetLineColor(7);
    gr2->SetLineWidth(3);
    TGraph *gr3 = new TGraph(2,QC1,countRange);
    gr3->Draw("same");
    gr3->SetLineColor(kBlack);
    gr3->SetLineWidth(3);
    TGraph *gr22 = new TGraph(2,QC2,countRange);
    gr22->Draw("same");
    gr22->SetLineColor(4);
    gr22->SetLineWidth(3);
    TGraph *gr4 = new TGraph(2,QC3,countRange);
    gr4->Draw("same");
    gr4->SetLineColor(6);
    gr4->SetLineWidth(3);
    //hCQ->Fit("shape", "R+");
    //hCQ->Fit("shape2", "R+");
	hCQ->GetXaxis()->SetRangeUser(-3,4);
	*/
    
}
