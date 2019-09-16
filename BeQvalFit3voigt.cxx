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

double shaper(double *x, double *p) {
    return p[0]*exp(-(x[0] - p[1])*(x[0] - p[1])/(2.*p[2]*p[2])) + p[3]*exp(-(x[0] - p[4])*(x[0] - p[4])/(2.*p[2]*p[2])) + p[5]*exp(-(x[0] - p[6])*(x[0] - p[6])/(2.*p[2]*p[2]));
}

double shaper2(double *x, double *p) {
    return p[0]*exp(-(x[0] - p[1])*(x[0] - p[1])/(2.*p[2]*p[2])) + p[3]*exp(-(x[0] - p[4])*(x[0] - p[4])/(2.*p[2]*p[2])) + p[5]*exp(-(x[0] - p[6])*(x[0] - p[6])/(2.*p[2]*p[2])) + 
        + p[7]*exp(-(x[0] - p[8])*(x[0] - p[8])/(2.*p[2]*p[2]));
}

double var_width(double *x, double *p) {
    return 2*p[0]/(1+exp(p[1]*(x[0]-p[2])));
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

double voig3(double *x, double *p) {
    return p[13]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[13])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[13]*(p[5]*exp(-pow((x[0]-p[6]),2)/(2*p[2]*p[2])))+(1-p[13])*(p[7]/(pow((x[0]-p[6]),2)+pow((p[8]/2),2)))+p[13]*(p[9]*exp(-pow((x[0]-p[10]),2)/(2*p[2]*p[2])))+(1-p[13])*(p[11]/(pow((x[0]-p[10]),2)+pow((p[12]/2),2)))+p[14]*TMath::Exp(p[15]*x[0]);
}

/*double voig3(double *x, double *p) {
    return p[5]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[5])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[11]*(p[6]*exp(-pow((x[0]-p[7]),2)/(2*p[8]*p[8])))+(1-p[11])*(p[9]/(pow((x[0]-p[7]),2)+pow((p[10]/2),2)))+p[17]*(p[12]*exp(-pow((x[0]-p[13]),2)/(2*p[14]*p[14])))+(1-p[17])*(p[15]/(pow((x[0]-p[13]),2)+pow((p[16]/2),2)))+p[18]*TMath::Exp(p[19]*x[0]);
}*/

double voig4(double *x, double *p) {
    return p[17]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[17])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[17]*(p[5]*exp(-pow((x[0]-p[6]),2)/(2*p[2]*p[2])))+(1-p[17])*(p[7]/(pow((x[0]-p[6]),2)+pow((p[8]/2),2)))+p[17]*(p[9]*exp(-pow((x[0]-p[10]),2)/(2*p[2]*p[2])))+(1-p[17])*(p[11]/(pow((x[0]-p[10]),2)+pow((p[12]/2),2)))+p[17]*(p[13]*exp(-pow((x[0]-p[14]),2)/(2*p[2]*p[2])))+(1-p[17])*(p[15]/(pow((x[0]-p[14]),2)+pow((p[16]/2),2)))+p[18]*TMath::Exp(p[19]*x[0]);
}

/*double voig4(double *x, double *p) {
    return p[5]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[5])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[11]*(p[6]*exp(-pow((x[0]-p[7]),2)/(2*p[8]*p[9])))+(1-p[11])*(p[9]/(pow((x[0]-p[7]),2)+pow((p[10]/2),2)))+p[17]*(p[12]*exp(-pow((x[0]-p[13]),2)/(2*p[14]*p[14])))+(1-p[17])*(p[15]/(pow((x[0]-p[13]),2)+pow((p[16]/2),2)))+p[23]*(p[18]*exp(-pow((x[0]-p[19]),2)/(2*p[20]*p[20])))+(1-p[23])*(p[21]/(pow((x[0]-p[19]),2)+pow((p[22]/2),2)))+p[24]*TMath::Exp(p[25]*x[0]);
}*/


void BeQvalFit3voigt(){
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
    
    TFile *f = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/C_pedestal_TRIUMF_DL_NoOffset_Shift0_0883_Yu_TargetDistance80_87mm.root","READ");
	//TFile *f = new TFile("/home/jerome/12Be_exp/Analysis/CarbonGain/C_pedestal_TRIUMF_DL_QvalMinGS_Yu.root","READ");
	TFile *g = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod.root","READ");

	TF1* fit_func2 = new TF1 ( "fit_func2","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp([6]*x)",-10,10);
    TF1* fit_func3 = new TF1 ( "fit_func3","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp([8]*x)",-10,10);
	TF1* fit_func4 = new TF1 ( "fit_func4","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp([10]*x)",-10,10);
    TF1* lor2 = new TF1 ("lor2","([0]*[2]*[2])/(pow((x-[1]),2)+pow(([2]/2),2))+([3]*[2]*[2])/(pow((x-[4]),2)+pow(([5]/2),2))+[6]*TMath::Exp([7]*x)",-10,10);
    TF1* lor3 = new TF1 ("lor3","([0]*[2]*[2])/(pow((x-[1]),2)+pow(([2]/2),2))+([3]*[2]*[2])/(pow((x-[4]),2)+pow(([2]/2),2))+([5]*[2]*[2])/(pow((x-[6]),2)+pow(([2]/2),2))+[7]*TMath::Exp([8]*x)",-10,10);
    
    TF1* voigt1 = new TF1("voigt1","[5]*([0]*exp(-pow((x[0]-[1]),2)/(2*[2]*[2])))+(1-[5])*([3]/(pow((x[0]-[1]),2)+pow(([4]/2),2)))+[6]*TMath::Exp([7]*x[0])",-10,10);
    
    TF1* voigt2 = new TF1("voigt2",voig2,-10,10,12);
    TF1* voigt3 = new TF1("voigt3",voig3,-10,10,16); //20
    TF1* voigt4 = new TF1("voigt4",voig4,-10,10,26);
    
    TF1* Sfunc = new TF1("shape", shaper, -10, 10, 7);
    TF1* Sfunc2 = new TF1("shape2", shaper2, -10, 10, 9);
	TH1D *hCQ = (TH1D*)f->Get("hQval");
    TH1D *hBeQ = (TH1D*)g->Get("hQval");
    
    
    ////////////////////////
	////	CARBON		////
	////////////////////////
    
    /*TF1* fit_func1a = new TF1 ( "fit_func1a","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",-2,-0.5);
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
	
	fit_func1a->SetParameters(35,-1.1,0.22);

	//fit_func1b->SetParameters(12,-0.36,0.22);
    fit_func1b->SetParLimits (1,-0.6,0.4);
    fit_func1b->SetParameter (2,-0.36);//(2,0.2,0.28);
	
	fit_func1c->SetParameters(15,2.7,0.22);
    
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
		
	hCQ->Fit ( fit_func1a, "R");
	hCQ->Fit ( fit_func1b, "R+");
	hCQ->Fit ( fit_func1c, "R+");
	hCQ->Fit ( "fit_func3","R+","",-3,4);
    hCQ->GetFunction("fit_func3")->SetLineColor(kGreen);
    hCQ->SetTitle("12C(d,p)13C Q value, TRIUMF DL");
    TGraph *gr2 = new TGraph(2,QC0,countRange);
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
	hCQ->GetXaxis()->SetRangeUser(-3,4);*/
	
    
    ////////////////////////
	////	BERILIUM	////
	////////////////////////
	
	TF1* fit_func1a = new TF1 ( "fit_func1a","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",-4.5,-3.8);
	TF1* fit_func1b = new TF1 ( "fit_func1b","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",-3.35,-2.8);
	TF1* fit_func1c = new TF1 ( "fit_func1c","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",-3,-2);
    
    //Double_t p[32];
    /*TF1* voigt1a = new TF1("voigt1a",voig1,-4.3,-4,6);
    TF1* voigt1b = new TF1("voigt1b",voig1,-3.2,-2.82,6);
    TF1* voigt1c = new TF1("voigt1c",voig1,-2.75,-2.2,6);*/
    
    //TF1* voigt1c = new TF1("voigt1c",voig1,-2.9,-2.45,6);
    //TF1* voigt1d = new TF1("voigt1d",voig1,-2.45,-2,6);
    
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
    voigt3->SetParLimits (1,-4.25,-4.15); //( 1,-4.18,-4.05);//centroid of first gaussian and lorentzian
	voigt3->SetParLimits (2,0.01,0.26);//SD of all the gaussians
	voigt3->SetParLimits (3,0.01,1);//Amplitude of first lorentzian
    voigt3->SetParLimits (4,0.01,0.25);//width of first lorentzian
    voigt3->SetParLimits (5,0.01,10);//Amplitude of second gaussian
	voigt3->SetParLimits (6,-3.1,-2.8);//centroid of second gaussian and lorentzian
    voigt3->SetParLimits (7,0.001,0.05);//Amplitude of second lorentzian
    voigt3->SetParLimits (8,0.1,0.5);//width of second lorentzian
    //voigt3->SetParLimits (9,0.01,10);//Amplitude of third gaussian
    voigt3->SetParLimits (10,-2.9,-2.55);//centroid of third gaussian and lorentzian
    //voigt3->SetParLimits (11,0.001,0.05);//Amplitude of third lorentzian
    voigt3->SetParLimits (12,0.01,0.5);//width of third lorentzian
	voigt3->SetParLimits (13,0.01,0.5);//voigt fraction
    
    
    // Individual Voigt Fits for Beryllium
    
    /*voigt1a->SetParLimits (1,-4.25,-4.1);
    voigt1a->SetParLimits (2,0.01,0.3);//width of the gaussian
    voigt1a->SetParLimits (4,0.05,0.5);//width of the lorentzian
    voigt1a->SetParLimits (5,0.01,1);//voigt fraction
    
    voigt1b->SetParLimits (1,-3.15,-2.9);
    voigt1b->SetParLimits (2,0.01,0.4);//width of the gaussian
    voigt1b->SetParLimits (4,0.05,0.5);//width of the lorentzian
    voigt1b->SetParLimits (5,0.01,1);//voigt fraction
    
    voigt1c->SetParLimits (1,-2.8,-2.5);
    voigt1c->SetParLimits (2,0.01,0.4);//width of the gaussian
    voigt1c->SetParLimits (4,0.01,0.5);//width of the lorentzian
    voigt1c->SetParLimits (5,0.01,1);//voigt fraction
    */
    
    /*voigt1c->SetParLimits (1,-2.7,-2.5);
    voigt1c->SetParLimits (2,0.01,0.4);//width of the gaussian
    voigt1c->SetParLimits (4,0.01,0.5);//width of the lorentzian
    voigt1c->SetParLimits (5,0.01,1);//voigt fraction
    */
    
    
    /*voigt1d->SetParLimits (1,-2.4,-2.2);
    voigt1d->SetParLimits (2,0.01,0.3);//width of the gaussian
    voigt1d->SetParLimits (4,0.05,0.5);//width of the lorentzian
    voigt1d->SetParLimits (5,0.01,1);//voigt fraction
    */
    
    /*hBeQ->Fit ( voigt1a, "R");
	hBeQ->Fit ( voigt1b, "R+");
	hBeQ->Fit ( voigt1c, "R+");*/
    //hBeQ->Fit ( voigt1d, "R+");
    
    /*Double_t p[20];
    voigt1a->GetParameters(&p[0]);
    voigt1b->GetParameters(&p[6]);
    voigt1c->GetParameters(&p[12]);
    //voigt1d->GetParameters(&p[18]);
    
    voigt3->SetParameters(p);*/
    

    
	//hBeQ->Fit ( "total","R+");
    hBeQ->Fit ( "voigt3","+","",-4.6,-2.35);
    hBeQ->GetFunction("voigt3")->SetLineColor(kBlack);
    hBeQ->GetFunction("voigt3")->SetLineWidth(3);
    
    
    TF1* voigt1a = new TF1("voigt1a",voig1,-4.35,-3.95,6);
    voigt1a->SetParameters(voigt3->GetParameter(0),voigt3->GetParameter(1),voigt3->GetParameter(2),voigt3->GetParameter(3),voigt3->GetParameter(4),voigt3->GetParameter(13));
    TF1* voigt1b = new TF1("voigt1b",voig1,-3.2,-2.8,6);
    voigt1b->SetParameters(voigt3->GetParameter(5),voigt3->GetParameter(6),voigt3->GetParameter(2),voigt3->GetParameter(7),voigt3->GetParameter(8),voigt3->GetParameter(13));
    TF1* voigt1c = new TF1("voigt1c",voig1,-2.74,-2.4,6);
    voigt1c->SetParameters(voigt3->GetParameter(9),voigt3->GetParameter(10),voigt3->GetParameter(2),voigt3->GetParameter(11),voigt3->GetParameter(12),voigt3->GetParameter(13));
    
    
    voigt1a->Draw("same");
    voigt1a->SetLineWidth(3);
    voigt1b->Draw("same");
    voigt1b->SetLineWidth(3);
    voigt1c->Draw("same");
    voigt1c->SetLineWidth(3);
    
    
	hBeQ->GetXaxis()->SetRangeUser(-5,-2);
    hBeQ->SetTitle("");
    //hBeQ->SetTitle("{}^{12}Be(d,p)^{13}Be Q Value Spectrum");
    hBeQ->GetXaxis()->SetTitle("Q Value [MeV]");
    hBeQ->GetYaxis()->SetTitle("Counts");
    hBeQ->GetXaxis()->SetLabelSize(0.05);
    hBeQ->GetXaxis()->SetTitleSize(0.05);
    hBeQ->GetXaxis()->SetTitleOffset(0.9);
    hBeQ->GetYaxis()->SetLabelSize(0.05);
    hBeQ->GetYaxis()->SetTitleSize(0.05);
    hBeQ->GetYaxis()->SetTitleOffset(0.9);
}
