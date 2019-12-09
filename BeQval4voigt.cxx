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


double exp_background(double *x, double *p) {
    return p[0]*TMath::Exp(p[1]*x[0]);
}

double poly5_background(double *x, double *p) {
    return p[0]*pow(x[0],5);
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

double voig1_poly5Background(double *x, double *p) {
    return p[5]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[5])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[6]*pow(x[0],5);
}

double voig1_poly6Background(double *x, double *p) {
    return p[5]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[5])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[6]*pow(x[0],6);
}


double voig4(double *x, double *p) {
    return p[17]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[17])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[18]*(p[5]*exp(-pow((x[0]-p[6]),2)/(2*p[2]*p[2])))+(1-p[18])*(p[7]/(pow((x[0]-p[6]),2)+pow((p[8]/2),2)))+p[19]*(p[9]*exp(-pow((x[0]-p[10]),2)/(2*p[2]*p[2])))+(1-p[19])*(p[11]/(pow((x[0]-p[10]),2)+pow((p[12]/2),2)))+p[20]*(p[13]*exp(-pow((x[0]-p[14]),2)/(2*p[21]*p[21])))+(1-p[20])*(p[15]/(pow((x[0]-p[14]),2)+pow((p[16]/2),2)))+p[22]*pow(x[0],5);
}

void BeQval4voigt(){
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
    

    TFile *g1 = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod.root","READ");
	TFile *g = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n2.root","READ");
    TFile *fg = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2.root","READ");
    TFile *g4 = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1.root","READ");
    TFile *gh1 = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2_FirstHalf.root","READ");
    TFile *gh2 = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2_SecondHalf.root","READ");


	TF1* fit_func2 = new TF1 ( "fit_func2","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp([6]*x)",-10,10);
    TF1* fit_func3 = new TF1 ( "fit_func3","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp([8]*x)",-10,10);
	TF1* fit_func4 = new TF1 ( "fit_func4","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp([10]*x)",-10,10);
    TF1* lor2 = new TF1 ("lor2","([0]*[2]*[2])/(pow((x-[1]),2)+pow(([2]/2),2))+([3]*[2]*[2])/(pow((x-[4]),2)+pow(([5]/2),2))+[6]*TMath::Exp([7]*x)",-10,10);
    TF1* lor3 = new TF1 ("lor3","([0]*[2]*[2])/(pow((x-[1]),2)+pow(([2]/2),2))+([3]*[2]*[2])/(pow((x-[4]),2)+pow(([2]/2),2))+([5]*[2]*[2])/(pow((x-[6]),2)+pow(([2]/2),2))+[7]*TMath::Exp([8]*x)",-10,10);
    
    TF1* voigt1 = new TF1("voigt1",voig1,-10,10,6);
    
    TF1* voigt4 = new TF1("voigt4",voig4,-10,10,23);
    TF1* voigt4I = new TF1("voigt4I",voig4,-10,10,23);

    
    TH1D *hBeQ = (TH1D*)g1->Get("hQval"); //Interchanged fg and g1. commentend out rebin(8)
    TH1D *hBeQ1 = (TH1D*)fg->Get("hQval");
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
	
	// 4-Voigt Fit for Beryllium
    voigt4->SetParLimits (0,0.01,15);//Amplitude of first gaussian
    voigt4->SetParLimits (1,-4.3,-4.13); //( 1,-4.18,-4.05);//centroid of first gaussian and lorentzian
	voigt4->SetParLimits (2,0.01,0.26);//SD of all the gaussians except for the fourth
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
    voigt4->SetParLimits (16,0.01,0.8);//width of fourth lorentzian
    voigt4->SetParLimits (17,0.001,1);//voigt1 fraction
    voigt4->SetParLimits (18,0.001,0.9);//voigt2 fraction
    voigt4->SetParLimits (19,0.001,0.9);//voigt3 fraction
    voigt4->SetParLimits (20,0.001,0.9);//voigt4 fraction
	voigt4->SetParLimits (21,0.001,0.2);//SD of the fourth gaussian
    voigt4->SetParLimits (22,0,1);//Background amplitude
    
    
    ////////////////////////BACKGROUND///////////////////////////////
    
    double simQval[100], dataQval[100], simBackground[100], NormSimBackground[100], dataCounts[100], newCounts[99];
    TH1D *gh = new TH1D("NewQval","",100,-7,0);
    TH1D *gd = new TH1D("DataQval","",100,-7,0);
    TH1D *gb = new TH1D("BackgroundQval","",100,-7,0);
    
    ifstream backgroundQval;
    backgroundQval.open("/home/jerome/12Be_exp/scripts/backgroundQValue1.out");
    
    if(backgroundQval.is_open()){
        for(int i=0;i<100;i++){
            backgroundQval >> simQval[i] >> simBackground[i];
            NormSimBackground[i] = simBackground[i]/(18242/35);//Background normalized to a maximum value of 20 by comparing to the data
            dataQval[i] = hBeQ->GetXaxis()->GetBinCenter(i);
            dataCounts[i] = hBeQ->GetBinContent(i);
            //cout << dataCounts[i] << endl;
            //cout << NormSimBackground[i] << endl;
        }
    }else{
        cout << "No file found " << endl;
    }
    backgroundQval.close();
    
    for (int i=0;i<99;i++){
        newCounts[i]=dataCounts[i+1]-NormSimBackground[i]; //i+1 and i because there is a difference in the indexes between data and simulation.
        
        for (int j=0;j<newCounts[i];j++){
            gh->Fill(simQval[i]);
        }
        
        for (int j=0;j<dataCounts[i+1];j++){
            gd->Fill(simQval[i]);
        }
        
        for (int j=0;j<NormSimBackground[i];j++){
            gb->Fill(simQval[i]);
        }
        //cout << newCounts[i] << endl;
    }    
    
    TCanvas *c1 = new TCanvas ( "c1" );
    gh->Draw();
    gb->Draw("same");
    gd->Draw("same");
    gh->SetLineWidth(3);
    gb->SetLineWidth(2);
    gd->SetLineWidth(2);
    gb->SetLineColor(kBlack);
    gd->SetLineColor(kRed);
    gh->GetYaxis()->SetRangeUser(0,50);
    gh->GetXaxis()->SetRangeUser(-5,-2);
    gh->GetXaxis()->SetTitle("Q Value [MeV]"); gh->GetXaxis()->CenterTitle();
    gh->GetYaxis()->SetTitle("Counts"); gh->GetYaxis()->CenterTitle();
    gh->SetStats(0);
    
    auto legend1 = new TLegend(0.7,0.7,0.7,0.7);
	legend1->AddEntry(gh,"Data with background subtracted","l");
	legend1->AddEntry(gb,"Background","l");
    legend1->AddEntry(gd,"Data","l");
    legend1->Draw();
    
    ////////////////////////BACKGROUND///////////////////////////////
        
    TCanvas *c2 = new TCanvas ( "c2" ); 
    /*
    voigt4->FixParameter (0,0.01);//Amplitude of first gaussian
    voigt4->FixParameter (1,-4.13); //( 1,-4.18,-4.05);//centroid of first gaussian and lorentzian
	voigt4->FixParameter (2,0.0653099);//SD of all the gaussians except for the fourth
	voigt4->FixParameter (3,0.337106);//Amplitude of first lorentzian
    voigt4->FixParameter (4,0.3);//width of first lorentzian
    voigt4->FixParameter (5,5.0674);//Amplitude of second gaussian
	voigt4->FixParameter (6,-2.974);//centroid of second gaussian and lorentzian
    voigt4->FixParameter (7,3.26862e-10);//Amplitude of second lorentzian
    voigt4->FixParameter (8,0.146949);//width of second lorentzian
    voigt4->FixParameter (9,16.7364);//Amplitude of third gaussian
    voigt4->FixParameter (10,-2.65129);//centroid of third gaussian and lorentzian
    voigt4->FixParameter (11,6.86478e-09);//Amplitude of third lorentzian
    voigt4->FixParameter (12,0.990387);//width of third lorentzian
	voigt4->FixParameter (13,5.92357);//Amplitude of fourth gaussian
    voigt4->FixParameter (14,-2.4);//centroid of fourth gaussian and lorentzian
    voigt4->FixParameter (15,1.27345e-07);//Amplitude of fourth lorentzian
    voigt4->FixParameter (16,0.397545);//width of fourth lorentzian
    voigt4->FixParameter (17,0.0100369);//voigt1 fraction
    voigt4->FixParameter (18,0.885955);//voigt2 fraction
    voigt4->FixParameter (19,0.572596);//voigt3 fraction
    voigt4->FixParameter (20,0.731534);//voigt4 fraction
	voigt4->FixParameter (21,0.153469);//SD of the fourth gaussian
    //voigt4->SetParLimits (22,0.00267,0.005);//Background amplitude
    */
    
    hBeQ->Draw();
	hBeQ->Rebin(8);
    hBeQ->SetLineWidth(2);
    
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
    ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(1000000);
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.001);
    
    hBeQ->Fit ( "voigt4","+","",-4.35,-2.15);
    hBeQ->GetFunction("voigt4")->SetLineColor(kBlack);
    hBeQ->GetFunction("voigt4")->SetLineWidth(3);
    
    
    Double_t p[23];
    voigt4->GetParameters(&p[0]);
    
    TF1* voigt1a = new TF1("voigt1a",voig1,-4.35,-3.9,6);
    voigt1a->SetParameters(p[0],p[1],p[2],p[3],p[4],p[17]);
    TF1* voigt1b = new TF1("voigt1b",voig1,-3.2,-2.8,6);
    voigt1b->SetParameters(p[5],p[6],p[2],p[7],p[8],p[18]);
    TF1* voigt1c = new TF1("voigt1c",voig1,-2.9,-2.55,6);
    voigt1c->SetParameters(p[9],p[10],p[2],p[11],p[12],p[19]);
    TF1* voigt1d = new TF1("voigt1d",voig1,-2.5,-2.15,6);
    voigt1d->SetParameters(p[13],p[14],p[21],p[15],p[16],p[20]);
    TF1* background = new TF1("background",poly5_background,-4.35,-2.15,2);
    background->SetParameter(0,p[22]);
    
    TF1* totala = new TF1("totala",voig1_poly6Background,-4.35,-3.9,7);
    totala->SetParameters(p[0],p[1],p[2],p[3],p[4],p[17],p[22]);
    TF1* totalb = new TF1("totalb",voig1_poly6Background,-3.2,-2.8,7);
    totalb->SetParameters(p[5],p[6],p[2],p[7],p[8],p[18],p[22]);
    TF1* totalc = new TF1("totalc",voig1_poly6Background,-2.9,-2.55,7);
    totalc->SetParameters(p[9],p[10],p[2],p[11],p[12],p[19],p[22]);
    TF1* totald = new TF1("totald",voig1_poly6Background,-2.5,-2.15,7);
    totald->SetParameters(p[13],p[14],p[21],p[15],p[16],p[20],p[22]);
    
    //voigt1a->Draw("same");
    voigt1a->SetLineWidth(3);
    //voigt1b->Draw("same");
    voigt1b->SetLineWidth(3);
    //voigt1c->Draw("same");
    voigt1c->SetLineWidth(3);
    //voigt1d->Draw("same");
    voigt1d->SetLineWidth(3);
    /*
    totala->Draw("same");
    totala->SetLineWidth(3);
    totala->SetLineColor(kBlue);
    
    totalb->Draw("same");
    totalb->SetLineWidth(3);
    totalb->SetLineColor(kBlue);
    
    totalc->Draw("same");
    totalc->SetLineWidth(3);
    totalc->SetLineColor(kBlue);
    
    totald->Draw("same");
    totald->SetLineWidth(3);
    totald->SetLineColor(kBlue);
    */
    
    //background->Draw("same");
    //background->SetLineWidth(3);
    //background->SetLineColor(3);
    
    
    hBeQ->GetXaxis()->SetRangeUser(-5,-2);
    hBeQ->SetTitle("");
    hBeQ->GetXaxis()->SetTitle("Q Value [MeV]");
    hBeQ->GetYaxis()->SetTitle("Counts");
    //hBeQ->GetXaxis()->SetLabelSize(0.05);
    hBeQ->GetXaxis()->SetTitleSize(0.05);
    hBeQ->GetXaxis()->SetTitleOffset(0.9);
    //hBeQ->GetYaxis()->SetLabelSize(0.05);
    hBeQ->GetYaxis()->SetTitleSize(0.05);
    hBeQ->GetYaxis()->SetTitleOffset(0.9);
    hBeQ->SetStats(0);
    
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
    
    // 4-Voigt Fit for Beryllium for separate ring groups
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
    
    /*TGraphErrors *gr1[4];
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
    */
}
