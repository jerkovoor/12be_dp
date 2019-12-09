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
#include<bits/stdc++.h> 
#include <Fit/Fitter.h>


string getString(char x) 
{ 
    string s(5, x); 
  
    return s;    
} 

double poly5_background(double *x, double *p) {
    return p[0]*pow(x[0],5);
}

double voig1(double *x, double *p) {
    return p[5]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[5])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)));
}


double voig1_poly5Background(double *x, double *p) {
    return p[5]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[5])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[6]*pow(x[0],5);
}

double voig4(double *x, double *p) { // Different standard deviation for the ground state
    return p[17]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[17])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[18]*(p[5]*exp(-pow((x[0]-p[6]),2)/(2*p[2]*p[2])))+(1-p[18])*(p[7]/(pow((x[0]-p[6]),2)+pow((p[8]/2),2)))+p[19]*(p[9]*exp(-pow((x[0]-p[10]),2)/(2*p[2]*p[2])))+(1-p[19])*(p[11]/(pow((x[0]-p[10]),2)+pow((p[12]/2),2)))+p[20]*(p[13]*exp(-pow((x[0]-p[14]),2)/(2*p[2]*p[2])))+(1-p[20])*(p[15]/(pow((x[0]-p[14]),2)+pow((p[16]/2),2)));
}

/*
double voig4_poly5Background(double *x, double *p) { // Different standard deviation for the ground state
    return p[17]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[17])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[18]*(p[5]*exp(-pow((x[0]-p[6]),2)/(2*p[2]*p[2])))+(1-p[18])*(p[7]/(pow((x[0]-p[6]),2)+pow((p[8]/2),2)))+p[19]*(p[9]*exp(-pow((x[0]-p[10]),2)/(2*p[2]*p[2])))+(1-p[19])*(p[11]/(pow((x[0]-p[10]),2)+pow((p[12]/2),2)))+p[20]*(p[13]*exp(-pow((x[0]-p[14]),2)/(2*p[21]*p[21])))+(1-p[20])*(p[15]/(pow((x[0]-p[14]),2)+pow((p[16]/2),2)))+p[22]*pow(x[0],5);
}
*/

double voig4_poly5Background(double *x, double *p) { // Different standard deviation for the ground state
    return p[17]*(p[0]*exp(-pow((x[0]-p[1]),2)/(2*p[2]*p[2])))+(1-p[17])*(p[3]/(pow((x[0]-p[1]),2)+pow((p[4]/2),2)))+p[18]*(p[5]*exp(-pow((x[0]-p[6]),2)/(2*p[2]*p[2])))+(1-p[18])*(p[7]/(pow((x[0]-p[6]),2)+pow((p[8]/2),2)))+p[19]*(p[9]*exp(-pow((x[0]-p[10]),2)/(2*p[2]*p[2])))+(1-p[19])*(p[11]/(pow((x[0]-p[10]),2)+pow((p[12]/2),2)))+p[20]*(p[13]*exp(-pow((x[0]-p[14]),2)/(2*p[2]*p[2])))+(1-p[20])*(p[15]/(pow((x[0]-p[14]),2)+pow((p[16]/2),2)))+p[21]*pow(x[0],5);
}

void BeQvalue_Test(){
    
    
    TFile *f1 = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/beryllium/Be_pedestal_TRIUMF_DL_BeamOffset_0_0_Shift_0_IC_Cut_TargetDistance80.00mm_TargetThickness43.00um_Yu_Mod_UnevenBinning_Sd1rAlphaCal_Sd2rNewInBeamCal.root","READ");
    TFile *f2 = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/beryllium/Be_pedestal_TRIUMF_DL_BeamOffset_0_0_Shift_0_IC_Cut_TargetDistance80.00mm_TargetThickness49.99um_Yu_Mod_UnevenBinning_Sd1rAlphaCal_Sd2rNewInBeamCal.root","READ");
    TFile *f3 = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/beryllium/Be_pedestal_TRIUMF_DL_BeamOffset_0_0_Shift_0.0000_CutIC_TargetDistance80.00mm_TargetThickness43.00um_Yu_Mod_UnevenBinning_Sd1rAlphaCal_Sd2rNewInBeamCal_NewSectorGeometry.root","READ");
    TFile *f4 = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/beryllium/Be_pedestal_TRIUMF_DL_BeamOffset_-3_4_Shift_0_IC_Cut_TargetDistance80.00mm_TargetThickness49.99um_Yu_Mod_UnevenBinning_Sd1rAlphaCal_Sd2rNewInBeamCal.root","READ");
    
    TFile *b3 = new TFile("/home/jerome/12Be_exp/Analysis/QvalBackground/calculated_qvalue_43um_offset.root","READ");
    TCanvas *c1 = (TCanvas*)b3->Get("c1");
    TH1D *hbg3 = (TH1D*)c1->GetPrimitive("calculatedQValue");
    
    
    TH1D *hBeQ1 = (TH1D*)f1->Get("hQval");
    TH1D *hBeQ2 = (TH1D*)f2->Get("hQval");
    TH1D *hBeQ3 = (TH1D*)f3->Get("hQval");
    TH1D *hBeQ4 = (TH1D*)f4->Get("hQval");

    
    TF1* voigt1_poly5Background = new TF1("voigt1",voig1_poly5Background,-10,10,6);
    TF1* voigt4_poly5Background = new TF1("voigt4_poly5Background",voig4_poly5Background,-10,10,23);
    TF1* voigt4 = new TF1("voigt4",voig4,-10,10,21);
    
    
    // 4-Voigt with 5 degree polynomial Fit for Beryllium
    
    //voigt4_poly5Background->SetParLimits (0,0.01,15);//Amplitude of first gaussian
    voigt4_poly5Background->SetParLimits (1,-4.2,-3.9); //( 1,-4.18,-4.05);//centroid of first gaussian and lorentzian
	//voigt4_poly5Background->FixParameter (2,0.240131);//SD of all the gaussians except for the fourth
    voigt4_poly5Background->SetParLimits (2,0.1,0.3);
    
	voigt4_poly5Background->SetParLimits (3,0.01,1);//Amplitude of first lorentzian
    voigt4_poly5Background->SetParLimits (4,0.01,0.5);//width of first lorentzian
    
    voigt4_poly5Background->SetParLimits (5,0.01,10);//Amplitude of second gaussian
	voigt4_poly5Background->SetParLimits (6,-3.3,-3.0);//centroid of second gaussian and lorentzian
    
    voigt4_poly5Background->SetParLimits (7,0,0.05);//Amplitude of second lorentzian
    voigt4_poly5Background->SetParLimits (8,0.1,10);//width of second lorentzian
    
    //voigt4_poly5Background->SetParLimits (9,0,10);//Amplitude of third gaussian
    voigt4_poly5Background->SetParLimits (10,-2.8,-2.6);//centroid of third gaussian and lorentzian
    
    //voigt4_poly5Background->SetParLimits (11,0,0.001);//Amplitude of third lorentzian
    voigt4_poly5Background->SetParLimits (12,0.01,5);//width of third lorentzian
    
	voigt4_poly5Background->SetParLimits (13,0,0.05);//Amplitude of fourth gaussian
    voigt4_poly5Background->SetParLimits (14,-2.7,-2.3);//centroid of fourth gaussian and lorentzian
    //voigt4_poly5Background->SetParLimits (21,0.001,0.3);//SD of the fourth gaussian
    
    voigt4_poly5Background->SetParLimits (15,0,0.5);//Amplitude of fourth lorentzian
    voigt4_poly5Background->SetParLimits (16,0.01,10);//width of fourth lorentzian
    
    voigt4_poly5Background->SetParLimits (17,0.001,1);//voigt1 fraction
    voigt4_poly5Background->SetParLimits (18,0.001,1);//voigt2 fraction
    voigt4_poly5Background->SetParLimits (19,0.001,1);//voigt3 fraction
    voigt4_poly5Background->SetParLimits (20,0.001,1);//voigt4 fraction
	
    //voigt4_poly5Background->SetParLimits (22,0,1);//Background amplitude
    
    
    TF1* voigt1a = new TF1("voigt1a",voig1,-4.35,-3.5,6);
    TF1* voigt1b = new TF1("voigt1b",voig1,-3.7,-3.3,6);
    TF1* voigt1c = new TF1("voigt1c",voig1,-3.3,-2.7,6);
    TF1* voigt1d = new TF1("voigt1d",voig1,-2.8,-2.15,6);
    
    TF1* voigt2a = new TF1("voigt2a",voig1,-4.35,-3.5,6);
    TF1* voigt2b = new TF1("voigt2b",voig1,-3.7,-3.3,6);
    TF1* voigt2c = new TF1("voigt2c",voig1,-3.3,-2.7,6);
    TF1* voigt2d = new TF1("voigt2d",voig1,-2.8,-2.15,6);
    
    TF1* voigt3a = new TF1("voigt3a",voig1,-4.3,-3.8,6);
    TF1* voigt3b = new TF1("voigt3b",voig1,-3.2,-2.9,6);
    TF1* voigt3c = new TF1("voigt3c",voig1,-2.8,-2.55,6);
    TF1* voigt3d = new TF1("voigt3d",voig1,-2.6,-2.25,6);
    
    TF1* voigt4a = new TF1("voigt4a",voig1,-4.35,-3.5,6);
    TF1* voigt4b = new TF1("voigt4b",voig1,-3.7,-3.3,6);
    TF1* voigt4c = new TF1("voigt4c",voig1,-3.3,-2.7,6);
    TF1* voigt4d = new TF1("voigt4d",voig1,-2.8,-2.15,6);
    TF1* background = new TF1("background",poly5_background,-4.35,-2.15,2);
    
    
    TCanvas *d1 = new TCanvas ( "d1" );
    
    d1->Divide(2,2);
    /*
    for(int i=0;i<4;i++){
        d1->cd(i+1);
        char *k = Form("hBeQ%i",i+1);
        cout << getString(&k) << endl;//->Draw();
    }
    */
    
    /*
    d1->cd(1);
    hBeQ1->Draw();
    hBeQ1->Fit ( "voigt4_poly5Background","+","",-4.1,-2.15);
    hBeQ1->GetFunction("voigt4_poly5Background")->SetLineColor(kBlack);
    hBeQ1->GetFunction("voigt4_poly5Background")->SetLineWidth(3);
    hBeQ1->SetTitle("Target distance 80mm, Thickness 43um, No offset");
    hBeQ1->GetYaxis()->SetTitle("Counts/75 keV");
    hBeQ1->SetStats(0);
    
    Double_t p1[23];
    voigt4_poly5Background->GetParameters(&p1[0]);
    
    
    voigt1a->SetParameters(p1[0],p1[1],p1[2],p1[3],p1[4],p1[17]);
    
    voigt1b->SetParameters(p1[5],p1[6],p1[2],p1[7],p1[8],p1[18]);
    
    voigt1c->SetParameters(p1[9],p1[10],p1[2],p1[11],p1[12],p1[19]);
    
    voigt1d->SetParameters(p1[13],p1[14],p1[21],p1[15],p1[16],p1[20]);
    
    background->SetParameter(0,p1[22]);
    
    
    voigt1a->Draw("same");
    voigt1a->SetLineWidth(3);
    voigt1b->Draw("same");
    voigt1b->SetLineWidth(3);
    voigt1c->Draw("same");
    voigt1c->SetLineWidth(3);
    voigt1d->Draw("same");
    voigt1d->SetLineWidth(3);
    //background->Draw("same");
    
    
    
    
    
    d1->cd(2);
    hBeQ2->Draw();
    hBeQ2->Fit ( "voigt4_poly5Background","+","",-4.35,-2.15);
    hBeQ2->GetFunction("voigt4_poly5Background")->SetLineColor(kBlack);
    hBeQ2->GetFunction("voigt4_poly5Background")->SetLineWidth(3);
    hBeQ2->SetTitle("Target distance 80mm, Thickness 49.99um, No offset");
    hBeQ2->GetYaxis()->SetTitle("Counts/75 keV");
    hBeQ2->SetStats(0);
    
    Double_t p2[23];
    voigt4_poly5Background->GetParameters(&p2[0]);
    

    voigt2a->SetParameters(p2[0],p2[1],p2[2],p2[3],p2[4],p2[17]);

    voigt2b->SetParameters(p2[5],p2[6],p2[2],p2[7],p2[8],p2[18]);

    voigt2c->SetParameters(p2[9],p2[10],p2[2],p2[11],p2[12],p2[19]);

    voigt2d->SetParameters(p2[13],p2[14],p2[21],p2[15],p2[16],p2[20]);
    
    background->SetParameter(0,p2[22]);

    
    voigt2a->Draw("same");
    voigt2a->SetLineWidth(3);
    voigt2b->Draw("same");
    voigt2b->SetLineWidth(3);
    voigt2c->Draw("same");
    voigt2c->SetLineWidth(3);
    voigt2d->Draw("same");
    voigt2d->SetLineWidth(3);
    //background->Draw("same");
    
    
    
    
    
    d1->cd(3);
    
    hBeQ3->Draw();
    hBeQ3->Fit ( "voigt4_poly5Background","+","",-4.3,-2.15);
    hBeQ3->GetFunction("voigt4_poly5Background")->SetLineColor(kBlack);
    hBeQ3->GetFunction("voigt4_poly5Background")->SetLineWidth(3);
    hBeQ3->SetTitle("Target distance 80mm, Thickness 43um, With offset");
    hBeQ3->GetYaxis()->SetTitle("Counts/50 keV");
    hBeQ3->SetStats(0);
    
    Double_t p3[23];
    voigt4_poly5Background->GetParameters(&p3[0]);
    
    Double_t fitError3[23];
    //voigt4_poly5Background->GetParErrors(&fitError3[0]);
    

    voigt3a->SetParameters(p3[0],p3[1],p3[2],p3[3],p3[4],p3[17]);

    voigt3b->SetParameters(p3[5],p3[6],p3[2],p3[7],p3[8],p3[18]);

    voigt3c->SetParameters(p3[9],p3[10],p3[2],p3[11],p3[12],p3[19]);

    //voigt3d->SetParameters(p3[13],p3[14],p3[21],p3[15],p3[16],p3[20]);
    
    voigt3d->SetParameters(p3[13],p3[14],p3[2],p3[15],p3[16],p3[20]);
    
    background->SetParameter(0,p3[22]);

    
    voigt3a->Draw("same");
    voigt3a->SetLineWidth(3);
    voigt3b->Draw("same");
    voigt3b->SetLineWidth(3);
    voigt3c->Draw("same");
    voigt3c->SetLineWidth(3);
    voigt3d->Draw("same");
    voigt3d->SetLineWidth(3);
    //background->Draw("same");
    
    
    
    
    
    d1->cd(4);
    hBeQ4->Draw();
    hBeQ4->Fit ( "voigt4_poly5Background","+","",-4.35,-2.15);
    hBeQ4->GetFunction("voigt4_poly5Background")->SetLineColor(kBlack);
    hBeQ4->GetFunction("voigt4_poly5Background")->SetLineWidth(3);
    hBeQ4->SetTitle("Target distance 80mm, Thickness 49.99um, With offset");
    hBeQ4->GetYaxis()->SetTitle("Counts/75 keV");
    hBeQ4->SetStats(0);
    
    Double_t p4[23];
    voigt4_poly5Background->GetParameters(&p4[0]);
    

    voigt4a->SetParameters(p4[0],p4[1],p4[2],p4[3],p4[4],p4[17]);

    voigt4b->SetParameters(p4[5],p4[6],p4[2],p4[7],p4[8],p4[18]);

    voigt4c->SetParameters(p4[9],p4[10],p4[2],p4[11],p4[12],p4[19]);

    voigt4d->SetParameters(p4[13],p4[14],p4[21],p4[15],p4[16],p4[20]);
    
    background->SetParameter(0,p4[22]);

    
    voigt4a->Draw("same");
    voigt4a->SetLineWidth(3);
    voigt4b->Draw("same");
    voigt4b->SetLineWidth(3);
    voigt4c->Draw("same");
    voigt4c->SetLineWidth(3);
    voigt4d->Draw("same");
    voigt4d->SetLineWidth(3);
    //background->Draw("same");
    
    //cout << p1[1] << "  " << p1[6] << "  " << p1[10] << "  " << p1[14] << endl;
    //cout << p2[1] << "  " << p2[6] << "  " << p2[10] << "  " << p2[14] << endl;
    cout << p3[1] << "  " << p3[6] << "  " << p3[10] << "  " << p3[14] << endl;
    //cout << p4[1] << "  " << p4[6] << "  " << p4[10] << "  " << p4[14] << endl;
    
    
    */
    
    
    hbg3->Rebin(2);
    //hbg3->Draw();
    
    
    
    
    

    Int_t hbg3_size = hbg3->GetSize();
    
    double En[hbg3_size], Counts[hbg3_size];
    
    float scale_factor = 105;
    
    
    TH1D *hB = new TH1D("hB",Form("Scaled Background (Scale factor = %.0f)",scale_factor),50,-5,-2);
    hB->Sumw2();
    
    //cout << hbg3_size << endl;
    
    for (int i=0;i<hbg3_size;i++){
        En[i] = hbg3->GetBinCenter(i);
        Counts[i] = hbg3->GetBinContent(i);
        for (int j=0;j<floor(Counts[i]/scale_factor);j++){
            hB->Fill(En[i]);
        }
        //cout << En[i] << "      " << Counts[i] << endl;
    }
    
    
    //hBeQ3->SetLineColor(kRed);
    
    
    /////////////////////Scaled background///////////////////////////
    
    Int_t hB_size = hB->GetSize();
    
    double En_hB[hB_size], Counts_hB[hB_size];
    
    for (int i=0;i<hB_size;i++){
        En_hB[i] = hB->GetBinCenter(i);
        Counts_hB[i] = hB->GetBinContent(i);
        //cout << En_hB[i] << "      " << Counts_hB[i] << endl;
    }
    
    ////////////////////Original Q value spectrum/////////////////////
    
    
    Int_t hBeQ3_size = hBeQ3->GetSize();
    
    double En_hBeQ3[hBeQ3_size], Counts_hBeQ3[hBeQ3_size];
    
    TGraphErrors *grQ = new TGraphErrors();
    
    for (int i=0;i<hBeQ3_size;i++){
        En_hBeQ3[i] = hBeQ3->GetBinCenter(i);
        Counts_hBeQ3[i] = hBeQ3->GetBinContent(i);
        
        grQ->SetPoint(i,En_hBeQ3[i],Counts_hBeQ3[i]);
        grQ->SetPointError(i,0,hBeQ3->GetBinError(i));
        
        //cout << En_hBeQ3[i] << "      " << Counts_hBeQ3[i] << "      " << hBeQ3->GetBinError(i) << endl;
    }
    
    
    TCanvas *c2 = new TCanvas ( "c2" );
    //hBeQ3->Draw();
    grQ->Draw("*AP");
    grQ->SetTitle(Form("13Be Q value spectrum (Blue) with a scaled background (Red)(Scale factor = %.0f)",scale_factor));
    grQ->GetYaxis()->SetRangeUser(0.,45.);
    grQ->GetXaxis()->SetRangeUser(-4.5,-2.);
    hB->SetLineWidth(3);
    hBeQ3->SetLineWidth(3);
    hB->Draw("same");
    hBeQ3->SetTitle(Form("13Be Q value spectrum (Blue) with a scaled background (Red)(Scale factor = %.0f)",scale_factor));
    hB->SetLineColor(kRed);
    hB->GetYaxis()->SetRangeUser(0.,45.);


    
    
    ////////////////////Background subtracted Q value spectrum//////////////////
    
    TH1D *hBSQValue = new TH1D("hB",Form("Background Subtracted Q Value Spectrum (Scale factor = %.0f)",scale_factor),50,-5,-2);
    
    TGraphErrors *grQBS = new TGraphErrors();
    
    double En_hBeQ3BS[hBeQ3_size], Counts_hBeQ3BS[hBeQ3_size];
    
    for (int i=0;i<hBeQ3_size;i++){
        
        Counts_hBeQ3BS[i] = fabs(Counts_hBeQ3[i]-Counts_hB[i]);
        
        grQBS->SetPoint(i,En[i],Counts_hBeQ3BS[i]);
        grQBS->SetPointError(i,0,sqrt(Counts_hBeQ3[i]+Counts_hB[i]));//sum of the quadratures of the errors
        
        for (int j=0;j<Counts_hBeQ3BS[i];j++){
            hBSQValue->Fill(En[i]);
        }
        //cout << En[i] << "      " << Counts_hBeQ3BS[i] << endl;
    }
    
    
    // 4-Voigt Fit for Beryllium
    
    //voigt4->SetParLimits (0,1,20);//Amplitude of first gaussian
    //voigt4->FixParameter (0,60);
    voigt4->SetParLimits (1,-4.2,-3.9); //( 1,-4.18,-4.05);//centroid of first gaussian and lorentzian
	//voigt4->FixParameter (2,0.240131);//SD of all the gaussians except for the fourth
    voigt4->SetParLimits (2,0.2,0.3);
    
	voigt4->SetParLimits (3,0,10);//Amplitude of first lorentzian
    //voigt4->FixParameter (3,50);
    voigt4->SetParLimits (4,0.01,10);//width of first lorentzian
    
    //voigt4->SetParLimits (5,1,10);//Amplitude of second gaussian
	voigt4->SetParLimits (6,-3.1,-2.9);//centroid of second gaussian and lorentzian
    
    voigt4->SetParLimits (7,0.01,10);//Amplitude of second lorentzian
    voigt4->SetParLimits (8,0.0,1);//width of second lorentzian
    
    voigt4->SetParLimits (9,0.1,20);//Amplitude of third gaussian
    voigt4->SetParLimits (10,-2.7,-2.6);//centroid of third gaussian and lorentzian
    //voigt4->FixParameter (10,-2.65);
    
    voigt4->SetParLimits (11,0.01,10);//Amplitude of third lorentzian
    voigt4->SetParLimits (12,0.0,1);//width of third lorentzian
    
	//voigt4->SetParLimits (13,0,10);//Amplitude of fourth gaussian
    voigt4->SetParLimits (14,-2.5,-2.3);//centroid of fourth gaussian and lorentzian
    //voigt4->FixParameter (14,-2.35);
    //voigt4->SetParLimits (21,0.001,0.3);//SD of the fourth gaussian
    
    //voigt4->SetParLimits (15,0.1,0.5);//Amplitude of fourth lorentzian
    voigt4->SetParLimits (16,0.0,1);//width of fourth lorentzian
    
    voigt4->SetParLimits (17,0,1);//voigt1 fraction
    voigt4->SetParLimits (18,0,1);//voigt2 fraction
    voigt4->SetParLimits (19,0,1);//voigt3 fraction
    voigt4->SetParLimits (20,0,1);//voigt4 fraction
    
    TCanvas *c3 = new TCanvas ( "c3" );
    //hBSQValue->SetLineWidth(3);
    grQBS->Draw("*AP");
    
    grQBS->SetTitle("13Be Q value spectrum after background subtraction");
    
    grQBS->GetYaxis()->SetRangeUser(0.,35.);
    grQBS->Fit ( "voigt4","LR","",-4.4,-2);
    grQBS->GetFunction("voigt4")->SetLineWidth(3);
    grQBS->GetXaxis()->SetRangeUser(-4.5,-2.);
    
    
    
    
    Double_t pBS[21];
    voigt4->GetParameters(&pBS[0]);
    
    
    
    
    ///////////////////////////////////////////////////////
    ///////                                         ///////
    ///////Calculation of energy above the threshold///////
    ///////                                         ///////
    ///////////////////////////////////////////////////////
    
    float amu=931.5;
    float mp=938.28;
    float mn=939.57;
    float me=0.511;
    float m12Be=12.02692;
    float m1H=1.00782;
    float m2H=2.0141;
    
    float Qgs=pBS[14];
    float Q1=pBS[10];
    float Q2=pBS[6];
    float Q3=pBS[1];
    
    double B_12Be=((4*mp)+(8*mn)-((m12Be*amu)-(4*me)));
    double B_2H=((1*mp)+(1*mn)-((m2H*amu)-(1*me)));
    double B_1H=((1*mp)+(0*mn)-((m1H*amu)-(1*me)));
    
    double Q_Zero_Sn=-(B_2H-B_1H);
    double B_13Be=B_12Be+B_2H-B_1H+Qgs;
    double Sn_13Be=B_13Be-B_12Be;
    
    double ET_gs=(Qgs-Qgs)-Sn_13Be;
    double ET_1=(Qgs-Q1)-Sn_13Be;
    double ET_2=(Qgs-Q2)-Sn_13Be;
    double ET_3=(Qgs-Q3)-Sn_13Be;
    
    //cout << Sn_13Be << endl;
    cout << ET_gs << "  " << ET_1 << "  " << ET_2 << "  " << ET_3 << endl;
    cout << pBS[16] << "  " << pBS[12] << "  " << pBS[8] << "  " << pBS[4] << endl;
    cout << pBS[2] << endl;
    
    
    ////////////////////////////////////////////////////////////
}
