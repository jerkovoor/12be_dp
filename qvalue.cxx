using namespace std;

#include <TFile.h> 
#include <TMath.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>

void qvalue() {

	TFile *f = new TFile("../Analysis/Be/AlphaOnly/Be_target_Qvalue.root","READ");
	//TH2D *h_telescope = (TH2D*)f->Get("hCSd1rSd2rICCut");
	TH1D *h_qval_loss = (TH1D*)f->Get("hQval");
	TH1D *h_qval_noloss = (TH1D*)f->Get("hQval1");
	
	TFile *g = new TFile("../Analysis/Be_notarget/Be_notarget_Qvalue.root","READ");
	TH2D *h_telescope_notarget = (TH2D*)g->Get("hCSd1rSd2rICCut");
	TH1D *h_qval_notarget = (TH1D*)g->Get("hQval");
	
	TFile *h = new TFile("../Analysis/Be_notarget/Be_notarget_Qvalue_test.root","READ");
	TH1D *h_qvaltest_notarget = (TH1D*)h->Get("hQval");
	
	TFile *fg = new TFile("../Analysis/Be/AlphaOnly/Be_target_Qvalue_test.root","READ");
	TH1D *h_qvaltest_target = (TH1D*)fg->Get("hQval");
	
	

	int MinBin1=h_qval_loss->GetXaxis()->FindBin(-4.3);
	int MaxBin1=h_qval_loss->GetXaxis()->FindBin(-3.8);
	int MinBin2=h_qval_noloss->GetXaxis()->FindBin(-4.3);
	int MaxBin2=h_qval_noloss->GetXaxis()->FindBin(-3.8);
	int MinBin3=h_qvaltest_notarget->GetXaxis()->FindBin(-4.3);
	int MaxBin3=h_qvaltest_notarget->GetXaxis()->FindBin(-3.8);
	int MinBin5=h_qval_notarget->GetXaxis()->FindBin(-4.3);
	int MaxBin5=h_qval_notarget->GetXaxis()->FindBin(-3.8);
	
	int length1=MaxBin1-MinBin1+1;
	int length2=MaxBin2-MinBin2+1;
	int length3=MaxBin3-MinBin3+1;
	int length5=MaxBin5-MinBin5+1;
	//cout << MinBin1 << "	" << MaxBin1 << "	" << length1 << endl;
	//cout << MinBin2 << "	" << MaxBin2 << "	" << length2 << endl;
	double QValue1[length1];
	double Counts1[length1];
	double error1[length1];
	double QValue2[length2];
	double Counts2[length2];
	double error2[length1];
	double CountsSub[length3];
	double QValue3[length3];
	double Counts3[length3];
	double error3[length1];
	double QValue4[length3];
	double Counts4[length3];
	double error4[length1];
	double QValue5[length5];
	double Counts5[length5];
	double error5[length1];
	
	ofstream qval1;
	qval1.open("/home/jerome/12Be_exp/Analysis/Be/QValue_NoTarget.txt"); //open a .txt file to store the results of the fit; change the path and name acordingly
	qval1 << " Q Value / Counts  \n";
	
	ofstream qval2;
	qval2.open("/home/jerome/12Be_exp/Analysis/Be/QValue_Target.txt"); //open a .txt file to store the results of the fit; change the path and name acordingly
	qval2 << " Q Value / Counts  \n";

	for (int j=MinBin1;j<=MaxBin1;j++){
		QValue1[j-MinBin1]=h_qval_loss->GetXaxis()->GetBinCenter(j);
		Counts1[j-MinBin1]=h_qval_loss->GetBinContent(j);
		error1[j-MinBin1]=h_qval_loss->GetBinError(j);
		//cout << QValue1[j-MinBin1] << "	" << Counts1[j-MinBin1] << "	" << error1[j-MinBin1] << endl;
	/*}
	
	for (int j=MinBin2;j<=MaxBin2;j++){*/
		QValue2[j-MinBin1]=h_qval_notarget->GetXaxis()->GetBinCenter(j);
		Counts2[j-MinBin1]=3.714*h_qval_notarget->GetBinContent(j);
		error2[j-MinBin1]=h_qval_notarget->GetBinError(j);
		CountsSub[j-MinBin1]=Counts1[j-MinBin1]-Counts2[j-MinBin1];
		//cout << QValue2[j-MinBin2] << "	" << Counts2[j-MinBin2] << endl;
	}
	
	for (int j=MinBin3;j<=MaxBin3;j++){
		QValue3[j-MinBin3]=h_qvaltest_notarget->GetXaxis()->GetBinCenter(j);
		Counts3[j-MinBin3]=3.714*h_qvaltest_notarget->GetBinContent(j);
		error3[j-MinBin3]=h_qvaltest_notarget->GetBinError(j);
		qval1 << QValue3[j-MinBin3] << "	" << Counts3[j-MinBin3] << endl;
	
		QValue4[j-MinBin3]=h_qvaltest_target->GetXaxis()->GetBinCenter(j);
		Counts4[j-MinBin3]=h_qvaltest_target->GetBinContent(j);
		error4[j-MinBin3]=h_qvaltest_target->GetBinError(j);
		qval2 << QValue4[j-MinBin3] << "	" << Counts4[j-MinBin3] << endl;
	}
	
	
	
	TF1 *fit_gauss2 = new TF1 ( "fit_gauss2","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))",-10,10 );
	
	TCanvas *c1 = new TCanvas ( "c1" ); //create a canvas
	c1->SetTitle("Different energy loss");
	c1->Divide ( 1,2 );
	c1->cd(1);
	TGraph *gr1 = new TGraphErrors(length3,QValue3,Counts3,0,error3);
	gr1->SetTitle("Q value without target");
	gr1->Draw("AP*");
	//gr1->Fit(fit_gauss2,"E","",QValue1[0],QValue1[length1]);
	c1->cd(2);
	TGraph *gr2 = new TGraphErrors(length3,QValue4,Counts4,0,error4);
	gr2->SetTitle("Q value with target");
	gr2->Draw("AP*");
	
	/*TCanvas *c2 = new TCanvas ( "c2" );
	c2->Divide(2,1);
	c2->cd(1);
	h_telescope_notarget->Draw("colz");
	c2->cd(2);
	h_telescope->Draw("colz");*/
	
	TCanvas *c3 = new TCanvas ( "c3" );
	c3->SetTitle("Different energy loss");
	c3->Divide(2,1);
	c3->cd(1);
	h_qvaltest_notarget->Draw();
	h_qvaltest_notarget->GetXaxis()->SetRangeUser(-4.4,-3.8);
	h_qvaltest_notarget->SetTitle("Q value without target");
	c3->cd(2);
	h_qvaltest_target->Draw();
	h_qvaltest_target->GetXaxis()->SetRangeUser(-4.4,-3.8);
	h_qvaltest_target->SetTitle("Q value with target");
	
	TCanvas *c4 = new TCanvas ( "c4" );
	c4->SetTitle("Same energy loss");
	c4->Divide(2,1);
	c4->cd(1);
	h_qval_notarget->Draw();
	h_qval_notarget->GetXaxis()->SetRangeUser(-4.4,-3.8);
	h_qval_notarget->SetTitle("Q value without target");
	c4->cd(2);
	h_qval_loss->Draw();
	h_qval_loss->GetXaxis()->SetRangeUser(-4.4,-3.8);
	h_qval_loss->SetTitle("Q value with target");
	
	TCanvas *c5 = new TCanvas ( "c5" ); //create a canvas
	c5->SetTitle("Same energy loss");
	c5->Divide ( 1,2 );
	c5->cd(1);
	TGraph *gr5 = new TGraphErrors(length2,QValue2,Counts2,0,error2);
	gr5->SetTitle("Q value without target");
	gr5->Draw("AP*");
	//gr1->Fit(fit_gauss2,"E","",QValue1[0],QValue1[length1]);
	c5->cd(2);
	TGraph *gr6 = new TGraphErrors(length1,QValue1,Counts1,0,error1);
	gr6->SetTitle("Q value with target");
	gr6->Draw("AP*");
	
	TCanvas *c6 = new TCanvas ( "c6" );
	c6->SetTitle("Same energy loss, subtracted");
	TGraph *gr7 = new TGraph(length2,QValue2,CountsSub);
	gr7->SetTitle("Background subtracted");
	gr7->Draw("AP*");
}