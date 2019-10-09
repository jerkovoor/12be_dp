#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCanvas.h"
using namespace std;

void lossfit(){
	int entries=22, entries1=9, entries2=9, entries3=10;
	//double endetB[16][entries1], lossB[16][entries1], endetAl[16][entries2], lossAl[16][entries2], endetD[16][entries3], lossD[16][entries3];
	double energy[16][entries], target[16][entries], tar_det[16][entries],  Al_layer[16][entries], Al_det[16][entries], B_layer[16][entries], B_det[16][entries], remained[16][entries];
	TGraph* f[16];
	TGraph* g[16];
	TGraph* h[16];


	/*for(int i=1;i<17;i++)
	{
		ifstream test1;
		string f_name1=Form("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/B_Dead_layer/B_%i.txt",i);
		test1.open(f_name1.c_str());
		if(test1.is_open())
		{
		for(int j=0;j<entries1;j++)
		{
		test1 >> lossB[i-1][j] >> endetB[i-1][j];
		//cout << lossB[i-1][j] << endl;
		}
		}
		test1.clear();
	}
	
	
	for(int i=1;i<17;i++)
	{
		ifstream test2;
		string f_name2=Form("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Al_Dead_layer/Al_%i.txt",i);
		test2.open(f_name2.c_str());
		if(test2.is_open())
		{
		for(int j=0;j<entries2;j++)
		{
		test2 >> lossAl[i-1][j] >> endetAl[i-1][j];
		}
		}
		test2.clear();
	}
	
	
	for(int i=1;i<17;i++)
	{
		ifstream test3;
		string f_name3=Form("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Target/D2_%i.txt",i);
		test3.open(f_name3.c_str());
		if(test3.is_open())
		{
		for(int j=0;j<entries3;j++)
		{
		test3 >> lossD[i-1][j] >> endetD[i-1][j];
		}
		}
		test3.clear();
	}*/
  
	for(int i=0;i<16;i++){
		ifstream test;
		string f_name=Form("/mnt/c/Users/Jerome/OneDrive/Work/12Be/elast/outputs/out_TRIUMF_%i_1H.dat",i);
		test.open(f_name.c_str());
		if(test.is_open()){
			for(int j=0;j<entries;j++){
				test >> energy[i][j] >> target[i][j] >> Al_layer[i][j] >> B_layer[i][j] >> B_det[i][j];
				tar_det[i][j]=energy[i][j]-target[i][j];
				Al_det[i][j]=energy[i][j]-target[i][j]-Al_layer[i][j];
				//cout << tar_det[i][j] << Al_det[i][j] << endl;
			}
		}
		test.clear();
	}
   
 
 
	TF1* f1 = new TF1 ( "f1","[0]+[1]*x^(-1)+[2]*x^(-2)",0,7); //[1]*x+[4]*TMath::Exp(-[5]*x)



	ofstream bor;
	bor.open ( "/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/B_Dead_layer/B_fitparameters_TRIUMF.txt" );//change the path and the name of the file acordingly
	bor << "Ring	a		b		c\n";
		
		
	TCanvas *c1 = new TCanvas ( "c1" ); //create a canvas
	c1->Divide ( 4,4 ); //devide the Canvas in 16
	
		
	for(int i=0;i<16;i++){
		c1->cd ( i+1 );
		//f[i] = new TGraph(entries,B_det[i],B_layer[i]);
		f[i] = new TGraph(entries,B_layer[i],B_det[i]);
		f[i]->Draw ( "AP*" );
		f[i]->Fit(f1,"","",0,6);
		bor << i << "	"  << f1->GetParameter ( 0 ) << "	" << f1->GetParameter ( 1 ) << "	"<< f1->GetParameter ( 2 ) << endl;
	}
	
	
	
	ofstream alum;
	alum.open ( "/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Al_Dead_layer/Al_fitparameters_TRIUMF.txt" );//change the path and the name of the file acordingly
	alum << "Ring	d		e		f\n";
	
	TCanvas *c2 = new TCanvas ( "c2" ); //create a canvas
	c2->Divide ( 4,4 ); //devide the Canvas in 16
		
		
	for(int i=0;i<16;i++){
		c2->cd ( i+1 );
		//g[i] = new TGraph(entries,Al_det[i],Al_layer[i]);
		g[i] = new TGraph(entries,Al_layer[i],Al_det[i]);
		g[i]->Draw ( "AP*" );
		g[i]->Fit(f1,"","",0,6);
		alum << i << "	"  << f1->GetParameter ( 0 ) << "	" << f1->GetParameter ( 1 ) << "	"<<f1->GetParameter ( 2 ) << endl;
	}

	
	
	ofstream tar;
	tar.open ( "/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Target/D2_fitparameters_TRIUMF.txt" );//change the path and the name of the file acordingly
	tar << "Ring	g		h		k\n";
	
	TCanvas *c3 = new TCanvas ( "c3" ); //create a canvas
	c3->Divide ( 4,4 ); //devide the Canvas in 16
		
		
	for(int i=0;i<16;i++){
		c3->cd ( i+1 );
		//h[i] = new TGraph(entries,tar_det[i],target[i]);
		h[i] = new TGraph(entries,target[i],tar_det[i]);
		h[i]->Draw ( "AP*" );
		h[i]->Fit(f1,"","",0,6);
		tar << i << "	"  << f1->GetParameter ( 0 ) << "	" << f1->GetParameter ( 1 ) << "	"<<f1->GetParameter ( 2 ) << endl;
	}
	

	/*mg->Add(g);
	mg->Draw("a");
	//mg->GetXaxis()->SetTitle("Effective Thickness (microns)");
	//mg->GetYaxis()->SetTitle("Energy loss (MeV/u)");
	//mg->GetXaxis()->SetLimits(0,1);
	c1->BuildLegend();*/
  
  
}