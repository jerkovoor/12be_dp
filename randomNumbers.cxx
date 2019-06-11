using namespace std;

#include "TFile.h" 
#include "TMath.h"
#include <cmath>
#include "iostream"
#include "fstream"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <TRandom3.h>

void randomNumbers(){
	// Open output file

	TFile* file = new TFile("RandomNumbers.root", "recreate");

	// Book histograms

	TH1D* h_Uni = new TH1D("h_Uni", "uniform random numbers",  100,  0, 1.0);
	TH1D* h_Exp = new TH1D("h_Exp", "exponential random numbers",  100,  0, 5.0);
	
	// Create a TRandom3 object to generate random numbers

	int seed = 123;
	TRandom3* ran = new TRandom3(seed);

	// Generate some random numbers and fill histograms

	const int numValues = 10000;
	const double xi = 1.0;                // mean of exponential pdf
	double Yubins[17]={0};
	double YuTheta[16];
	
	for (int i=0; i<17; i++){
		Yubins[i]=180-TMath::RadToDeg()*TMath::ATan((50+(16-i)*4.94)/(85));
	}
	
	for (int i=0; i<16; i++){
		YuTheta[i]=ran->Uniform(Yubins[i],Yubins[i+1]);
		
	}
	
	for(int i=0;i<10;i++){
		cout << ran->Uniform(Yubins[1],Yubins[2]) << endl;
	}
	
	/*for (int i=0; i<32; i++){
		cout << YuTheta[i%16] << endl;
	}*/


	for (int i=0; i<numValues; ++i){
		double r = ran->Uniform(0,2);             // uniform in ]0,1]
		double x = - xi * log(r);
		h_Uni->Fill(r);
		h_Exp->Fill(x);
	}

	// Store all histograms in the output file and close up
	//h_Uni->Draw();
	//h_Exp->Draw();
	file->Write();
	file->Close();
	

	return 0;
}
