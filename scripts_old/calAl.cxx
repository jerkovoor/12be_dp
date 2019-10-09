//M. Vostinar
//The purpose of this code is to read in the positions of the pedestals as well as the alpha peaks
//and use that to create graphs which are further fitted with a linear fit 
//the output of the fit gives the gain and the offset necessary for the calibration of the corresponding detector
//the outputs are .txt and .root files

using namespace std;

#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "iostream"
#include "fstream"

double calAl(){
	
	TGraph* gr[128];
	
	//definition of variables for YY detectors; strip number, pedestal and 3 peaks positions
	double strip[128], peakX1[128], peakX2[128], peakX3[128], ped[128], strip2[128];
	//double strip[32], peakX1[32], ped[32], strip2[32];
	
	ifstream Alpha_peaks; //open the YY detector .txt file containing the positions of three alpha peaks
	Alpha_peaks.open("/home/jerome/12Be_exp/Alpha_peak_position/Yu_AlphaPeaks.txt"); //adjust the location and file name acordingly
		
	if(Alpha_peaks.is_open()){
		Alpha_peaks.ignore(256,'\n');
		for(int i=0;i<128;i++){ //adjust to the number coresponding to the detector you are looking at
			Alpha_peaks >> strip[i] >> peakX1[i] >> peakX2[i] >> peakX3[i];
			//Alpha_peaks >> strip[i] >> peakX1[i];
			cout << "Strip " <<  strip[i] << " " << peakX1[i] << " " << peakX2[i] << " " << peakX3[i] << endl;
			//cout << "Strip " <<  strip[i] << endl;
		}
	}
	
	else{
		cout << "No Alpha peaks file found " << endl;
	}
	Alpha_peaks.close();
 
	ifstream Pedestals; //open the YY detector .txt file containing the positions of pedestals
	Pedestals.open("/home/jerome/12Be_exp/Alpha_peak_position/4815Yu_Pedestals.txt");//adjust the location and file name acordingly
	if(Pedestals.is_open()){
		Pedestals.ignore(256,'\n'); //adjust to the number coresponding to the detector you are looking at
		for(int i=0;i<128;i++){
			Pedestals >> strip2[i] >> ped[i];
			cout << strip2[i] << " " << ped[i] << endl;
		}
	}
	
	else{
		cout << "No pedestal position file found " << endl;
	}
	
	Pedestals.close();
 
	//double Etheo[4]= {0,5156,5486,5805}; //definition of the array with known values
	double Etheo[3]= {5156,5486,5805}; // without pedestal
	//double Etheo[2]= {0,21423.2}; //21423.2 for Sd1	82977.6 for Sd2
 
 
	double alpha_loss1[16];
	ifstream loss1;
	loss1.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Alpha_calibration/alpha1.txt");
	if(loss1.is_open()){
		for(int i=0; i<16; i++){
			loss1 >> alpha_loss1[i];
		}
	}

	double alpha_loss2[16];
	ifstream loss2;
	loss2.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Alpha_calibration/alpha2.txt");
	if(loss2.is_open()){
		for(int i=0; i<16; i++){
			loss2 >> alpha_loss2[i];
		}
	}

	double alpha_loss3[16];
	ifstream loss3;
	loss3.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Alpha_calibration/alpha3.txt");
	if(loss3.is_open()){
		for(int i=0; i<16; i++){
			loss3 >> alpha_loss3[i];
		}
	}
 
	double alphaE1[128], alphaE2[128], alphaE3[128]; 
 
	for(int i=0; i<128;i++){
		for(int j=0; j<8; j++){
			for(int k=1; k<16; k++){
				if(i==j*16){
					alphaE1[i]=alpha_loss1[0];
				}
				else if (i==16*j+k){
					alphaE1[i]=alpha_loss1[k];
				}
			}
		}
		//cout << alphaE1[i] << endl;
	}
   
	for(int i=0; i<128;i++){
		for(int j=0; j<8; j++){
			for(int k=1; k<16; k++){
				if(i==j*16){
					alphaE2[i]=alpha_loss2[0];
				}
				else if(i==16*j+k){
					alphaE2[i]=alpha_loss2[k];
				}
			}
		}
	}
   
	for(int i=0; i<128;i++){
		for(int j=0; j<8; j++){
			for(int k=1; k<16; k++){
				if(i==j*16){
					alphaE3[i]=alpha_loss3[0];
				}
				else if(i==16*j+k){
					alphaE3[i]=alpha_loss3[k];
				}
			}
		}
	}
	
 
	for(int i=0; i<128;i++){ //loop over the strips of the detector
		double Eexp[4]={ped[i],peakX1[i],peakX2[i],peakX3[i]}; //create the array with the values read out from .txt files
		//double Eexp[3]={peakX1[i],peakX2[i],peakX3[i]};
		//double Eexp[2]={ped[i],peakX1[i]};
		//double Ethe[3]={1000*alphaE1[i],1000*alphaE2[i],1000*alphaE3[i]};
		double Ethe[4]={0,1000*alphaE1[i],1000*alphaE2[i],1000*alphaE3[i]};
		gr[i] = new TGraph(4,Eexp,Ethe); //fill in the graphs for each strip of the detector
	}
	
	
	//open a txt file to store the positions of the peaks
	ofstream Alpha_cal;
	Alpha_cal.open("/home/jerome/12Be_exp/Analysis/with_pedestal/Yu_pedestal_fit.txt"); //open a .txt file to store the results of the fit; change the path and name acordingly
	Alpha_cal << " strip / gain / offset  \n";
	
	TF1 *fit_lin = new TF1 ( "fit_lin","[0]*x+[1]",0,4000 ); //definition of a linear fuction in a range from 0 to 4000

	//Adjust the number of canvas according to the amount of strips you have in the corresponding detector
	
	
	
	
	
	/*TCanvas *c1 = new TCanvas ( "c1" ); //create a canvas
	c1->Divide ( 4,4 ); //devide the Canvas in 16
		
		
		for ( int i=0; i<16; i++ ){ //loop over the first 16 graphs
			c1->cd ( i+1 );
			gr[i]->Draw ( "AP*" );
			gr[i]->Fit ( "fit_lin","","",0,4000 ); //fir the first 16 graphs
			Alpha_cal << i << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl; //fill in the .txt file
		}

		TCanvas *c2 = new TCanvas ( "c2" ); //create a new canvas
		c2->Divide ( 4,4 ); //devide in 16
		for ( int i=0; i<16; i++ ){ //loop over them
			c2->cd ( i+1 );
			gr[i+16]->Draw ( "AP*" );
			gr[i+16]->Fit ( "fit_lin","","",0,4000 );
			Alpha_cal << i+16 << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;  
		}

	
	
	TCanvas *c3 = new TCanvas ( "c3" );
		c3->Divide ( 4,3 );
		for ( int i=0; i<12; i++ )
		{
			c3->cd ( i+1 );
			gr[i]->Draw ( "AP*" );
			gr[i]->Fit ( "fit_lin","","",0,4000 );
			Alpha_cal << i << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;  
		}

		TCanvas *c4 = new TCanvas ( "c4" );
		c4->Divide ( 4,3 );
		for ( int i=0; i<12; i++ )
		{
			c4->cd ( i+1 );
			gr[i+12]->Draw ( "AP*" );
			gr[i+12]->Fit ( "fit_lin","","",0,4000 );
			Alpha_cal << i+12 << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;  
		}*/
		
	TCanvas *c1 = new TCanvas ( "c1" ); //create a canvas
	c1->Divide ( 4,4 ); //devide the Canvas in 16
		
		
		for ( int i=0; i<16; i++ ) //loop over the first 16 graphs
		{
			c1->cd ( i+1 );
			gr[i]->Draw ( "AP*" );
			gr[i]->Fit ( "fit_lin","","",0,4000 ); //fir the first 16 graphs
			Alpha_cal << i << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl; //fill in the .txt file
		
		}

		TCanvas *c2 = new TCanvas ( "c2" ); //create a new canvas
		c2->Divide ( 4,4 ); //devide in 16
		for ( int i=0; i<16; i++ ) //loop over them
		{
			c2->cd ( i+1 );
			gr[i+16]->Draw ( "AP*" );
			gr[i+16]->Fit ( "fit_lin","","",0,4000 );
			Alpha_cal << i+16 << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;  
		}

		TCanvas *c3 = new TCanvas ( "c3" );
		c3->Divide ( 4,4 );
		for ( int i=0; i<16; i++ )
		{
			c3->cd ( i+1 );
			gr[i+32]->Draw ( "AP*" );
			gr[i+32]->Fit ( "fit_lin","","",0,4000 );
			Alpha_cal << i+32 << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;  
		}

		TCanvas *c4 = new TCanvas ( "c4" );
		c4->Divide ( 4,4 );
		for ( int i=0; i<16; i++ )
		{
			c4->cd ( i+1 );
			gr[i+48]->Draw ( "AP*" );
			gr[i+48]->Fit ( "fit_lin","","",0,4000 );
			Alpha_cal << i+48 << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;  
		}

		TCanvas *c5 = new TCanvas ( "c5" );
		c5->Divide ( 4,4 );
		for ( int i=0; i<16; i++ )
		{
			c5->cd ( i+1 );
			gr[i+64]->Draw ( "AP*" );
			gr[i+64]->Fit ( "fit_lin","","",0,4000 );
			Alpha_cal << i+64 << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;  
		}

		TCanvas *c6 = new TCanvas ( "c6" );
		c6->Divide ( 4,4 );
		for ( int i=0; i<16; i++ )
		{
			c6->cd ( i+1 );
			gr[i+80]->Draw ( "AP*" );
			gr[i+80]->Fit ( "fit_lin","","",0,4000 );
			Alpha_cal << i+80 << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;  
		}

		TCanvas *c7 = new TCanvas ( "c7" );
		c7->Divide ( 4,4 );
		for ( int i=0; i<16; i++ )
		{
			c7->cd ( i+1 );
			gr[i+96]->Draw ( "AP*" );
			gr[i+96]->Fit ( "fit_lin","","",0,4000 );
			Alpha_cal << i+96 << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;  
		}

		TCanvas *c8 = new TCanvas ( "c8" );
		c8->Divide ( 4,4 );
		for ( int i=0; i<16; i++ )
		{
			c8->cd ( i+1 );
			gr[i+112]->Draw ( "AP*" );
			gr[i+112]->Fit ( "fit_lin","","",0,4000 );
			Alpha_cal << i+112 << " "  << fit_lin->GetParameter ( 0 ) << " " << fit_lin->GetParameter ( 1 ) << endl;  
		}
		
		
		TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/with_pedestal/Yu_pedestal_fit.root","RECREATE"); //create an output .root file; change the path and name acordingly
		
		f_out->cd();
		
		for ( int i=0; i<128; i++ ) //adjust the sieze of the loop according to the number of strips in the interesting detector
	{
		string name2 = Form ( "4points_%i",i );
		gr[i]->Write(name2.c_str());
	}
		
		f_out->Close(); //close the output file
	
	return 0;
	} //end of script