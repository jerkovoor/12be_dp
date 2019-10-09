//M. Vostinar
//The code takes the files which were already treated with decode.cxx
//It is used for calibration purposes
//Fils in the corresponding histograms and uses TSpectrum to search for peaks
//After TSpectrum has found the peaks, depending on the judgment of person fiting one can choose a 3,6 or 9 gauss fit to use on the histogrames
//The histogrames with the fit are saved in a .root file and 2 txt. files contain the position of the tree interesting peaks and the amplitude of these peaks
//The code can be used on all detectors

using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TSpectrum.h"
#include "iostream"
#include "fstream"
#include "TH1.h"
#include "TF1.h"


double readAlpha(){
	//definition of variables necessary for reading the input root file
	//YYd detector
	int YdMul;
	vector<int>* YdChannel = new vector<int>;
	vector<double>* YdEnergyRaw= new vector<double>;
	vector<int>* YdRing= new vector<int>;
	vector<int>* YdSector= new vector<int>;
	
	//YYU detector
	int YuMul;
	vector<int>* YuChannel= new vector<int>;
	vector<double>* YuEnergyRaw= new vector<double>;
	vector<int>* YuRing= new vector<int>;
	vector<int>* YuSector= new vector<int>;
		
	//Purpose of two readouts is to increase the range of operation (multiply it by 2 in total)
	//CsI 1 detector
	int CsI1Mul;
	vector<int>* CsI1Channel= new vector<int>; //Channel corresponds to one christal - 1 means that the gain is set to one value
	vector<double>* CsI1EnergyRaw= new vector<double>;

	//CsI 2 detector
	int CsI2Mul;
	vector<int>* CsI2Channel= new vector<int>;//Channel corresponds to one christal - 2 means that the gain is set to a different value than 1
	vector<double>* CsI2EnergyRaw= new vector<double>;

	//IC chamber
	int ICChannel;
	double ICEnergyRaw;

	//Sdd1 detector
	int Sd1rMul;
	vector<int>* Sd1rChannel= new vector<int>; //Ring Channels
	vector<double>* Sd1rEnergyRaw= new vector<double>;
	
	//Sd1s detector
	int Sd1sMul;
	vector<int>* Sd1sChannel= new vector<int>; //Sector Channels
	vector<double>* Sd1sEnergyRaw= new vector<double>;
	
	//Sd2r detector
	int Sd2rMul;
	vector<int>* Sd2rChannel= new vector<int>; //Ring Channels
	vector<double>* Sd2rEnergyRaw= new vector<double>;
	vector<double>* Sd2rEnergy= new vector<double>;
	
	//Sd2s detector
	int Sd2sMul;
	vector<int>* Sd2sChannel= new vector<int>; //Sector Channels
	vector<double>* Sd2sEnergyRaw= new vector<double>;
	vector<double>* Sd2sEnergy= new vector<double>;

	//Sur
	int SurMul;
	vector<int>* SurChannel= new vector<int>; //Ring Channels
	vector<double>* SurEnergyRaw= new vector<double>;

	//Sus
	int SusMul;
	vector<int>* SusChannel= new vector<int>; //Sector Channels
	vector<double>* SusEnergyRaw= new vector<double>;
	
	TChain* chain = new TChain ( "AutoTree" );

	//Open the calibration files 
	// Yd calibration before the experiment
	//chain->Add ( "/home/jerome/12Be_exp/Processed_files/C_notarget/decode4992.root" ); //change the path to the files and the file name accordingly
	/*chain->Add ( "/home/jerome/12Be_exp/Analysis/Su_calibration/decode_4972.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/Su_calibration/decode_4973.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/Su_calibration/decode_4974.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/Su_calibration/decode_4975.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/Su_calibration/decode_4977.root" );*/
	
	chain->Add ( "/home/jerome/12Be_exp/Analysis/Su_calibration/decode_5225_1.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/Su_calibration/decode_5225_2.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/Su_calibration/decode_5225_3.root" );

	//the 5225 Yu calibration run
	/* chain->Add("/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Alpha/Decode_5225_1.root");
	chain->Add("/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Alpha/Decode_5225_2.root");
	chain->Add("/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Alpha/Decode_5225_3.root");*/

	//read out the input tree
	chain->SetBranchAddress ( "YdMul",&YdMul );
	chain->SetBranchAddress ( "YdChannel",&YdChannel );
	chain->SetBranchAddress ( "YdEnergyRaw",&YdEnergyRaw );
	chain->SetBranchAddress ( "YdRing",&YdRing );
	chain->SetBranchAddress ( "YdSector",&YdSector );

	chain->SetBranchAddress ( "YuMul",&YuMul );
	chain->SetBranchAddress ( "YuChannel",&YuChannel );
	chain->SetBranchAddress ( "YuEnergyRaw",&YuEnergyRaw );
	chain->SetBranchAddress ( "YuRing",&YuRing );
	chain->SetBranchAddress ( "YuSector",&YuSector );
	
	chain->SetBranchAddress ( "CsI1Mul",&CsI1Mul );
	chain->SetBranchAddress ( "CsI1Channel",&CsI1Channel );
	chain->SetBranchAddress ( "CsI1EnergyRaw",&CsI1EnergyRaw );

	chain->SetBranchAddress ( "CsI2Mul",&CsI2Mul );
	chain->SetBranchAddress ( "CsI2Channel",&CsI2Channel );
	chain->SetBranchAddress ( "CsI2EnergyRaw",&CsI2EnergyRaw );
	
	chain->SetBranchAddress ( "ICChannel",&ICChannel );
	chain->SetBranchAddress ( "ICEnergyRaw",&ICEnergyRaw );

	chain->SetBranchAddress ( "Sd1rMul",&Sd1rMul );
	chain->SetBranchAddress ( "Sd1rChannel",&Sd1rChannel );
	chain->SetBranchAddress ( "Sd1rEnergyRaw",&Sd1rEnergyRaw );

	chain->SetBranchAddress ( "Sd1sMul",&Sd1sMul );
	chain->SetBranchAddress ( "Sd1sChannel",&Sd1sChannel );
	chain->SetBranchAddress ( "Sd1sEnergyRaw",&Sd1sEnergyRaw );

	chain->SetBranchAddress ( "Sd2rMul",&Sd2rMul );
	chain->SetBranchAddress ( "Sd2rChannel",&Sd2rChannel );
	chain->SetBranchAddress ( "Sd2rEnergyRaw",&Sd2rEnergyRaw );
	chain->SetBranchAddress ( "Sd2rEnergy",&Sd2rEnergy );

	chain->SetBranchAddress ( "Sd2sMul",&Sd2sMul );
	chain->SetBranchAddress ( "Sd2sChannel",&Sd2sChannel );
	chain->SetBranchAddress ( "Sd2sEnergyRaw",&Sd2sEnergyRaw );
	chain->SetBranchAddress ( "Sd2sEnergy",&Sd2sEnergy );
	
	chain->SetBranchAddress ( "SurMul",&SurMul );
	chain->SetBranchAddress ( "SurChannel",&SurChannel );
	chain->SetBranchAddress ( "SurEnergyRaw",&SurEnergyRaw );
	
	chain->SetBranchAddress ( "SusMul",&SusMul );
	chain->SetBranchAddress ( "SusChannel",&SusChannel );
	chain->SetBranchAddress ( "SusEnergyRaw",&SusEnergyRaw );
		
	//definition of various histogrames
	//one histograme is created for each strip to be able to perform the peak search for calibration
	TH1D *hYd[128];
	TH1D *hYu[128];
	TH1D* hCsI1[16];
	TH1D* hCsI2[16];
	TH1D* hSd1r[24];
	TH1D* hSd2r[24];
	TH1D* hSd1s[32];
	TH1D* hSd2s[32];
	TH1D* hSur[24];
	TH1D* hSus[32];

	for ( int i=0; i<128; i++ ){
		string Yd_name = Form ( "Yd_Ch%i",i );
		hYd[i] =  new TH1D ( Yd_name.c_str(),"Yd_Ch",1024,0,4096 );

		string Yu_name = Form ( "Yu_Ch%i",i );
		hYu[i] =  new TH1D ( Yu_name.c_str(),"Yu_Ch",175,2400,3800); //modified to see if the calib is better 

	}
		
	TH1D *hYuSectorRing[8][16]; //These histogrames indicate the sector and the ring of the Yu detector. The purpose is to investigate trends conected to sectors/rings
		
	for(int i=0; i<8; i++){
		for(int j=0; j<16; j++){
				string nameSR = Form("Yu_Sec%i_Ring%i", i,j);
				hYuSectorRing[i][j] = new TH1D ( nameSR.c_str(),"Yu", 350,2400,3800 );
		}
	}
		
	for(int i=0; i<16; i++){
		string CsI1_name = Form("CsI1_%i",i); 
		hCsI1[i] = new TH1D (CsI1_name.c_str(),"CsI1",1024,0,4096);
		
		string CsI2_name = Form("CsI2_%i",i); 
		hCsI2[i] = new TH1D (CsI2_name.c_str(),"CsI2",1024,0,4096);
		
	}
		
	for(int i=0; i<24; i++){
		string Sd1r_name = Form ( "Sd1r_Ch%i",i );
		hSd1r[i] =  new TH1D ( Sd1r_name.c_str(),"Sd1r_Ch",1024,0,4096);

		string Sd2r_name = Form ( "Sd2r_Ch%i",i );
		hSd2r[i] =  new TH1D ( Sd2r_name.c_str(),"Sd2r_Ch",1024,0,4096 );

		string Sur_name = Form ( "Sur_Ch%i",i );
		hSur[i] =  new TH1D ( Sur_name.c_str(),"Sur_Ch",1024,0,4096 );
	}
		
	for(int i=0; i<32; i++){
		string Sd1s_name = Form ( "Sd1s_Ch%i",i );
		hSd1s[i] =  new TH1D ( Sd1s_name.c_str(),"Sd1s_Ch",1024,0,4096 );

		string Sd2s_name = Form ( "Sd2s_Ch%i",i );
		hSd2s[i] =  new TH1D ( Sd2s_name.c_str(),"Sd2s_Ch",1024,0,4096 );

		string Sus_name = Form ( "Sus_Ch%i",i );
		hSus[i] =  new TH1D ( Sus_name.c_str(),"Sus_Ch",1024,0,4096 );
	}
		
	TFile *f_out = new TFile ( "/home/jerome/12Be_exp/Analysis/Su_calibration/Sur_AlphaCalibration.root","RECREATE" ); //Creating the output .root file put the corect path and/or name
		
	int ev = chain->GetEntries(); //get the total number of entries
	cout << "Total number of events =" << ev << endl;
	int ev_num=0;
		
	while ( ev_num<=ev ){ //the main while loop, loops over all the events in the tree and does whatever needs to be done 
		
		if ( ev_num%10000==0 ) cout << "Current event = " << ev_num << "\r"<< flush;
		//cout << "in the while loop" << endl;
		chain->GetEntry ( ev_num ); //get the current event
		//the filing of the histograms is performed in the same way for all different detectors, requireing that a channel should or shouldn't be >500 to evade confusion of the TSpectrum 
		if ( YdEnergyRaw->size() >0){
			if ( YdEnergyRaw->at(0)>500 ){
				hYd[YdChannel->at(0)]->Fill ( YdEnergyRaw->at(0) );
			}
		} //requires that the vector is acctually filled and that the channel is >500 so not to take into acount the pedestal, and estimate only the position of the alpha peaks
		
		if ( YuEnergyRaw->size() >0){
			if(YuEnergyRaw->at(0)>0){
				hYu[YuChannel->at ( 0 )]->Fill ( YuEnergyRaw->at ( 0 ) );
				hYuSectorRing[YuSector->at(0)][YuRing->at(0)]->Fill(YuEnergyRaw->at ( 0 ) ); 
			}
		}
			
		//if ( CsI1EnergyRaw->size() >0 && CsI1Mul==1 ){if ( CsI1EnergyRaw->at ( 0 ) >0 ){hCsI1[CsI1Channel->at ( 0 )]->Fill ( CsI1EnergyRaw->at ( 0 ) );}} //CsI is not calibrated using the Alpha runs
		//if ( CsI2EnergyRaw->size() >0 && CsI2Mul==1 ){if ( CsI2EnergyRaw->at ( 0 ) >0 ){hCsI2[CsI2Channel->at ( 0 )]->Fill ( CsI2EnergyRaw->at ( 0 ) );}}
		if ( Sd1rEnergyRaw->size() >0){
			if ( Sd1rEnergyRaw->at ( 0 ) >0 ){
				hSd1r[Sd1rChannel->at ( 0 )]->Fill ( Sd1rEnergyRaw->at ( 0 ) );
			}
		}
		
		if ( Sd1sEnergyRaw->size() >0){
			if ( Sd1sEnergyRaw->at ( 0 ) >0 ){
				hSd1s[Sd1sChannel->at ( 0 )]->Fill ( Sd1sEnergyRaw->at ( 0 ) );
			}
		}
		
		if ( Sd2rEnergyRaw->size() >0){
			if ( Sd2rEnergyRaw->at ( 0 ) >0 ){
				hSd2r[Sd2rChannel->at ( 0 )]->Fill ( Sd2rEnergyRaw->at ( 0 ) );
			}
		}
		
		if ( Sd2sEnergyRaw->size() >0){
			if ( Sd2sEnergyRaw->at ( 0 ) >0 ){
				hSd2s[Sd2sChannel->at ( 0 )]->Fill ( Sd2sEnergyRaw->at ( 0 ) );
			}
		}
		
		if ( SurEnergyRaw->size() >0){
			if ( SurEnergyRaw->at ( 0 ) >500 ){
				hSur[SurChannel->at ( 0 )]->Fill ( SurEnergyRaw->at ( 0 ) );
			}
		}
		
		if ( SusEnergyRaw->size() >0){
			if ( SusEnergyRaw->at ( 0 ) >500 ){
				hSus[SusChannel->at ( 0 )]->Fill ( SusEnergyRaw->at ( 0 ) );
			}
		}
		
		ev_num++;
	}//end of the while loop
	
	//definition of convoluted gaussian functions
	//1 gaussian
	TF1* fit_func1 = new TF1 ( "fit_func1","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",0,5000);
	TF1* fit_func2 = new TF1 ( "fit_func2","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",0,5000);
	//6 gaussians 
	TF1* fit_func6 = new TF1 ( "fit_func6","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp(-pow((x-[10]),2)/(2*[2]*[2]))+[11]*TMath::Exp(-pow((x-[12]),2)/(2*[2]*[2]))",500,5000 ); //removed the liner part for comparison
	//3 gaussians
	TF1* fit_func3 = new TF1 ( "fit_func3","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*x+[8]",0,5000 );
	//9 gaussians with linear background
	TF1* fit_func9 = new TF1 ( "fit_func9","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp(-pow((x-[10]),2)/(2*[2]*[2]))+[11]*TMath::Exp(-pow((x-[12]),2)/(2*[2]*[2]))+[13]*TMath::Exp(-pow((x-[14]),2)/(2*[2]*[2]))+[15]*TMath::Exp(-pow((x-[16]),2)/(2*[2]*[2]))+[17]*TMath::Exp(-pow((x-[18]),2)/(2*[2]*[2]))+[19]+[20]*x",2000,5000 );
	
	//Use TSpectrum to find the position of the peaks and limit the fit parameters to extract the exact positions of max for each peak
	//Search for peaks in Yd strips
	int npeaks=1; //number of peaks expected to see
	TSpectrum* s = new TSpectrum ( 2*npeaks );
	int nfound=0;
	double* txp = s->GetPositionX(); //returns the channel in which the peak has a max value
    
	////////////////////////////////////
	////		SU DETECTORS		////
	////////////////////////////////////
	
	
	//find peaks and for the Sdr detectors
	double Surxp[3];
	double Sura[24],Surb[24],Surc[24];
	
	
	for ( int j=0; j<24; j++ ){
		nfound = s->Search ( hSur[j],2,"",0.1 ); //search for peaks in the hSdr detectors change the name of the histograme according to which peaks you want to search for Sd1r, Sd2r or Sur
       
		for ( int p=0; p<nfound; p++ ){Surxp[p] = txp[p];} 
        
		Sura[j] = min ( min ( Surxp[0],Surxp[1] ), Surxp[2] );
		Surc[j] = max ( max ( Surxp[0],Surxp[1] ), Surxp[2] );

		if(Surxp[0]!=Sura[j] && Surxp[0]!=Surc[j]) {Surb[j] = Surxp[0];} 
		else if(Surxp[1]!=Sura[j] && Surxp[1]!=Surc[j]) {Surb[j] = Surxp[1];}
		else if(Surxp[2]!=Sura[j] && Surxp[2]!=Surc[j]) {Surb[j] = Surxp[2];}  
    }
  
	//open a txt file to store the positions of the peaks
	ofstream Sur_alpha;
	Sur_alpha.open ( "/home/jerome/12Be_exp/Analysis/Su_calibration/Sur_AlphaCalibration_peak.txt" ); //change the path and the name of the file acordingly
	Sur_alpha << " Strip / Peak 1 / 2 / 3  \n";
	//Sur_alpha << " Strip / Peak 1 \n";
	
	ofstream Sur_alphaCounts;
	Sur_alphaCounts.open ( "/home/jerome/12Be_exp/Analysis/Su_calibration/Sur_AlphaCalibration_counts.txt" ); //change the path and the name of the file acordingly
	Sur_alphaCounts << " Strip / Amplitude 1 / 2 / 3  \n";
	//Sur_alphaCounts << " Strip / Amplitude 1 \n";
	
	for ( int i=0; i<24; i++ ){ //loop over all the rings of S3 detectors to set the parameters and fit
		
		//if(i<6){
			fit_func3->SetParLimits ( 1,Sura[i]-5,Sura[i]+5 ); //remove here
			fit_func3->SetParLimits ( 4,Surb[i]-5,Surb[i]+5 );
			fit_func3->SetParLimits ( 6,Surc[i]-5,Surc[i]+5 );
			fit_func3->SetParLimits ( 2,0,10 );			//remove here
			
			//fit_func1->SetParLimits ( 1,Surxp[i]-200,Surxp[i]+200 ); //2000,2500
			//fit_func1->SetParLimits ( 2,0,200 ); //0,200
	//	}else if(i>5){
			//fit_func3->SetParLimits ( 1,Sura[i]-3,Sura[i]+3 );
			//fit_func3->SetParLimits ( 4,Surb[i]-3,Surb[i]+3 );
			//fit_func3->SetParLimits ( 6,Surc[i]-3,Surc[i]+3 );
			//fit_func3->SetParLimits ( 2,0,3 );
		//}
		hSur[i]->Fit ( "fit_func3","","",Surxp[i]-500,Surxp[i]+500 ); //change the histograme accordingly // 1000,3000

		//fill in the .txt files
		Sur_alpha << i << " "  << fit_func3->GetParameter ( 1 ) << " " << fit_func3->GetParameter ( 4 ) << " " << fit_func3->GetParameter ( 6 ) << endl;
		//Sur_alpha << i << " "  << fit_func1->GetParameter ( 1 ) << endl;
		Sur_alphaCounts << i << " " << fit_func3->GetParameter ( 0 ) << " " <<  fit_func3->GetParameter ( 3 ) << " " <<  fit_func3->GetParameter ( 5 ) << endl;
		//Sur_alphaCounts << i << " " << fit_func1->GetParameter ( 0 ) << endl;
	}//end of fit on Sur detector
    
    
  
	//find peaks and for the Sus detectors
	/*double Susxp[1];
	double Susa[32],Susb[32],Susc[32],Susped[32];
	
	for ( int j=0; j<32; j++ ){ //loop over the sectors of the S3 detectors
		nfound = s->Search ( hSus[j],2,"",0.1 ); //change the name of the hitogram according to which detector you want to look at: Sd1s, Sd2s or Sus
		
		for ( int p=0; p<nfound; p++ ){Susxp[p] = txp[p];}
		
		Susa[j] = min ( min ( Susxp[0],Susxp[1] ), Susxp[2] );
		Susc[j] = max ( max ( Susxp[0],Susxp[1] ), Susxp[2] );

		if(Susxp[0]!=Susa[j] && Susxp[0]!=Susc[j]) {Susb[j] = Susxp[0];}
		else if(Susxp[1]!=Susa[j] && Susxp[1]!=Susc[j]) {Susb[j] = Susxp[1];}
		else if(Susxp[2]!=Susa[j] && Susxp[2]!=Susc[j]) {Susb[j] = Susxp[2];}
                           
	}
  
	//open a txt file to store the positions of the peaks
	ofstream Sus_alpha;
	Sus_alpha.open ( "/home/jerome/12Be_exp/Analysis/Su_calibration/Sus_AlphaCalibration_peak.txt" );//change the path and the name of the file acordingly
	Sus_alpha << " Strip / Peak 1 / 2 / 3  \n";
	//Sus_alpha << " Strip / Peak 1 \n";
  
	ofstream Sus_alphaCounts;
	Sus_alphaCounts.open ( "/home/jerome/12Be_exp/Analysis/Su_calibration/Sus_AlphaCalibration_counts.txt" );//change the path and the name of the file acordingly
	Sus_alphaCounts << " Strip / Amplitude 1 / 2 / 3  \n";
	//Sus_alphaCounts <<  " Strip / Amplitude 1 \n";
    
	for ( int i=0; i<32; i++ ){ //loop over the sectors of S3 detectors to see the parameters and fit
		//if(i<16){
			fit_func3->SetParLimits ( 1,Susa[i]-5,Susa[i]+5 );
			fit_func3->SetParLimits ( 4,Susb[i]-5,Susb[i]+5 );
			fit_func3->SetParLimits ( 6,Susc[i]-5,Susc[i]+5 );
			fit_func3->SetParLimits ( 2,0,10 );
			//fit_func2->SetParLimits ( 1,Susxp[i]-200,Susxp[i]+200 );
			//fit_func2->SetParLimits ( 2,0,200 ); 
		//}else if(i>15){
			//fit_func3->SetParLimits ( 1,Susa[i]-3,Susa[i]+3 );
			//fit_func3->SetParLimits ( 4,Susb[i]-3,Susb[i]+3 );
			//fit_func3->SetParLimits ( 6,Susc[i]-3,Susc[i]+3 );
			//fit_func3->SetParLimits ( 2,0,3 );
		//}
  
		hSd1s[i]->Fit ( "fit_func2","","",Susxp[i]-500,Susxp[i]+500 ); //change the histograme name acordingly

		//fill in the .txt files
		Sus_alpha << i << " "  << fit_func3->GetParameter ( 1 ) << " " << fit_func3->GetParameter ( 4 ) << " " << fit_func3->GetParameter ( 6 ) << endl;
		//Sus_alpha << i << " "  << fit_func2->GetParameter ( 1 ) << endl;
		Sus_alphaCounts << i << " " << fit_func3->GetParameter ( 0 ) << " " <<  fit_func3->GetParameter ( 3 ) << " " <<  fit_func3->GetParameter ( 5 ) << endl;
		//Sus_alphaCounts << i << " " << fit_func2->GetParameter ( 0 ) <<  endl;
	}//end of fit on Sd1s detector
  */
	f_out->cd();

	//Turn on the histogrames you want to write in the .root file
	//for ( int i=0; i<128; i++ ){hYd[i]->Write();}
	//for ( int i=0; i<128; i++ ){hYu[i]->Write();}
	//for(int i=0; i<8; i++){for(int j=0; j<16; j++){hYuSectorRing[i][j]->Write();}}
	
	//for ( int i=0; i<16; i++ ){hCsI1[i]->Write();}
	//for ( int i=0; i<16; i++ ){hCsI2[i]->Write();}
	
	//for ( int i=0; i<24; i++ ){hSd1r[i]->Write();}
	//for ( int i=0; i<24; i++ ){hSd2r[i]->Write();}
	
	//for ( int i=0; i<32; i++ ){hSd1s[i]->Write();}
	//for ( int i=0; i<32; i++ ){hSd2s[i]->Write();}
	
	for ( int i=0; i<24; i++ ){hSur[i]->Write();}
	for ( int i=0; i<32; i++ ){hSus[i]->Write();}
	

	f_out->Close();

	return 0;
} //end of readAlpha.cxx