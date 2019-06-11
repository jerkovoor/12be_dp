//M. Vostinar
using namespace std;

#include <TFile.h>
#include <TTree.h>
#include "/home/jerome/Software/TRIUMF/treeIris/include/IDet.h" //the treeIris library dependancy - You need to change the path to IDet.h depending on your computer!
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>

//double decode ( int run_num, int index)
double decodeYu ( int run_num ){
	
	//string f_name = Form ( "/home/jerome/12Be_exp/Raw_root/tree%i_%d.root",run_num,index);   
	string f_name = Form ( "/home/jerome/12Be_exp/Raw_root/12Be/tree%i.root",run_num ); //open the raw root file you want to convert to the format used by my codes; put the correct path and file name
	TFile* f_in = TFile::Open ( f_name.c_str(),"READ" );
	TTree* tr_in = ( TTree* ) f_in->Get ( "Iris" ); //reading the Iris tree; 
	//only the ADC leafs of the tree are used in the code. The rest of the data has a vector-array configuration or something like that and presents reliability which introduces unnecessary uncertanties
	
	IDet* det = new IDet();
	
	// Yu Sector and Ring Map
	std::map<int, std::pair<int, int> > YuSectorRingMap;
	for(int i = 0; i < 128; i++) {
		int sector = i/16;
		int ring = i%16;
		YuSectorRingMap[i] = std::make_pair(sector, ring);
	}
	

	Int_t* TYuMul= & ( det->TYuMul );
	std::vector<Int_t>* TYuChannel= & ( det->TYuChannel );
	std::vector<Double_t>* TYuEnergy= & ( det->TYuEnergy );
	std::vector<Int_t>* TYuADC= & ( det->TYuADC );
	std::vector<Int_t>* TYuNo= & ( det->TYuNo );
	std::vector<Int_t>* TYuRing= & ( det->TYuRing );
	std::vector<Double_t>* TYuTheta= & ( det->TYuTheta ); // Yd theta angle




	tr_in->SetBranchAddress ( "det", &det );

	//output variables
	
	//YYU detector
	int YuMul;
	vector<int> YuChannel;
	vector<double> YuEnergyRaw;
	vector<double> YuEnergy;
	vector<int> YuRing;
	vector<int> YuSector;
	/*std::vector<double> YuHitRawEnergy;
	std::vector<double> YuHitEnergy;
	std::vector<int> YuHitChannel;
	std::vector<int> YuHitSector;
	std::vector<int> YuHitRing;*/
	

	//open the output file
	
	string fOut_name = Form ( "/home/jerome/12Be_exp/Analysis/Be/CarbCalibAlpha/decodeBe_noYupedestal_1_2Peaks_CCalib_test%i.root",run_num ); //Open the output file; put the correct path and file name
	TFile* f_out = new TFile ( fOut_name.c_str(),"RECREATE" );
	TTree* tr_out = new TTree ( "AutoTree","AutoTree" ); //create the new root tree
	//Branches created in the new file; the leafs are removed from the three as well as arrays of vectors or whatever that was

	tr_out->Branch ( "YuMul",&YuMul );
	tr_out->Branch ( "YuChannel",&YuChannel );
	tr_out->Branch ( "YuEnergyRaw",&YuEnergyRaw );
	tr_out->Branch ( "YuEnergy",&YuEnergy );
	tr_out->Branch ( "YuRing",&YuRing );
	tr_out->Branch ( "YuSector",&YuSector );
	

	int ev = tr_in->GetEntries(); //get the total number of events in the tree
	cout << "Total number of events =" << ev << endl;
	int ev_num=0; 

	//open the calibration files

	//Opening and reading the Yu calibration file, same as for Yd
	double stripYu[128], Yup0[128],Yup1[128], Yup2[128];
	ifstream Alpha_peaksYu;
	//Alpha_peaksYu.open ( "/home/jerome/12Be_exp/Analysis/C_calib/Yu_C_calib_nopedestal_1_2Peaksfit.txt" );
	Alpha_peaksYu.open ( "/home/jerome/12Be_exp/Analysis/C_calib/Yu_C_calib_nopedestal_1_2Peaksfit.txt" ); //put the correct path and file name
	if ( Alpha_peaksYu.is_open() ){
		Alpha_peaksYu.ignore(256,'\n');
		for ( int i=0; i<128; i++ ){
			Alpha_peaksYu >> stripYu[i] >> Yup0[i] >> Yup1[i] ;
			//   cout << " Strip " <<  stripYu[i] << " P0 " << Yup0[i] << " p1 " << Yup1[i] << endl;
		}
	}else{
		cout << " No calib file for Yu detector "<< endl;
	}
	Alpha_peaksYu.close();
	

	//end of calibration filles

	while( ev_num<=ev ){ //entering the main while loop and looping over all the events in the Iris tree
		if ( ev_num%10000==0 ) cout << "Current event = " << ev_num << "\r"<< flush;
		tr_in->GetEntry ( ev_num ); //pulling up the current event from the Iris tree
		
		

		//////////////////////////////////////////////////////////
		//						YU DETECTOR						//
		//Treat the Yu detector, same as all the previous ones	//
		//////////////////////////////////////////////////////////

		if(det->TYuADC.size()>0){
			YuMul=0;
			for ( unsigned int adci=0; adci < det->TYuADC.size(); adci++ ){
				if ( det->TYuADC.at ( adci ) >0 ){
					YuEnergyRaw.push_back(det->TYuADC.at(adci));
					YuChannel.push_back(adci);
					YuMul++;
				}
			}
						
			//  cout << " Entered the Yu loop " << endl;
			
			for(size_t i = 0; i < YuEnergyRaw.size(); i++) {
				YuEnergy.push_back((YuEnergyRaw[i]*Yup0[YuChannel[i]] + Yup1[YuChannel[i]])/1000);
				std::pair<int, int> sectorRing = YuSectorRingMap[YuChannel[i]];
				YuSector.push_back(sectorRing.first);
				YuRing.push_back(sectorRing.second);
			}
			
			/*
			if(maxYu>0){
				YuEnergyRaw.push_back (maxYu );
				YuChannel.push_back ( channelYu );
				
				double tempYuE = maxYu *Yup0[channelYu]+Yup1[channelYu];
				YuEnergy.push_back ( tempYuE/1000. ); //convert from keV to MeV in the calibrated energy of the Yu detector
				//  YuEnergy.push_back ( tempYuE); //if you are trying to check the calibration do not divide by 1000
				// cout << "Yu channel " << channelYu << " energy raw " << maxYu << " Calibrated " << tempYuE << endl;
				
				if ( channelYu<16 ){
					YuSector.push_back ( 0 ); YuRing.push_back ( channelYu );
				} //sector 0
				else if ( channelYu>=16 && channelYu<32 ){
					YuSector.push_back ( 1 ); YuRing.push_back ( channelYu-16 );
				} //sector 1
				else if ( channelYu>=32 && channelYu<48 ){
					YuSector.push_back ( 2 ); YuRing.push_back ( channelYu-32 );
				} //sector 2
				else if ( channelYu>=48 && channelYu<64 ){
					YuSector.push_back ( 3 ); YuRing.push_back ( channelYu-48 );
				} //sector 3
				else if ( channelYu>=64 && channelYu<80 ){
					YuSector.push_back ( 4 ); YuRing.push_back ( channelYu-64 );
				} //sector 4
				else if ( channelYu>=80 && channelYu<96 ){
					YuSector.push_back ( 5 ); YuRing.push_back ( channelYu-80 );
				} //sector 5
				else if ( channelYu>=96 && channelYu<112){
					YuSector.push_back ( 6 ); YuRing.push_back ( channelYu-96 );
				} //sector 6
				else if ( channelYu>=112 && channelYu<128){
					YuSector.push_back ( 7 );YuRing.push_back ( channelYu-112 );
				} //sector 7
			}
			*/
		}
		
			
		


		tr_out->Fill(); //fill out the tree

		
		YuChannel.clear();
		YuEnergyRaw.clear();
		YuEnergy.clear();
		YuRing.clear();
		YuSector.clear();
		
		

		ev_num++; //looping over the while loop
	} //end of the main while loop

	f_out->cd(); //access the output file
	tr_out->Write(); //write the output file
	f_out->Close(); //close the output file

	return 0;


} //end of the script
