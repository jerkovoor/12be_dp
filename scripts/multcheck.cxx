using namespace std;

#include "TFile.h" 
#include "TChain.h"
#include "TCutG.h"
#include "TVector3.h"
#include "TMath.h"
#include <cmath>
#include "iostream"
#include "fstream"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"

typedef struct YuDet{
	int mult;
	int channel;
	int ring;
	int sector;
	double energyRaw;
	double energy;
} YuDet;
	
struct sortByEnergy {
	inline bool operator() (const YuDet& En1,
							const YuDet& En2){
		return (En1.energy>En2.energy);
	}
};


double multcheck(){
	//open the output file
	TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/Be_newdecode/Be_Yupedestal_multcheck2.root","RECREATE"); //change the path and name accordingly
  
	//Open the input files
	TChain *chain = new TChain ( "AutoTree" ); //chain the desired input files
	
	for(int run_num=5021;run_num<5113;run_num++){
		if(run_num==5026||run_num==5040||run_num==5043||run_num==5046||run_num==5047||run_num==5059||run_num==5062||run_num==5063||run_num==5093||run_num==5099||run_num==5101){
			continue;
		}else{
			string f_name=Form("/home/jerome/12Be_exp/Analysis/Be_newdecode/decodeBe_Yupedestal_%i.root",run_num);
			chain->Add(f_name.c_str());
		}
	}
	
	/*for(int run_num=5115;run_num<5133;run_num++){
		string f_name=Form("/home/jerome/12Be_exp/Analysis/Be/CarbCalibAlpha/decodeBe_noYupedestal_1_2Peaks_CCalib_%i.root",run_num);
		chain->Add(f_name.c_str());
	}
	
	for(int run_num=5180;run_num<5223;run_num++){
		if(run_num==5187||run_num==5200){
			continue;
		}else{
			string f_name=Form("/home/jerome/12Be_exp/Analysis/Be/CarbCalibAlpha/decodeBe_noYupedestal_1_2Peaks_CCalib_%i.root",run_num);
			chain->Add(f_name.c_str());
		}
	}*/
	
	//define the input variables
	
	//YYD detector
	int YdMul;
	vector<int>* YdChannel = new vector<int>();
	vector<double>* YdEnergyRaw= new vector<double>();
	vector<double>* YdEnergy= new vector<double>();
	vector<int>* YdRing= new vector<int>();
	vector<int>* YdSector= new vector<int>();
	
	//YYU detector
	/*int YuMul;
	vector<int> YuChannel;
	vector<double> YuEnergyRaw;
	vector<double> YuEnergy;
	vector<int> YuRing;
	vector<int> YuSector;*/
	
	int YuMul;
	vector<int>* YuChannel= new vector<int>();
	vector<double>* YuEnergyRaw= new vector<double>();
	vector<double>* YuEnergy= new vector<double>();
	vector<int>* YuRing= new vector<int>();
	vector<int>* YuSector= new vector<int>();
	
	//Purpose of two readouts is to increase the range of operation (multiply it by 2 in total)
	//CsI 1 detector
	int CsI1Mul;
	vector<int>* CsI1Channel= new vector<int>(); //Channel corresponds to one christal - 1 means that the gain is set to one value
	vector<double>* CsI1EnergyRaw= new vector<double>();
	
	//CsI 2 detector
	int CsI2Mul;
	vector<int>* CsI2Channel= new vector<int>();//Channel corresponds to one christal - 2 means that the gain is set to a different value than 1
	vector<double>* CsI2EnergyRaw= new vector<double>();
	
	//IC chamber
	int ICChannel;
	double ICEnergyRaw;
	
	//Sdd1 detector
	int Sd1rMul;
	vector<int>* Sd1rChannel= new vector<int>(); //Ring Channels
	vector<double>* Sd1rEnergyRaw= new vector<double>();
	vector<double>* Sd1rEnergy= new vector<double>();
 
	//Sd1s detector
	int Sd1sMul;
	vector<int>* Sd1sChannel = new vector<int>(); //Sector Channels
	vector<double>* Sd1sEnergyRaw= new vector<double>();
	vector<double>* Sd1sEnergy= new vector<double>();
	
	//Sd2r detector
	int Sd2rMul;
	vector<int>* Sd2rChannel= new vector<int>(); //Ring Channels
	vector<double>* Sd2rEnergyRaw= new vector<double>();
	vector<double>* Sd2rEnergy= new vector<double>();
	
	//Sd2s detector
	int Sd2sMul;
	vector<int>* Sd2sChannel= new vector<int>(); //Sector Channels
	vector<double>* Sd2sEnergyRaw= new vector<double>();
	vector<double>* Sd2sEnergy= new vector<double>();
	
	//Sur
	int SurMul;
	vector<int>* SurChannel= new vector<int>(); //Ring Channels
	vector<double>* SurEnergyRaw= new vector<double>();
	
	//Sus
	int SusMul;
	vector<int>* SusChannel= new vector<int>(); //Sector Channels
	vector<double>* SusEnergyRaw= new vector<double>();
	
	//reading the input tree
	chain->SetBranchAddress ( "YdMul",&YdMul );
	chain->SetBranchAddress ( "YdChannel",&YdChannel );
	chain->SetBranchAddress ( "YdEnergyRaw",&YdEnergyRaw );
	chain->SetBranchAddress ( "YdEnergy",&YdEnergy );
	chain->SetBranchAddress ( "YdRing",&YdRing );
	chain->SetBranchAddress ( "YdSector",&YdSector );
	
	chain->SetBranchAddress ( "YuMul",&YuMul );
	chain->SetBranchAddress ( "YuChannel",&YuChannel );
	chain->SetBranchAddress ( "YuEnergyRaw",&YuEnergyRaw );
	chain->SetBranchAddress ( "YuEnergy",&YuEnergy );
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
	chain->SetBranchAddress ( "Sd1rEnergy",&Sd1rEnergy );
	
	chain->SetBranchAddress ( "Sd1sMul",&Sd1sMul );
	chain->SetBranchAddress ( "Sd1sChannel",&Sd1sChannel );
	chain->SetBranchAddress ( "Sd1sEnergyRaw",&Sd1sEnergyRaw );
	chain->SetBranchAddress ( "Sd1sEnergy",&Sd1sEnergy );
	
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
	float Ydr0,Ydr1,Ydz, Yuz, Sd1z, Sd2z, Sdr0, Sdr1; // Ydz distance from the target, Ydr0 inner radius, Ydz outer radius, same for the Sd

	while ( getline ( geometry,read_geometry ) ){ //in the calibration file named geometry start reading the lines
		
		if ( read_geometry.find ( "YD_DISTANCE",0 ) !=string::npos ){ //if you find the "YD_DISTANCE" before the end of the file
			iss.clear();
			iss.str ( read_geometry );
			iss >> name >> dummy >> Ydz;
			// cout << " name " << name << " / " << dummy << " number " << Ydz << endl;
        }

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

		if ( read_geometry.find ( "SD1_DISTANCE",0 ) !=string::npos ){
			iss.clear();
			iss.str ( read_geometry );
			iss >> name >> dummy >> Sd1z;
		}

		if ( read_geometry.find ( "SD2_DISTANCE",0 ) !=string::npos ){
			iss.clear();
			iss.str ( read_geometry );
			iss >> name >> dummy >> Sd2z;
		}

		if ( read_geometry.find ( "SD_INNER_RADIUS",0 ) !=string::npos ){
			iss.clear();
			iss.str ( read_geometry );
			iss >> name >> dummy >> Sdr0;
		}

		if ( read_geometry.find ( "SD_OUTER_RADIUS",0 ) !=string::npos ){
			iss.clear();
			iss.str ( read_geometry );
			iss >> name >> dummy >> Sdr1;
		}
	}
	//end of the geometry file
	
	TFile *f_cut = TFile::Open("/home/jerome/12Be_exp/Analysis/Be/BeCut1.root"); //PID cut for Be
	TCutG *pidcut = (TCutG*) f_cut->Get("BeCut_CalSd1rSd2r");  // for Be
	
	double Yubins[17]={0};
	
	for (int i=1; i<18; i++){
		Yubins[i-1]=180-TMath::RadToDeg()*TMath::ATan((50+(16-i)*4.94)/(0-Yuz));
	}
	TH1D *hYuMul1 = new TH1D("hYuMul1","hYuMul1",10,0,10);
	TH1D *hYuMul2 = new TH1D("hYuMul2","hYuMul2",10,0,10);
	TH1D *hYuEn1 = new TH1D("hYuEn1","hYuEn1",2400,0,2.4); //energy singles in Yu
	TH1D *hYuEn2 = new TH1D("hYuEn2","hYuEn2",2400,0,2.4);
	TH1D *hYuEnIC1 = new TH1D("hYuEnIC1","hYuEnIC1",2400,0,2.4);
	TH1D *hYuEnIC2 = new TH1D("hYuEnIC2","hYuEnIC2",2400,0,2.4);
	TH1D *hYuEnPID1 = new TH1D("hYuEnPID1","hYuEnPID1",2400,0,2.4);
	TH1D *hYuEnPID2 = new TH1D("hYuEnPID2","hYuEnPID2",2400,0,2.4);
	TH2D *hYuAn1 = new TH2D("hYuAn1","YuE vs Angle1",16,Yubins,2500,0,10); //Energy vs angle in Yu no gates
	TH2D *hYuAn2 = new TH2D("hYuAn2","YuE vs Angle2",16,Yubins,2500,0,10);
	TH2D *hYuAnIC1 = new TH2D("hYuAnIC1","YuE vs Angle with a gate on the IC1",16,Yubins,2400,0,2.4); //Energy vs angle in Yu with a gate on IC
	TH2D *hYuAnIC2 = new TH2D("hYuAnIC2","YuE vs Angle with a gate on the IC2",16,Yubins,2400,0,2.4);
	TH2D *hYuAnPID1 = new TH2D("hYuAnPID1","YuE vs Angle with an IC and PID gate1",16,Yubins,2400,0,2.4); //Energy vs angle in Yu with a gate on IC and PID
	TH2D *hYuAnPID2 = new TH2D("hYuAnPID2","YuE vs Angle with an IC and PID gate2",16,Yubins,2400,0,2.4);
	TH2D *hYuEnM = new TH2D("hYuEnM","YuE1 vs YuE2 for multiplicity 2",2400,0,2.4,2400,0,2.4);
	TH2D *hYuEnMICPID = new TH2D("hYuEnMICPID","YuE1 vs YuE2 for multiplicity 2 with IC and PID gates",2400,0,2.4,2400,0,2.4);
    
	//Q values
	TH1D *hQval1 = new TH1D("hQval1","Q values1",100,-10,10);
	TH1D *hQval2 = new TH1D("hQval2","Q values2",100,-10,10);
    
	//start reading the tree
	int ev = chain->GetEntries(); //get the total number of entries
	cout << "Total number of events =" << ev << endl;
	int ev_num=0;
	TH1D *hYuEnSize = new TH1D("hYuEnSize","hYuEnSize",10,0,10);
	//variable calculation for the YY1 detectors used in angle calculations
	float YChWidth = ( Ydr1 - Ydr0 ) /16.;
	float YChMiddle = YChWidth/2;
	float yuM=0, ydM=0;
	double YuthetaM, YdthetaM; //angle for Yu/Yd
	TVector3 beam(0,0,1); //beam vector
    
	//variable calculation for the S3 detectors used in angle calculations
	float SdChWidth = ( Sdr1 - Sdr0 ) /24.;
	float SdChMiddle = SdChWidth/2;
	float SdM=0;
	double SdthetaM; //angle for Sd1
	
	//variable definition for the Q value calculations
	float amu = 931.5; // atomic mass unit in MeV
	float massEjec = 938.28; //mass of the proton in MeV/c2 
	float kBeam = 114; //put the correct value; beam energy
	float mbeam = 12 * amu;  //mass of the beam (12Be or 12C) in MeV
	float mrecoil = 13 * amu;  //mass of the recoil (13Be or 13C) in MeV
	float mejec = 1 * amu; //mass of the proton
    

	double Qval; //Q value variable
	vector<double> *YuAngle = new vector<double>; //Yu angle variable definition
    
    
	double ringB[16], ringAl[16], ringD[16], a[16], b[16], c[16], d[16], e[16], f[16], g[16], h[16], k[16], loss, loss1, loss2, loss3;

	ifstream Bfit;
	Bfit.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/B_Dead_layer/B_fitparameters.txt");

	if(Bfit.is_open()){
		Bfit.ignore(256,'\n');
		
		for(int i=0;i<16;i++){
			Bfit >> ringB[i] >> a[i] >> b[i] >> c[i];
		}
	}

	ifstream Alfit;
	Alfit.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Al_Dead_layer/Al_fitparameters.txt");

	if(Alfit.is_open()){
		Alfit.ignore(256,'\n');
		
		for(int i=0;i<16;i++){
			Alfit >> ringAl[i] >> d[i] >> e[i] >> f[i];
		}
	}
    
    
	ifstream Dfit;
	Dfit.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Target/D2_fitparameters.txt");
	
	if(Dfit.is_open()){
		Dfit.ignore(256,'\n');
    
		for(int i=0;i<16;i++){
			Dfit >> ringD[i] >> g[i] >> h[i] >> k[i];
		}
	}
	
	////////////////////////////////
	// Event by event starts here //
	////////////////////////////////
	        
	while(ev_num<=ev){
		if ( ev_num%10000==0 ) cout << "Current event = " << ev_num << "\r"<< flush;
		chain->GetEntry ( ev_num ); //get the current entry
		
		vector<YuDet> YuDetector;
		for(size_t i = 0; i < YuEnergy->size(); i++) {
			if(YuChannel->at(i) == 82 || YuChannel->at(i) == 96 || YuChannel->at(i) == 106 || YuChannel->at(i) == 111) continue;
			YuDet hit = {YuMul, YuChannel->at(i), YuRing->at(i), YuSector->at(i), YuEnergyRaw->at(i), YuEnergy->at(i)};
			// YuDet hit;
			// hit.mult = YuMul;
			// hit.channel = YuChannel->at(i);
			//hit.energy = YuEnergy->at(i);
			YuDetector.push_back(hit);
		}
		
		
		//if(!YuEnergy->empty()) {
		//	//cout << YuEnergy->at(0) << '\t' << YuDetector[0].energy << endl;
		//	//cout << YuEnergy->at(0) << '\t' << YuDetector[0].energy << endl;
		//	cout << YuMul << '\t' << YuChannel->at(0) << '\t' << YuRing->at(0) << '\t' << YuSector->at(0) << '\t' << YuEnergyRaw->at(0) << '\t' << YuEnergy->at(0) << endl;
		//	cout << YuDetector[0].mult << '\t' << YuDetector[0].channel << '\t' << YuDetector[0].ring << '\t' << YuDetector[0].sector << '\t' << YuDetector[0].energyRaw <<
		//		'\t' << YuDetector[0].energy << endl;
		//}
		
		//fill in the Yu detector
		if(!YuDetector.empty()){
			std::sort(YuDetector.begin(),YuDetector.end(),
			sortByEnergy());
		}
		

		if(YuDetector.size() == 2) {
			if(YuDetector[0].energy>0){
				hYuEnM->Fill(YuDetector[0].energy,YuDetector[1].energy);
				if(Sd2rEnergy->size()>0 && Sd1rEnergy->size()>0){
					if (pidcut->IsInside(Sd2rEnergy->at ( 0 ),Sd1rEnergy->at ( 0 ) ) && ICEnergyRaw>620 && ICEnergyRaw<1100 ){
						hYuEnMICPID->Fill(YuDetector[0].energy,YuDetector[1].energy);
					}
				}
			}
		}
		
		//for(auto hit : YuDetector) {
		/*for(size_t i = 0; i < YuDetector.size(); i++) {
			if(YuDetector[i].mult != 2) continue;
			yuM = Ydr0 + (hit.ring*YChWidth) - YChMiddle;
			TVector3 YuposM(0, yuM, Yuz); // shifting the detector
			YuposM.SetPhi ( TMath::Pi() /2-YuDetector[i].sector *TMath::Pi() /4 ); //Pi/2 because the center of the sector is at 90degrees, Pi/4 is because there is 8 sectors so it is 2Pi/8
			YuthetaM = beam.Angle ( YuposM ) *TMath::RadToDeg();
			YuAngle->push_back(YuthetaM);
			hYuAn1->Fill(YuthetaM,YuDetector[0].energy); //all, no cuts
			hYuEnM->Fill(YuDetector[0].energy,YuDetector[1].energy);
			if(ICEnergyRaw>620 && ICEnergyRaw<1100){
				hYuEnIC1->Fill(YuEnergy->at (0));
				hYuAnIC1->Fill(YuthetaM,YuEnergy->at (0));
			} //IC gate for C 1500 < E < 2200; for Be 620< E < 1100
			hYuEn1->Fill(YuEnergy->at(0));
			hYuMul1->Fill(YuMul);
		}*/
		
		
		if(YuDetector.size()>0){ //Looking at channel number 28 (Ring no. 12)
			//if (YuEnergy->at(0)>0.20 && ( YuChannel->at(0)!=82 && YuChannel->at(0)!=96 && YuChannel->at(0)!=106 && YuChannel->at(0)!=111 ))
			if (YuDetector[0].energy>0.1){
				//calculate the angles for the Yu detector
				yuM = Ydr0+ ( YuDetector[0].ring*YChWidth )-YChMiddle;
				TVector3 YuposM ( 0,yuM,Yuz); //shifting the detecor
				YuposM.SetPhi ( TMath::Pi() /2-YuDetector[0].sector *TMath::Pi() /4 ); //Pi/2 because the center of the sector is at 90degrees, Pi/4 is because there is 8 sectors so it is 2Pi/8
				YuthetaM = beam.Angle ( YuposM ) *TMath::RadToDeg();
				YuAngle->push_back(YuthetaM);
				hYuAn1->Fill(YuthetaM,YuDetector[0].energy); //all, no cuts
				if(ICEnergyRaw>620 && ICEnergyRaw<1100){
					hYuEnIC1->Fill(YuDetector[0].energy);
					hYuAnIC1->Fill(YuthetaM,YuDetector[0].energy);
				} //IC gate for C 1500 < E < 2200; for Be 620< E < 1100
				hYuEn1->Fill(YuDetector[0].energy);
				hYuMul1->Fill(YuMul);
				
				if(Sd2rEnergy->size()>0 && Sd1rEnergy->size()>0){
					if (pidcut->IsInside(Sd2rEnergy->at ( 0 ),Sd1rEnergy->at ( 0 ) ) && ICEnergyRaw>620 && ICEnergyRaw<1100 ){
						hYuEnPID1->Fill(YuDetector[0].energy);
						hYuAnPID1->Fill ( YuthetaM,YuDetector[0].energy );
						Qval = ( 1+mejec/mrecoil ) * ( YuDetector[0].energy ) - ( 1 - mbeam/mrecoil ) * ( kBeam ) - 2 * TMath::Sqrt ( mbeam*mejec* ( YuDetector[0].energy ) * ( kBeam ) ) / mrecoil * TMath::Cos ( YuthetaM * TMath::Pi() / 180. );
						hQval1->Fill ( Qval );
					}
				}
			}
		}
		ev_num++;
	} //end of the main while loop
    
  
	YuAngle->clear();
    
	f_out->cd();
    
	//write the histograms in the output root file
	//singles in the Sd detectors
	
	//Yu
	hYuAn1->Write();
	hYuAn2->Write();
	hYuAnIC1->Write();
	hYuAnIC2->Write();
	hYuAnPID1->Write();
	hYuAnPID2->Write();
	hYuEn1->Write();
	hYuEn2->Write();
	hYuEnIC1->Write();
	hYuEnIC2->Write();
	hYuMul1->Write();
	hYuMul2->Write();
	hYuEnPID1->Write();
	hYuEnPID2->Write();
	hQval1->Write();
	hQval2->Write();
	hYuEnSize->Write();
	hYuEnM->Write();
	hYuEnMICPID->Write();
	hYuEnM->Draw("colz");
	hYuEnMICPID->Draw("colz");
	f_out->Close();
	
	return 0;
	
} //end of the program