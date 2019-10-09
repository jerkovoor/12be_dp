//M. Vostinar
//This program is used to actually produce the histograms for analysis of the data
//There are a lot of things in it, because I was testing different things and we needed to check for them 
//One needs to turn on/off the things one wants to look at

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

#include <map>

#include "EnergyLoss.h"

double analysis(){
	//open the output file
	TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/Be/CarbCalibAlpha/Be_nopedestal_1_2Peaks_CCalib_test.root","RECREATE"); //change the path and name accordingly
  
	//Open the input files
	TChain *chain = new TChain ( "AutoTree" ); //chain the desired input files
	
	EnergyLoss* protonEL = new EnergyLoss("Proton_DeuteriumTarget.dat");
	
	std::map<int, double> YuAngleMap;
	YuAngleMap[0] = 31.68434901;
    
	//no target C data
	//chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Carbon/Decode4815Ped_NoTar_4992.root"); 
    
	//no target Be data
	/*chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5142.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5143.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5144.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5145.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5146.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5147.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5148.root");*/
	//chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5149.root");
	//chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5150.root");
	
	/*chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5151.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5152.root");
	//chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5153.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5154.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5155.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5156.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5157.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5158.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5159.root");
	
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5161.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5162.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5163.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5164.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5165.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5166.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5167.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5168.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5169.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5170.root");
	
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5171.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5172.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5173.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5174.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5176.root");
	chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Be/Decode4815Ped_5177.root");*/
	
	//Be data with target
	
	for(int run_num=5021;run_num<5113;run_num++){
		if(run_num==5026||run_num==5040||run_num==5043||run_num==5046||run_num==5047||run_num==5059||run_num==5062||run_num==5063||run_num==5093||run_num==5099||run_num==5101){
			continue;
		}else{
			string f_name=Form("/home/jerome/12Be_exp/Analysis/Be/CarbCalibAlpha/decodeBe_noYupedestal_1_2Peaks_CCalib_%i.root",run_num);
			chain->Add(f_name.c_str());
		}
	}
	
	for(int run_num=5115;run_num<5133;run_num++){
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
	}
	
	//Carbon data with target
	/*chain->Add ( "/home/jerome/12Be_exp/Analysis/C_calib/decode_noYupedestal_1_2Peaks_CCalib_5001.root" ); //change the path to the files and the file name accordingly
	chain->Add ( "/home/jerome/12Be_exp/Analysis/C_calib/decode_noYupedestal_1_2Peaks_CCalib_5002.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/C_calib/decode_noYupedestal_1_2Peaks_CCalib_5003.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/C_calib/decode_noYupedestal_1_2Peaks_CCalib_5004.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/C_calib/decode_noYupedestal_1_2Peaks_CCalib_5006.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/C_calib/decode_noYupedestal_1_2Peaks_CCalib_5007.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/C_calib/decode_noYupedestal_1_2Peaks_CCalib_5009.root" );*/
	
	
	//define the input variables
	//YYD detector
	int YdMul;
	vector<int>* YdChannel = new vector<int>();
	vector<double>* YdEnergyRaw= new vector<double>();
	vector<double>* YdEnergy= new vector<double>();
	vector<int>* YdRing= new vector<int>();
	vector<int>* YdSector= new vector<int>();
	
	//YYU detector
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
    
	//load the geometrical cuts used to look at different detectors; for example PID cut to look at Yu detector. 
	//change the path and name accordingly
	/* TFile *f_carbon = TFile::Open ( "/home/jerome/12Be_exp/scripts/scripts/CarbonNew.root");
	TCutG *cutC = (TCutG*) f_carbon->Get("C");*/
   
	/*TFile *f_beNoTC = TFile::Open ( "/home/jerome/12Be_exp/scripts/scripts/Be_noTarCut.root");
	TCutG *cutBeNT = (TCutG*) f_beNoTC->Get("Be");*/
   
   
	/*TFile *f_cut = TFile::Open("/home/jerome/12Be_exp/scripts/13Ccut.root"); //PID cut for C
	TCutG *pidcut = (TCutG*) f_cut->Get("CUTG"); */ // for C
	
	TFile *f_cut = TFile::Open("/home/jerome/12Be_exp/Analysis/Be/BeCut1.root"); //PID cut for Be
	TCutG *pidcut = (TCutG*) f_cut->Get("BeCut_CalSd1rSd2r");  // for Be
      
    
	//definition of histograms
	//the Sd1 rings face the beam and the Sd2 sectors face the beam
	//the plots seem better when one looks at ring vs ring but TRIUMF insists to use ring vs sector, so we have all of the figures
	TH1D *hSd1rMul = new TH1D("hSd1rMul","hSd1rMul",10,0,10);
	TH1D *hSd1sMul = new TH1D("hSd1sMul","hSd1sMul",10,0,10);
	TH1D *hSd2rMul = new TH1D("hSd2rMul","hSd2rMul",10,0,10);
	TH1D *hSd2sMul = new TH1D("hSd2sMul","hSd2sMul",10,0,10);
	TH1D *hSdr1 = new TH1D("hSdr1","Sdr1 calibrated singles",400,0,40); //Sdr1 calibrated and plotted to look at the peak for additional calibration
	TH1D *hSdr2 = new TH1D("hSdr2","Sdr2 calibrated singles",1000,0,200); //Sdr2 calibrated and plotted to look at the peak for additional calibration
	TH1D *hSds1 = new TH1D("hSds1","Sds1 calibrated singles",400,0,40); //Sds1 calibrated and plotted to look at the peak for additional calibration
	TH1D *hSds2 = new TH1D("hSds2","Sds2 calibrated singles",1000,0,200); //Sds2 calibrated and plotted to look at the peak for additional calibration
   
	//ploting channel vs channel for the two detectors to check where to put the requirements so to have less background in them     
	TH2D *hSdr1r2Ch = new TH2D("hSdr1r2Ch","Ring 1 vs ring 2",24,0,24,24,0,24); 
	TH2D *hSds1s2Ch = new TH2D("hSds1s2Ch","Sector 1 vs sector 2",32,0,32,32,0,32);
    
	//PID plotts
	TH2D *hSd1rSd2r = new TH2D("hSd1rSd2r","Sd1r vs Sd2r",2048,0,8192,2048,0,8192); //non calibrated Sd1r vs Sd2r
	TH2D *hSd1rSd2rCut = new TH2D("hSd1rSd2rCut","Sd1r vs Sd2r with correlation requirenment & IC gate",2048,0,8192,2048,0,8192); //non calibrated Sd1r vs Sd2r with a ring correlation requirenment and an IC gate
	TH2D *hSd1rSd2s = new TH2D("hSd1rSd2s","Sd1r vs Sd2s",2048,0,8192,2048,0,8192); //non calibrated Sd1r vs Sd2s
       
	TH2D *hCSd1rSd2r = new TH2D("hCSd1rSd2r","Calibrated Sd1r vs Sd2r",1500,0,150,600,0,60); // calibrated Sd1r vs Sd2r
	TH2D *hCSd1rSd2s = new TH2D("hCSd1rSd2s","Calibrated Sd1r vs Sd2s",1500,0,150,600,0,60); // calibrated Sd1r vs Sd2s
	TH2D *hCSd1sSd2r = new TH2D("hCSd1sSd2r","Calibrated Sd1s vs Sd2r",1500,0,150,600,0,60); // calibrated Sd1s vs Sd2r
	TH2D *hCSd1sSd2s = new TH2D("hCSd1sSd2s","Calibrated Sd1s vs Sd2s",1500,0,150,600,0,60); // calibrated Sd1s vs Sd2s
    
	TH2D *hCSd1rSd2rIC = new TH2D("hCSd1rSd2rIC","Cal Sd1r vs Sd2r, gated by IC",1500,0,150,600,0,60); // calibrated Sd1r vs Sd2r and gated by IC
	TH2D *hCSd1rSd2rICCut = new TH2D("hCSd1rSd2rICCut","Cal Sd1r vs Sd2r, gated by IC gated by the channel corr",1500,0,150,600,0,60); // calibrated Sd1r vs Sd2r gated by IC and required to be correlated with each other
	TH2D *hCSd1rSd2sIC = new TH2D("hCSd1rSd2sIC","Cal Sd1r vs Sd2s, gated by IC",1500,0,150,600,0,60); // calibrated Sd1r vs Sd2s and gated by IC
    
	TH2D *hCSd1rSd2rYdIC = new TH2D("hCSd1rSd2rYdIC","Cal Sd1r vs Sd2r, gated by Yd & IC",1500,0,150,600,0,60); // calibrated Sd1r vs Sd2r and gated by Yd & IC
	TH2D *hCSd1rSd2rYdICCut = new TH2D("hCSd1rSd2rYdICCut","Cal Sd1r vs Sd2r, gated by Yd & IC & chann corr",1500,0,150,600,0,60); // calibrated Sd1r vs Sd2r gated by Yd, IC and required to be correlated with each other
	TH2D *hCSd1rSd2sYdIC = new TH2D("hCSd1rSd2sYdIC","Cal Sd1r vs Sd2s, gated by Yd & IC",1500,0,150,600,0,60); // calibrated Sd1r vs Sd2s and gated by Yd and IC
            
	//IC non calibrated
	TH1D *hIC = new TH1D("hIC","IC raw energy",2048,0,8192);
    
	//calculating the variable bin size for the Sd1 detector to plot the energy vs angle
	double Sd1bins[25]={0}, Sdbintemp[25]={0};
	for(int i=0; i<25;i++){
		Sdbintemp[i]=TMath::RadToDeg()*TMath::ATan((Sdr0+(24-i)*((Sdr1-Sdr0)/24))/Sd1z);
	}
    
	for(int k=24; k>=0; k--){
		int a = TMath::Abs(k-24);
		Sd1bins[a]=Sdbintemp[k];
	}
	
	TH2D *hSd1An = new TH2D("hSd1An","S1dE vs Angle", 24,Sd1bins,400,0,40); //energy vs angle for Sd1
    
	//calculate the variable bins for the Yd detector to be able to plot energy vs angle
	double Ydbins[17]={0}, Ydbintemp[17]={0};
	for(int i=0; i<17;i++){
		Ydbintemp[i]=TMath::RadToDeg()*TMath::ATan((Ydr0+(16-i)*((Ydr1-Ydr0)/16))/Ydz);
	}
    
	for(int k=16; k>=0; k--){
		int a = TMath::Abs(k-16);
		Ydbins[a]=Ydbintemp[k];
	}
    
	TH2D *hYdAn = new TH2D("hYdAn","YdE vs Angle", 16,Ydbins,10000,0,5); //Yd energy vs angle no gates
	TH2D *hYdAnIC = new TH2D("hYdAnIC","YdE vs Angle with an IC gate", 16,Ydbins,10000,0,5); //Yd energy vs angle with gate on the IC peak
	TH2D *hYdAnPID = new TH2D("hYdAnPID","YdE vs Angle with an PID gate", 16,Ydbins,10000,0,5); //Yd energy vs angle with gate on the PID
    
	TH2D *hYdCsIch = new TH2D("hYdCsIch","Ring Yd vs CsI1",16,0,16,16,0,16); //Energy in Yd vs CsI1
    
	TH2D *hYdCsI[16]; //Energy in Yd vs CsI gain1
	TH2D *hYdCsI2[16];//Energy in Yd vs CsI gain2
    
	for(int i=0; i<16; i++){
		string nameYdCsI = Form("YdCsI_%i",i);
		hYdCsI[i] = new TH2D(nameYdCsI.c_str(),"Yd vs CsI1",2048,0,4096,10000,0,5);
     
		string nameYdCsI2 = Form("YdCsI2_%i",i);
		hYdCsI2[i] = new TH2D(nameYdCsI2.c_str(),"Yd vs CsI2",2048,0,4096,10000,0,5);
	}
        
	//calculate the variable bins for the Yu detector to plot the angles
	double Yubins[17]={0};
	
	for (int i=1; i<18; i++){
		Yubins[i-1]=180-TMath::RadToDeg()*TMath::ATan((50+(16-i)*4.94)/(0-Yuz));
	}
  	TH1D *hYuMul = new TH1D("hYuMul","hYuMul",10,0,10);
	TH1D *hYuEn = new TH1D("hYuEn","hYuEn",2500,0,10); //energy singles in Yu
	TH1D *hYuEnPID = new TH1D("hYuEnPID","hYuEnPID",2400,0,2.4);
	TH1D *hYuEn1 = new TH1D("hYuEn1","hYuEn1",2500,0,10);
	TH1D *hYuEnIC = new TH1D("hYuEnIC","hYuEnIC",2400,0,2.4);
	TH1D *hYuEnIC1 = new TH1D("hYuEnIC1","hYuEnIC1",2400,0,2.4);
	TH2D *hYuAn = new TH2D("hYuAn","YuE vs Angle",16,Yubins,250,0,10); //Energy vs angle in Yu no gates
	TH2D *hYuAn1 = new TH2D("hYuAn1","YuE vs Angle",16,Yubins,250,0,10);
	TH2D *hYuAnIC = new TH2D("hYuAnIC","YuE vs Angle with a gate on the IC",16,Yubins,240,0,2.4); //Energy vs angle in Yu with a gate on IC
	TH2D *hYuAnIC1 = new TH2D("hYuAnIC1","YuE vs Angle with a gate on the IC",16,Yubins,240,0,2.4);
	TH2D *hYuAnPID = new TH2D("hYuAnPID","YuE vs Angle with an IC and PID gate",16,Yubins,240,0,2.4); //Energy vs angle in Yu with a gate on IC and PID
	TH2D *hYuAnPID1 = new TH2D("hYuAnPID1","YuE vs Angle with an IC and PID gate",16,Yubins,240,0,2.4);
    
	//Q values
	TH1D *hQval = new TH1D("hQval","Q values",100,-10,10);
	TH1D *hQval1 = new TH1D("hQval1","Q values",100,-10,10);
    
	//start reading the tree
	int ev = chain->GetEntries(); //get the total number of entries
	cout << "Total number of events =" << ev << endl;
	int ev_num=0;
	
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
	
	        
	while(ev_num<=ev){
		if ( ev_num%10000==0 ) cout << "Current event = " << ev_num << "\r"<< flush;
		chain->GetEntry ( ev_num ); //get the current entry
		
		//fill the Sd singles histogrames
		if(Sd1rEnergy->size()>0 && Sd1rEnergy->at(0)>0){
			hSd1rMul->Fill(Sd1rMul);
			hSdr1->Fill(Sd1rEnergy->at(0));
		}
		
		if(Sd2rEnergy->size()>0 && Sd2rEnergy->at(0)>0){
			hSd2rMul->Fill(Sd2rMul);
			hSdr2->Fill(Sd2rEnergy->at(0));
		}
		
		if(Sd1sEnergy->size()>0 && Sd1sEnergy->at(0)>0){
			hSd1sMul->Fill(Sd1sMul);
			hSds1->Fill(Sd1sEnergy->at(0));
		}
		
		if(Sd2sEnergy->size()>0 && Sd2sEnergy->at(0)>0){
			hSd2sMul->Fill(Sd2sMul);
			hSds2->Fill(Sd2sEnergy->at(0));
		}
      
		//fill in the PID plots
		//Sd1s vs Sd2s
		if(Sd1sChannel->size()>0 && Sd2sChannel->size()>0 && Sd1sMul==1 && Sd2sMul==1){
			hSds1s2Ch->Fill ( Sd2sChannel->at ( 0 ),Sd1sChannel->at ( 0 ) ); //Sector vs Sector for Sd1 and Sd2 
		} 
		
		//Sd1r vs Sd2r
		if(Sd1rChannel->size()>0 && Sd2rChannel->size()>0 && Sd1rMul==1 && Sd2rMul==1){
			hSdr1r2Ch->Fill(Sd2rChannel->at(0),Sd1rChannel->at(0)); //Ring vs Ring for Sd1 and Sd2 
			hSd1rSd2r->Fill(Sd2rEnergyRaw->at(0),Sd1rEnergyRaw->at(0)); //Non-calibrated Sd1r vs Sd2r energy - PID plot
			if(TMath::Abs(Sd2rChannel->at(0)-Sd1rChannel->at(0))==1 && ICEnergyRaw>620 && ICEnergyRaw<1100){hSd1rSd2rCut->Fill(Sd2rEnergyRaw->at(0),Sd1rEnergyRaw->at(0));} //Non-calibrated Sd1r vs Sd2r energy with an ring correlation requirenment and IC gate - PID plot
		}
		
		//Sd1r vs Sd2s
		if(Sd1rChannel->size()>0 && Sd2sChannel->size()>0 && Sd1rMul==1 && Sd2sMul==1){
			hSd1rSd2s->Fill(Sd2sEnergyRaw->at(0),Sd1rEnergyRaw->at(0)); //Non-calibrated Sd1r vs Sd2s energy - PID plot
		} 
		
		//Sd1r vs Sd2r calibrated
		if(Sd1rEnergy->size()>0 && Sd2rEnergy->size()>0 && Sd1rMul==1 && Sd2rMul==1){
			hCSd1rSd2r->Fill(Sd2rEnergy->at(0),Sd1rEnergy->at(0)); //Calibrated Sd1r vs Sd2r energy with no gates
			if ( ICEnergyRaw>620 && ICEnergyRaw<1100){
				hCSd1rSd2rIC->Fill ( Sd2rEnergy->at ( 0 ),Sd1rEnergy->at ( 0 ) );  //Calibrated Sd1r vs Sd2r energy with IC gate
				if (TMath::Abs(Sd2rChannel->at(0)-Sd1rChannel->at(0))==1){
					hCSd1rSd2rICCut->Fill ( Sd2rEnergy->at ( 0 ),Sd1rEnergy->at ( 0 ) );
				} //Calibrated Sd1r vs Sd2r energy with IC gate and ring vs ring correlation
				if ( YdEnergy->size() >0 && YdMul==1 && YdEnergy->at( 0 ) >0){
					hCSd1rSd2rYdIC->Fill(Sd2rEnergy->at(0),Sd1rEnergy->at(0));//Calibrated Sd1r vs Sd2r energy with IC gate and Yd gate
					if (TMath::Abs(Sd2rChannel->at(0)-Sd1rChannel->at(0))==1){
						hCSd1rSd2rYdICCut->Fill(Sd2rEnergy->at(0),Sd1rEnergy->at(0));
					} //Calibrated Sd1r vs Sd2r energy with IC and Yd gate and ring vs ring correlation
				} //end of Yd gates
			} //end of IC gates
		} //end of Sd1r vs Sd2r calibrated

		//Sd1r vs Sd2s calibrated
		if(Sd1rEnergy->size()>0 && Sd2sEnergy->size()>0 && Sd1rMul==1 && Sd2sMul==1){
			hCSd1rSd2s->Fill(Sd2sEnergy->at(0),Sd1rEnergy->at(0));//Calibrated Sd1r vs Sd2s energy with no gates
			if(ICEnergyRaw>620 && ICEnergyRaw<1100 ){
				hCSd1rSd2sIC->Fill(Sd2sEnergy->at(0),Sd1rEnergy->at(0));//Calibrated Sd1r vs Sd2r energy with IC gate
				if ( YdEnergy->size() >0 && YdMul==1 && YdEnergy->at( 0 ) >0){hCSd1rSd2sYdIC->Fill(Sd2sEnergy->at(0),Sd1rEnergy->at(0));
				}//Calibrated Sd1r vs Sd2r energy with IC and Yd gate 
			} //end of IC gate
		}//end of Sd1r vs Sd2s calibrated

		//Sd1s vs Sd2r calibrated
		if(Sd1sEnergy->size()>0 && Sd2rEnergy->size()>0 && Sd1sMul==1 && Sd2rMul==1 && ICEnergyRaw>620 && ICEnergyRaw<1100){hCSd1sSd2r->Fill(Sd2rEnergy->at(0),Sd1sEnergy->at(0));}//Calibrated Sd1s vs Sd2r energy with IC gate
		if(Sd1sEnergy->size()>0 && Sd2sEnergy->size()>0 && Sd1sMul==1 && Sd2sMul==1 && ICEnergyRaw>620 && ICEnergyRaw<1100){hCSd1sSd2s->Fill(Sd2sEnergy->at(0),Sd1sEnergy->at(0));}//Calibrated Sd1s vs Sd2s energy with IC gate
      
		//fill in the IC
		if(ICEnergyRaw>0){
			hIC->Fill(ICEnergyRaw);
		}
		
		//fill in the Yd detector
		if(YdChannel->size()>0 && YdMul==1 && YdEnergy->size()>0){
			if ( YdEnergy->at ( 0 ) >0 && ( YdChannel->at ( 0 ) !=20 && YdChannel->at ( 0 ) !=0 && YdChannel->at ( 0 ) !=55 &&  YdChannel->at ( 0 ) !=127 ) ){ //channels are excluded because they are not calibrated and don work
				//calculate the angles for Yd detector
				ydM = Ydr0+ ( YdRing->at ( 0 ) *YChWidth )-YChMiddle;
				TVector3 YdposM ( 0,ydM,Ydz );
				YdposM.SetPhi ( TMath::Pi() /2-YdSector->at ( 0 ) *TMath::Pi() /4 ); //Pi/2 because the center of the sector is at 90 degrees, Pi/4 is because there is 8 sectors so it is 2Pi/8
				YdthetaM = beam.Angle ( YdposM ) *TMath::RadToDeg();
				
				hYdAn->Fill ( YdthetaM,YdEnergy->at ( 0 ) ); //all, no cuts
				if(ICEnergyRaw>620 && ICEnergyRaw<1100){hYdAnIC->Fill(YdthetaM,YdEnergy->at ( 0 ));} //IC gate for C 1500 < E < 2200; for Be 620< E < 1100
				if(Sd1rEnergy->size()>0 && Sd2rEnergy->size()>0 && Sd1rMul==1 && Sd2rMul==1){if(pidcut->IsInside(Sd2rEnergy->at(0),Sd1rEnergy->at(0))){hYdAnPID->Fill ( YdthetaM,YdEnergy->at ( 0 ) );}} //put the proper PID cut in here
				
				if(CsI1Channel->size()>0)if(CsI1EnergyRaw->at(0)>0 && CsI1Mul==1){{
				hYdCsIch->Fill(YdChannel->at(0),CsI1Channel->at(0));
				hYdCsI[CsI1Channel->at(0)]->Fill(CsI1EnergyRaw->at(0),YdEnergy->at(0)); }}
				
				if(CsI2Channel->size()>0)if(CsI2EnergyRaw->at(0)>0 && CsI2Mul==1){{
				hYdCsI2[CsI2Channel->at(0)]->Fill(CsI2EnergyRaw->at(0),YdEnergy->at(0)); }}
			}
		}
      
		//fill in the Yu detector
		if(YuChannel->size()>0  && YuEnergy->size()>0 ){
			//if (YuEnergy->at(0)>0.20 && ( YuChannel->at(0)!=82 && YuChannel->at(0)!=96 && YuChannel->at(0)!=106 && YuChannel->at(0)!=111 ))
			if (YuEnergy->at(0)>0.1 && ( YuChannel->at(0)!=82 && YuChannel->at(0)!=96 && YuChannel->at(0)!=106 && YuChannel->at(0)!=111)){
				//calculate the angles for the Yu detector
				yuM = Ydr0+ ( YuRing->at(0) *YChWidth )-YChMiddle;
				TVector3 YuposM ( 0,yuM,Yuz); //shifting the detecor
				YuposM.SetPhi ( TMath::Pi() /2-YuSector->at(0) *TMath::Pi() /4 ); //Pi/2 because the center of the sector is at 90degrees, Pi/4 is because there is 8 sectors so it is 2Pi/8
				YuthetaM = beam.Angle ( YuposM ) *TMath::RadToDeg();
				YuAngle->push_back(YuthetaM);
							
				//Adding the energy lost by the protons through the target and the deadlayers of the Yu detector (11/15/2018)  
				for(int RingNumber=0; RingNumber<16; RingNumber++){	    
					for(int j=0; j<8; j++){ 
						
						if(YuChannel->at(0)==j*16+RingNumber){
							
							loss1=(-b[RingNumber]-TMath::Sqrt(pow(b[RingNumber],2)-4*(a[RingNumber]-YuEnergy->at(0))*c[RingNumber]))/(2*(a[RingNumber]-YuEnergy->at(0)));// Boron dead layer
							loss2=(-e[RingNumber]-TMath::Sqrt((e[RingNumber]*e[RingNumber])-4*(d[RingNumber]-(YuEnergy->at(0)+loss1))*f[RingNumber]))/(2*(d[RingNumber]-(YuEnergy->at(0)+loss1)));// Aluminium dead layer
							loss3=(-h[RingNumber]-TMath::Sqrt((h[RingNumber]*h[RingNumber])-4*(g[RingNumber]-(YuEnergy->at(0)+loss1+loss2))*k[RingNumber]))/(2*(g[RingNumber]-(YuEnergy->at(0)+loss1+loss2)));// Target
							if(RingNumber == 0) {
								double newEnergy = protonEL->AddBack(YuEnergy->at(0) + loss1 + loss2, 30./1000./cos(YuAngleMap[RingNumber]*M_PI/180.));
								std::cout << YuEnergy->at(0) + loss1 + loss2 << '\t' << YuEnergy->at(0) + loss1 + loss2 + loss3 << '\t' << newEnergy << std::endl;
							}
							loss=loss1+loss2+loss3; //
							//cout << ev_num << "	" << RingNumber << "		" << YuEnergy->at ( 0 ) << loss1 << loss2 << loss3 << loss << endl;
							hYuEn->Fill(YuEnergy->at(0)+loss);
							hYuAn->Fill(YuthetaM,YuEnergy->at(0)+loss);
							if(ICEnergyRaw>620 && ICEnergyRaw<1100){
								hYuAnIC->Fill(YuthetaM,YuEnergy->at (0)+loss);
								hYuEnIC->Fill(YuEnergy->at (0)+loss);
							} //IC gate for C 1500 < E < 2200; for Be 620< E < 1100
							
							if(Sd2rEnergy->size()>0 && Sd1rEnergy->size()>0){
								if(pidcut->IsInside(Sd2rEnergy->at(0),Sd1rEnergy->at(0)) && ICEnergyRaw>620 && ICEnergyRaw<1100 ){
									hYuEnPID->Fill(YuEnergy->at(0)+loss);
									hYuAnPID->Fill(YuthetaM,YuEnergy->at(0)+loss);
									Qval = ( 1+mejec/mrecoil ) * ( YuEnergy->at(0)+loss) - ( 1 - mbeam/mrecoil ) * ( kBeam ) - 2 * TMath::Sqrt ( mbeam*mejec* ( YuEnergy->at(0)+loss) * ( kBeam ) ) / mrecoil * TMath::Cos ( YuthetaM * TMath::Pi() / 180. );
									hQval->Fill ( Qval );
								}	
							}
						}
					}//End of j loop
				}//End of Ringnumber loop
				
				hYuAn1->Fill(YuthetaM,YuEnergy->at(0)); //all, no cuts
				if(ICEnergyRaw>620 && ICEnergyRaw<1100){
					hYuAnIC1->Fill(YuthetaM,YuEnergy->at (0));
					hYuEnIC1->Fill(YuEnergy->at (0));
				} //IC gate for C 1500 < E < 2200; for Be 620< E < 1100
				hYuEn1->Fill(YuEnergy->at(0));
				hYuMul->Fill(YuMul);
								
				if(Sd2rEnergy->size()>0 && Sd1rEnergy->size()>0){
					if (pidcut->IsInside(Sd2rEnergy->at ( 0 ),Sd1rEnergy->at ( 0 ) ) && ICEnergyRaw>620 && ICEnergyRaw<1100 ){
						hYuAnPID1->Fill ( YuthetaM,YuEnergy->at(0) );
						Qval = ( 1+mejec/mrecoil ) * ( YuEnergy->at(0) ) - ( 1 - mbeam/mrecoil ) * ( kBeam ) - 2 * TMath::Sqrt ( mbeam*mejec* ( YuEnergy->at(0) ) * ( kBeam ) ) / mrecoil * TMath::Cos ( YuthetaM * TMath::Pi() / 180. );
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
	hSdr1->Write();
	hSdr2->Write();
	hSds1->Write();
	hSds2->Write();
	hSd1rMul->Write();
	hSd1sMul->Write();
	hSd2rMul->Write();
	hSd2sMul->Write();
	
	//channel ring vs ring and sector vs sector for the Sd detectors 
	hSdr1r2Ch->Write();
	hSds1s2Ch->Write();
	
	hSd1rSd2r->Write(); //non calibrated ring vs ring
	hSd1rSd2rCut->Write(); //non calibrated ring vs ring but taking into account the diagonal
	hCSd1rSd2r->Write();
	hCSd1rSd2rIC->Write(); //calibrated ring vs ring with IC cut
	hCSd1rSd2rICCut->Write(); //calibrated ring vs ring with IC cut and diagonal taken into account
	hCSd1rSd2rYdIC->Write();
	hCSd1rSd2rYdICCut->Write();
	hCSd1rSd2sYdIC->Write();
	
	hSd1rSd2s->Write(); //non calibrated ring vs sector
	
	hCSd1rSd2s->Write(); //calibrated ring vs sector
	hCSd1rSd2sIC->Write(); //calibrated ring vs sector with IC cut
	
	hCSd1sSd2r->Write(); //calibrated sector vs ring
	hCSd1sSd2s->Write(); //calibrated sector vs sector
	
	hSd1An->Write();
	
	//IC
	hIC->Write(); 
	
	//Yd
	hYdAn->Write(); 
	hYdAnIC->Write(); 
	hYdAnPID->Write(); 
	hYdCsIch->Write();
	for(int i=0; i<16; i++){hYdCsI[i]->Write();}
	for(int i=0; i<16; i++){hYdCsI2[i]->Write();}
	
	//Yu
	hYuAn->Write();
	hYuAn1->Write();
	hYuAnIC->Write();
	hYuAnIC1->Write();
	hYuAnPID->Write();
	hYuAnPID1->Write();
	hYuEn->Write();
	hYuEnIC->Write();
	hYuEnIC1->Write();
	hYuMul->Write();
	hYuEnPID->Write();
	hYuEn1->Write();
	hQval->Write();
	hQval1->Write();
	
	f_out->Close();
	
	return 0;
	
} //end of the program
