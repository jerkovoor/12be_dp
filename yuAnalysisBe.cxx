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
#include <TRandom3.h>
#include "EnergyLoss.h"

typedef struct YuDet{
	int mult;
	int channel;
	int ring;
	int sector;
	double energyRaw;
	double energy;
} YuDet;
	
typedef struct S3Det{
	int mult;
	int channel;
	double energyRaw;
	double energy;
} Sd1rDet, Sd2rDet;

struct sortByEnergy {
	inline bool operator() (const YuDet& En1,
							const YuDet& En2){
		return (En1.energy>En2.energy);
	}
};

struct sortByEnergyS3 {
	inline bool operator() (const Sd1rDet& EnS3_1,
							const Sd1rDet& EnS3_2){
		return (EnS3_1.energy>EnS3_2.energy);
	}
};

double yuAnalysisBe(){
	//open the output file
	TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/Be/CarbonGain/Be_pedestal_TRIUMF_DL_CarbonGain_QvalMinGS_Random_NewDecode_Cut2_Yu.root","RECREATE"); //change the path and name accordingly
  
	//Open the input files
	TChain *chain = new TChain ( "AutoTree" ); //chain the desired input files
	
	//No target Be data
	
	/*chain->Add ("/home/jerome/12Be_exp/Analysis/Be_notarget/decodeBe_notarget_pedestal5018.root");
	chain->Add ("/home/jerome/12Be_exp/Analysis/Be_notarget/decodeBe_notarget_pedestal5135.root");
	
	for(int run_num=5140;run_num<5178;run_num++){
		if(run_num==5141||run_num==5149||run_num==5150||run_num==5153||run_num==5160||run_num==5175){
			continue;
		}else{
			string f_name=Form("/home/jerome/12Be_exp/Analysis/Be_notarget/decodeBe_notarget_pedestal%i.root",run_num);
			chain->Add(f_name.c_str());
		}
	}*/
	
		//Be data with target
	
	for(int run_num=5021;run_num<5113;run_num++){
		if(run_num==5026||run_num==5040||run_num==5043||run_num==5046||run_num==5047||run_num==5059||run_num==5062||run_num==5063||run_num==5093||run_num==5099||run_num==5101){
			continue;
		}else{
			//string f_name=Form("/home/jerome/12Be_exp/Analysis/Be/AlphaOnly/decode_Yupedestal_%i.root",run_num);
			//string f_name=Form("/home/jerome/12Be_exp/Analysis/TRIUMF_DL/decodeBe_Yupedestal_TDL_%i.root",run_num);
			string f_name=Form("/home/jerome/12Be_exp/Analysis/Be_newdecode/decodeBe_Yupedestal_%i.root",run_num);// new decode
			
			chain->Add(f_name.c_str());
		}
	}
	
	for(int run_num=5115;run_num<5133;run_num++){
		//string f_name=Form("/home/jerome/12Be_exp/Analysis/Be/AlphaOnly/decode_Yupedestal_%i.root",run_num);
		//string f_name=Form("/home/jerome/12Be_exp/Analysis/TRIUMF_DL/decodeBe_Yupedestal_TDL_%i.root",run_num);
		string f_name=Form("/home/jerome/12Be_exp/Analysis/Be_newdecode/decodeBe_Yupedestal_%i.root",run_num);// new decode
		chain->Add(f_name.c_str());
	}
	
	for(int run_num=5180;run_num<5223;run_num++){
		if(run_num==5187||run_num==5200){
			continue;
		}else{
			//string f_name=Form("/home/jerome/12Be_exp/Analysis/Be/AlphaOnly/decode_Yupedestal_%i.root",run_num);
			//string f_name=Form("/home/jerome/12Be_exp/Analysis/TRIUMF_DL/decodeBe_Yupedestal_TDL_%i.root",run_num);
			string f_name=Form("/home/jerome/12Be_exp/Analysis/Be_newdecode/decodeBe_Yupedestal_%i.root",run_num);// new decode
			chain->Add(f_name.c_str());
		}
	}
	
	
	
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
	
	
	/*TFile *f_cut = TFile::Open("/home/jerome/12Be_exp/scripts/13Ccut.root"); //PID cut for C
	TCutG *pidcut = (TCutG*) f_cut->Get("CUTG");*/ // for C
	
	/*TFile *f_cut = TFile::Open("/home/jerome/12Be_exp/Analysis/Be/BeCut1.root"); //PID cut for Be
	TCutG *pidcut = (TCutG*) f_cut->Get("BeCut_CalSd1rSd2r"); */ // for Be
	
	TFile *f_cut = TFile::Open("/home/jerome/12Be_exp/Analysis/Be/CarbonGain/BeCcutFull2.root"); //PID cut for Be
	TCutG *pidcut = (TCutG*) f_cut->Get("BeCut_CalSd1rSd2rFull2");
	
	//definition of histograms
	//calculate the variable bins for the Yu detector to plot the angles
	double Yubins[17]={0};
	
	for (int i=0; i<17; i++){
		Yubins[i]=180-TMath::RadToDeg()*TMath::ATan((50+(16-i)*4.94)/(0-Yuz));
		//cout << Yubins[i] << endl;
	}
	
	TH2D *hCSd1rSd2r = new TH2D("hCSd1rSd2r","Calibrated Sd1r vs Sd2r",1500,0,150,600,0,40); // calibrated Sd1r vs Sd2r
	TH2D *hCSd1rSd2rIC = new TH2D("hCSd1rSd2rIC","Cal Sd1r vs Sd2r, gated by IC",1500,0,150,600,0,40); // calibrated Sd1r vs Sd2r and gated by IC
	TH2D *hCSd1rSd2rICCut = new TH2D("hCSd1rSd2rICCut","Cal Sd1r vs Sd2r, gated by IC gated by the channel corr",1500,0,150,600,0,40); // calibrated Sd1r vs Sd2r gated by IC and required to be
	TH2D *hCSd1rSd2rYdIC = new TH2D("hCSd1rSd2rYdIC","Cal Sd1r vs Sd2r, gated by Yd & IC",1500,0,150,600,0,40); // calibrated Sd1r vs Sd2r and gated by Yd & IC
	TH2D *hCSd1rSd2rYdICCut = new TH2D("hCSd1rSd2rYdICCut","Cal Sd1r vs Sd2r, gated by Yd & IC & chann corr",1500,0,150,600,0,40); // calibrated Sd1r vs Sd2r gated by Yd, IC and required to be correlated with each other
	
  	TH1D *hYuMul = new TH1D("hYuMul","hYuMul",10,0,10);
	
	TH1D *hYuEn = new TH1D("hYuEn","hYuEn",2500,0,10); //energy singles in Yu
	TH1D *hYuEnT = new TH1D("hYuEnT","hYuEn with target energy loss",2500,0,10);
	TH1D *hYuEn1 = new TH1D("hYuEn1","hYuEn1",2500,0,10);
	
	TH1D *hYuEnPID = new TH1D("hYuEnPID","hYuEnPID",2400,0,2.4);//6000,0,6; 2400,0,2.4
	
	TH1D *hYuEnIC = new TH1D("hYuEnIC","hYuEnIC",2400,0,2.4);
	TH1D *hYuEnIC1 = new TH1D("hYuEnIC1","hYuEnIC1",2400,0,2.4);
	
	TH2D *hYuAn = new TH2D("hYuAn","YuE vs Angle",16,Yubins,250,0,10); //Energy vs angle in Yu no gates
	TH2D *hYuAnT = new TH2D("hYuAnT","YuE vs Angle with target energy loss",16,Yubins,250,0,10);
	TH2D *hYuAn1 = new TH2D("hYuAn1","YuE vs Angle",16,Yubins,250,0,10);
	
	TH2D *hYuAnIC = new TH2D("hYuAnIC","YuE vs Angle with a gate on the IC",16,Yubins,2400,0,2.4); //Energy vs angle in Yu with a gate on IC//240,0,2.4
	TH2D *hYuAnICT = new TH2D("hYuAnICT","YuE vs Angle with a gate on the IC and with target energy loss",16,Yubins,240,0,2.4);//240,0,2.4
	TH2D *hYuAnIC1 = new TH2D("hYuAnIC1","YuE vs Angle with a gate on the IC",16,Yubins,240,0,2.4);//240,0,2.4
	
	TH2D *hYuAnPID = new TH2D("hYuAnPID","YuE vs Angle with an IC and PID gate",16,Yubins,240,0,2.4); //Energy vs angle in Yu with a gate on IC and PID//240,0,2.4
	TH2D *hYuAnPIDT = new TH2D("hYuAnPIDT","YuE vs Angle with an IC and PID gate and with target energy loss",16,Yubins,240,0,2.4);//240,0,2.4
	TH2D *hYuAnPID1 = new TH2D("hYuAnPID1","YuE vs Angle with an IC and PID gate",16,Yubins,240,0,2.4);//240,0,2.4
	TH2D *hYuAnPIDYd = new TH2D("hYuAnPIDYd","YuE vs Angle with IC, Yd, Su, and PID gates",16,Yubins,240,0,2.4); //Energy vs angle in Yu with a gate on IC and PID//240,0,2.4
    
	//Q values
	TH1D *hSn = new TH1D("hSn","Sn values",800,-10,10);
	TH1D *hQval = new TH1D("hQval","Q values",800,-10,10);
	TH1D *hQvalT = new TH1D("hQvalT","Q values with target energy loss",400,-10,10);
	TH1D *hQval1 = new TH1D("hQval1","Q values",400,-10,10);
	TH2D *hQvalAn = new TH2D("hQvalAn","Q values vs Angle with energy loss",16,Yubins,400,-10,10);
	TH2D *hQvalAnT = new TH2D("hQvalAnT","Q values vs Angle with target energy loss",16,Yubins,400,-10,10);
	TH1D *hQvalYd = new TH1D("hQvalYd","Q values gated by Yd and Su",400,-10,10);
	
	TH2D *hYuEnM = new TH2D("hYuEnM","YuE1 vs YuE2 for multiplicity 2",2400,0,2.4,2400,0,2.4);
	TH2D *hYuEnMICPID = new TH2D("hYuEnMICPID","YuE1 vs YuE2 for multiplicity 2 with IC and PID gates",2400,0,2.4,2400,0,2.4);//2400,0,2.4
	
	TH2D *hYuAnPIDSec[8];//YuE vs Angle with an IC and PID gate sectorwise
    
	for(int i=0; i<8; i++){
		string namehYuAnPIDSec = Form("hYuAnPIDSec_%i",i);
		hYuAnPIDSec[i] = new TH2D(namehYuAnPIDSec.c_str(),"YuE vs Angle with an IC and PID gate",16,Yubins,240,0,2.4);//240,0,2.4
	}
    
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
	float kBeam = 112.21; //put the correct value; beam energy
	float mbeam = 12 * amu;  //mass of the beam (12Be or 12C) in MeV
	float mrecoil = 13 * amu;  //mass of the recoil (13Be or 13C) in MeV
	float mejec = 1 * amu; //mass of the proton
    

	double Qval, QvalT, Qval1, QvalYd; //Q value variable
	vector<double> *YuAngle = new vector<double>; //Yu angle variable definition
    
    
	double ringB[16], ringAl[16], ringD[16], a[16], b[16], c[16], d[16], e[16], f[16], g[16], h[16], k[16], loss, loss1, loss2, loss3;

	ifstream Bfit;
	Bfit.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/B_Dead_layer/B_fitparameters_TRIUMF.txt");
	//Bfit.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/B_Dead_layer/B_fitparameters.txt");//Micron DL

	if(Bfit.is_open()){
		Bfit.ignore(256,'\n');
		
		for(int i=0;i<16;i++){
			Bfit >> ringB[i] >> a[i] >> b[i] >> c[i];
		}
	}

	ifstream Alfit;
	Alfit.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Al_Dead_layer/Al_fitparameters_TRIUMF.txt");
	//Alfit.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Al_Dead_layer/Al_fitparameters.txt");//Micron DL

	if(Alfit.is_open()){
		Alfit.ignore(256,'\n');
		
		for(int i=0;i<16;i++){
			Alfit >> ringAl[i] >> d[i] >> e[i] >> f[i];
		}
	}
    
    
	ifstream Dfit;
	Dfit.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Target/D2_fitparameters_TRIUMF.txt");
	//Dfit.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Target/D2_fitparameters.txt");//MIcron DL
	
	if(Dfit.is_open()){
		Dfit.ignore(256,'\n');
    
		for(int i=0;i<16;i++){
			Dfit >> ringD[i] >> g[i] >> h[i] >> k[i];
		}
	}
	
	double stripGS[128], DiffSqGS[128], DiffGS[128], stripTES[128], DiffSqTES[128], DiffTES[128], DiffAvg[128];
	
	ifstream EShiftGS;
	EShiftGS.open("/home/jerome/12Be_exp/Analysis/CarbonGain/carbon_gain_GS_shift.txt");
	
	if(EShiftGS.is_open()){
		EShiftGS.ignore(256,'\n');
    
		for(int i=0;i<128;i++){
			EShiftGS >> stripGS[i] >> DiffSqGS[i] >> DiffGS[i];
		}
	}
	
	ifstream EShiftTES;
	EShiftTES.open("/home/jerome/12Be_exp/Analysis/CarbonGain/carbon_gain_TES_shift.txt");
	
	if(EShiftTES.is_open()){
		EShiftTES.ignore(256,'\n');
    
		for(int i=0;i<128;i++){
			EShiftTES >> stripTES[i] >> DiffSqTES[i] >> DiffTES[i];
			DiffAvg[i]=(DiffGS[i]+DiffTES[i])/2;
		}
	}
	
	int seed = 123;
	TRandom3* ran = new TRandom3(seed);
	double shift=0.0931323;//0.0917619;//0.0932167;
	float Q13Be0 = -2.2202;
	////////////////////////////////
	// Event by event starts here //
	////////////////////////////////
	
	        
	while(ev_num<=ev){
		if ( ev_num%10000==0 ) cout << "Current event = " << ev_num << "\r"<< flush;
		chain->GetEntry ( ev_num ); //get the current entry
		
		//Defining a strucuture YuDetector
		vector<YuDet> YuDetector;
		for(size_t i = 0; i < YuEnergy->size(); i++) {
			if(YuChannel->at(i) == 82 || YuChannel->at(i) == 96 || YuChannel->at(i) == 106 || YuChannel->at(i) == 111) continue;
			YuDet hit = {YuMul, YuChannel->at(i), YuRing->at(i), YuSector->at(i), YuEnergyRaw->at(i), YuEnergy->at(i)};
			YuDetector.push_back(hit);
		}
		
		//Defining a structure Sd1rDetector
		vector<Sd1rDet> Sd1rDetector;
		for(size_t i = 0; i < Sd1rEnergy->size(); i++) {
			Sd1rDet hit = {Sd1rMul, Sd1rChannel->at(i), Sd1rEnergyRaw->at(i), Sd1rEnergy->at(i)};
			Sd1rDetector.push_back(hit);
		}
		
		//Defining a structure Sd2rDetector
		vector<Sd2rDet> Sd2rDetector;
		for(size_t i = 0; i < Sd2rEnergy->size(); i++) {
			Sd2rDet hit = {Sd2rMul, Sd2rChannel->at(i), Sd2rEnergyRaw->at(i), Sd2rEnergy->at(i)};
			Sd2rDetector.push_back(hit);
		}
		
		//Sorting Yu
		if(!YuDetector.empty()){
			std::sort(YuDetector.begin(),YuDetector.end(),
			sortByEnergy());
		}
		
		//Sorting Sd1r
		if(!Sd1rDetector.empty()){
			std::sort(Sd1rDetector.begin(),Sd1rDetector.end(),
			sortByEnergyS3());
		}
		
		//Sorting Sd2r
		if(!Sd2rDetector.empty()){
			std::sort(Sd2rDetector.begin(),Sd2rDetector.end(),
			sortByEnergyS3());
		}
		
		if(YuDetector.size()>1) {
			if(YuDetector[0].energy>0){
				hYuEnM->Fill(YuDetector[0].energy,YuDetector[1].energy);
				if(Sd2rEnergy->size()>0 && Sd1rEnergy->size()>0){
					//if (pidcut->IsInside(Sd2rEnergyRaw->at ( 0 ),Sd1rEnergyRaw->at ( 0 ) ) && ICEnergyRaw>1500 && ICEnergyRaw<2200 ){
					if (pidcut->IsInside(Sd2rEnergy->at ( 0 ),Sd1rEnergy->at ( 0 ) ) && ICEnergyRaw>620 && ICEnergyRaw<1100 ){
						hYuEnMICPID->Fill(YuDetector[0].energy,YuDetector[1].energy);
					}
				}
			}
		}
		
		//Sd1r vs Sd2r calibrated
		if(Sd1rDetector.size()>0 && Sd2rDetector.size()>0 && Sd1rMul==1 && Sd2rMul==1){
			hCSd1rSd2r->Fill(Sd2rDetector[0].energy,Sd1rDetector[0].energy); //Calibrated Sd1r vs Sd2r energy with no gates
			if ( ICEnergyRaw>1500 && ICEnergyRaw<2200){
				hCSd1rSd2rIC->Fill ( Sd2rDetector[0].energy,Sd1rDetector[0].energy );  //Calibrated Sd1r vs Sd2r energy with IC gate
				if (TMath::Abs(Sd2rDetector[0].channel-Sd1rDetector[0].channel)<2){
					hCSd1rSd2rICCut->Fill ( Sd2rDetector[0].energy,Sd1rDetector[0].energy);
				} //Calibrated Sd1r vs Sd2r energy with IC gate and ring vs ring correlation
				if ( YdEnergy->size() >0 && YdMul==1 && YdEnergy->at( 0 ) >0){
					hCSd1rSd2rYdIC->Fill(Sd2rDetector[0].energy,Sd1rDetector[0].energy);//Calibrated Sd1r vs Sd2r energy with IC gate and Yd gate
					if (TMath::Abs(Sd2rDetector[0].channel-Sd1rDetector[0].channel)<2){
						hCSd1rSd2rYdICCut->Fill(Sd2rDetector[0].energy,Sd1rDetector[0].energy);
					} //Calibrated Sd1r vs Sd2r energy with IC and Yd gate and ring vs ring correlation
				} //end of Yd gates
			} //end of IC gates
		} //end of Sd1r vs Sd2r calibrated
		
		//fill in the Yu detector
		if(YuDetector.size()>0){
			//if (YuDetector[0].energy>0.20 && ( YuChannel->at(0)!=82 && YuChannel->at(0)!=96 && YuChannel->at(0)!=106 && YuChannel->at(0)!=111 ))
			if (YuDetector[0].energy>0.5){
				//calculate the angles for the Yu detector
				yuM = Ydr0+ ( YuDetector[0].ring*YChWidth )+YChMiddle;
				TVector3 YuposM ( 0,yuM,Yuz); //shifting the detecor
				YuposM.SetPhi ( TMath::Pi() /2-YuDetector[0].sector *TMath::Pi() /4 ); //Pi/2 because the center of the sector is at 90degrees, Pi/4 is because there is 8 sectors so it is 2Pi/8
				//YuthetaM = beam.Angle ( YuposM ) *TMath::RadToDeg();
				YuthetaM=ran->Uniform(Yubins[15-YuDetector[0].ring],Yubins[15-YuDetector[0].ring+1]);
				YuAngle->push_back(YuthetaM);
							
				//Adding the energy lost by the protons through the target and the deadlayers of the Yu detector (11/15/2018)  
				for(int RingNumber=0; RingNumber<16; RingNumber++){	    
					for(int j=0; j<8; j++){ 
						
						if(YuChannel->at(0)==j*16+RingNumber){
				
							loss1=(-b[RingNumber]-TMath::Sqrt(pow(b[RingNumber],2)-4*(a[RingNumber]-YuDetector[0].energy)*c[RingNumber]))/(2*(a[RingNumber]-YuDetector[0].energy));// Boron dead layer
							loss2=(-e[RingNumber]-TMath::Sqrt((e[RingNumber]*e[RingNumber])-4*(d[RingNumber]-(YuDetector[0].energy+loss1))*f[RingNumber]))/(2*(d[RingNumber]-(YuDetector[0].energy+loss1)));// Aluminium dead layer
							loss3=(-h[RingNumber]-TMath::Sqrt((h[RingNumber]*h[RingNumber])-4*(g[RingNumber]-(YuDetector[0].energy+loss1+loss2))*k[RingNumber]))/(2*(g[RingNumber]-(YuDetector[0].energy+loss1+loss2)));// Target
							loss=loss1+loss2+loss3; //
							//cout << ev_num << "	" << RingNumber << "		" << YuEnergy->at ( 0 ) << loss1 << loss2 << loss3 << loss << endl;
							hYuEn->Fill(YuDetector[0].energy+loss);
							hYuEnT->Fill(YuDetector[0].energy+loss3);
							hYuAn->Fill(YuthetaM,YuDetector[0].energy+loss);
							hYuAnT->Fill(YuthetaM,YuDetector[0].energy+loss3);
							//if(ICEnergyRaw>1500 && ICEnergyRaw<2200){
							if(ICEnergyRaw>620 && ICEnergyRaw<1100){
								hYuAnIC->Fill(YuthetaM,YuDetector[0].energy+loss);
								hYuAnICT->Fill(YuthetaM,YuDetector[0].energy+loss3);
								hYuEnIC->Fill(YuEnergy->at (0)+loss);
								
							} //IC gate for C 1500 < E < 2200; for Be 620< E < 1100
							
							if(Sd2rEnergy->size()>0 && Sd1rEnergy->size()>0){
								//if(pidcut->IsInside(Sd2rEnergyRaw->at(0),Sd1rEnergyRaw->at(0)) && ICEnergyRaw>1500 && ICEnergyRaw<2200 ){
								if(pidcut->IsInside(Sd2rEnergy->at(0),Sd1rEnergy->at(0)) && ICEnergyRaw>620 && ICEnergyRaw<1100 ){
									hYuEnPID->Fill(YuDetector[0].energy+loss);
									hYuAnPID->Fill(YuthetaM,YuDetector[0].energy+loss-shift);
									hYuAnPIDT->Fill(YuthetaM,YuDetector[0].energy+loss3);
									hYuAnPIDSec[YuDetector[0].sector]->Fill(YuthetaM,YuDetector[0].energy+loss-shift);
									Qval = ( 1+mejec/mrecoil ) * ( YuDetector[0].energy+loss-shift) - ( 1 - mbeam/mrecoil ) * ( kBeam ) - 2 * TMath::Sqrt ( mbeam*mejec* ( YuDetector[0].energy+loss-shift) * ( kBeam ) ) / mrecoil * TMath::Cos ( YuthetaM * TMath::Pi() / 180. );
									hSn->Fill(Q13Be0-Qval);
									hQval->Fill ( Qval );
									hQvalAn->Fill(YuthetaM,Qval);
									QvalT = ( 1+mejec/mrecoil ) * ( YuDetector[0].energy+loss3) - ( 1 - mbeam/mrecoil ) * ( kBeam ) - 2 * TMath::Sqrt ( mbeam*mejec* ( YuDetector[0].energy+loss3) * ( kBeam ) ) / mrecoil * TMath::Cos ( YuthetaM * TMath::Pi() / 180. );
									hQvalT->Fill ( QvalT);
									hQvalAnT->Fill(YuthetaM,QvalT);
									if(YdEnergy->size()==0 && SurEnergyRaw->size()==0 && SusEnergyRaw->size()==0){
										hYuAnPIDYd->Fill(YuthetaM,YuDetector[0].energy+loss-shift);
										QvalYd = ( 1+mejec/mrecoil ) * ( YuDetector[0].energy+loss-shift) - ( 1 - mbeam/mrecoil ) * ( kBeam ) - 2 * TMath::Sqrt ( mbeam*mejec* ( YuDetector[0].energy+loss-shift) * ( kBeam ) ) / mrecoil * TMath::Cos ( YuthetaM * TMath::Pi() / 180. );
										hQvalYd->Fill ( QvalYd );
									}
								}	
							}
						}
					}//End of j loop
				}//End of Ringnumber loop
				
				hYuAn1->Fill(YuthetaM,YuDetector[0].energy); //all, no cuts
				//if(ICEnergyRaw>1500 && ICEnergyRaw<2200){
				if(ICEnergyRaw>620 && ICEnergyRaw<1100){
					hYuAnIC1->Fill(YuthetaM,YuEnergy->at (0));
					hYuEnIC1->Fill(YuEnergy->at (0));
				} //IC gate for C 1500 < E < 2200; for Be 620< E < 1100
				hYuEn1->Fill(YuDetector[0].energy);
				hYuMul->Fill(YuMul);
								
				if(Sd2rEnergy->size()>0 && Sd1rEnergy->size()>0){
					//if (pidcut->IsInside(Sd2rEnergyRaw->at ( 0 ),Sd1rEnergyRaw->at ( 0 ) ) && ICEnergyRaw>1500 && ICEnergyRaw<2200 ){
					if (pidcut->IsInside(Sd2rEnergy->at ( 0 ),Sd1rEnergy->at ( 0 ) ) && ICEnergyRaw>620 && ICEnergyRaw<1100 ){
						hYuAnPID1->Fill ( YuthetaM,YuDetector[0].energy );
						Qval1 = ( 1+mejec/mrecoil ) * ( YuDetector[0].energy ) - ( 1 - mbeam/mrecoil ) * ( kBeam ) - 2 * TMath::Sqrt ( mbeam*mejec* ( YuDetector[0].energy ) * ( kBeam ) ) / mrecoil * TMath::Cos ( YuthetaM * TMath::Pi() / 180. );
						hQval1->Fill ( Qval1 );
					}
				}
			}
		} 
		ev_num++;
		
	} //end of the main while loop
    
  
	YuAngle->clear();
    
	f_out->cd();
    
	//write the histograms in the output root file
	//Yu
	hCSd1rSd2r->Write();
	hCSd1rSd2rIC->Write(); //calibrated ring vs ring with IC cut
	hCSd1rSd2rICCut->Write(); //calibrated ring vs ring with IC cut and diagonal taken into account
	hCSd1rSd2rYdIC->Write();
	hCSd1rSd2rYdICCut->Write();
	hYuAn->Write();
	hYuAnT->Write();
	hYuAn1->Write();
	hYuAnIC->Write();
	hYuAnICT->Write();
	hYuAnIC1->Write();
	hYuAnPID->Write();
	hYuAnPIDT->Write();
	hYuAnPID1->Write();
	hYuAnPIDYd->Write();
	hYuEn->Write();
	hYuEnT->Write();
	hYuEnIC->Write();
	hYuEnIC1->Write();
	hYuMul->Write();
	hYuEnPID->Write();
	hYuEn1->Write();
	hSn->Write();
	hQval->Write();
	hQvalT->Write();
	hQval1->Write();
	hQvalYd->Write();
	hQvalAn->Write();
	hQvalAnT->Write();
	hYuEnM->Write();
	hYuEnMICPID->Write();
	for(int i=0; i<8; i++){hYuAnPIDSec[i]->Write();}
	
	f_out->Close();
	
	return 0;
	
} //end of the program