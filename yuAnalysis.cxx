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
} S3Det;
	
struct sortByEnergyYY1 {
	inline bool operator() (const YuDet& EnYY1_1,
							const YuDet& EnYY1_2){
		return (EnYY1_1.energy>EnYY1_2.energy);
	}
};

struct sortByEnergyS3 {
	inline bool operator() (const S3Det& EnS3_1,
							const S3Det& EnS3_2){
		return (EnS3_1.energy>EnS3_2.energy);
	}
};

std::vector<YuDet> CheckChargeSharing(std::vector<YuDet> detect) {
	std::vector<YuDet> newDetector;
	
	int sectorHits[8];
	for(int i = 0; i < 8; i++) {
		sectorHits[i] = 0;
	}
	for(auto det : detect) {
		sectorHits[det.sector]++;
	}
	
	for(int i = 0; i < 8; i++) {
		if(sectorHits[i] > 1) printf("Sector %d has %d hits\n", i, sectorHits[i]);
	}
	
	return newDetector;
}


double yuAnalysis(){
    
    bool TT43 = 0;//Target thickness 43 um
    bool TT45 = 0;//Target thickness 45 um
    bool TT47 = 0;//Target thickness 47 um
    bool TT49_99 = 1;//Target thickness 49.99um
    
    float TThickness;
    
    if(TT43){
        TThickness = 43;
    }else if(TT45){
        TThickness = 45;
    }else if(TT47){
        TThickness = 47;
    }else if(TT49_99){
        TThickness = 49.99;
    }   
    
	//open the output file
    
    //TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/C_pedestal_TRIUMF_DL_NoOffset_Shift0_0883_Yu_TargetDistance80_88mm_TargetMiddle.root","RECREATE");
    //TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/C_pedestal_TRIUMF_DL_NoOffset_NoShift_Yu_TargetDistance80_88mm_TargetMiddle.root","RECREATE");
    //TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/C_pedestal_TRIUMF_DL_NoOffset_NoShift_Yu_TargetDistance80_88mm_TargetMiddle.root","RECREATE");

	//TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/C_pedestal_TRIUMF_DL_NoOffset_Shift0_0883_Yu_TargetDistance80_88mm_TargetThickness49.99um_SdrNewInBeamCal_TargetFront.root","RECREATE"); //change the path and name accordingly
    
    //TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/C_pedestal_TRIUMF_DL_NoOffset_NoShift_Yu_TargetDistance85mm_TargetThickness43um_SdrNewInBeamCal_TargetMiddle.root","RECREATE"); //change the path and name accordingly
    
    //TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withoutTarget/C_pedestal_TRIUMF_DL_NoTarget_NoOffset_Shift0_0883_Yu_TargetDistance80_88mm_TargetThickness_SdrNewInBeamCal_TargetMiddle.root","RECREATE"); //No target
    
    
    
    
	//Open the input files
	TChain *chain = new TChain ( "AutoTree" ); //chain the desired input files
	
    EnergyLoss* protonELB = new EnergyLoss("Proton_Boron.dat");
    EnergyLoss* protonELA = new EnergyLoss("Proton_Aluminum.dat");
	EnergyLoss* protonELD = new EnergyLoss("Proton_DeuteriumTarget.dat");
    
    protonELB->AddBackHigh(15.);
    protonELA->AddBackHigh(15.);
    protonELD->AddBackHigh(15.);
	
	std::map<int, double> YuAngleMap;
	YuAngleMap[0] = 31.68434901;
	
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
	
	/*for(int run_num=5021;run_num<5113;run_num++){
		if(run_num==5026||run_num==5040||run_num==5043||run_num==5046||run_num==5047||run_num==5059||run_num==5062||run_num==5063||run_num==5093||run_num==5099||run_num==5101){
			continue;
		}else{
			string f_name=Form("/home/jerome/12Be_exp/Analysis/Be/AlphaOnly/decode_Yupedestal_%i.root",run_num);
			chain->Add(f_name.c_str());
		}
	}
	
	for(int run_num=5115;run_num<5133;run_num++){
		string f_name=Form("/home/jerome/12Be_exp/Analysis/Be/AlphaOnly/decode_Yupedestal_%i.root",run_num);
		chain->Add(f_name.c_str());
	}
	
	for(int run_num=5180;run_num<5223;run_num++){
		if(run_num==5187||run_num==5200){
			continue;
		}else{
			string f_name=Form("/home/jerome/12Be_exp/Analysis/Be/AlphaOnly/decode_Yupedestal_%i.root",run_num);
			chain->Add(f_name.c_str());
		}
	}*/
    
    
	
	//Carbon data with target
	
	chain->Add("/home/jerome/12Be_exp/Analysis/ReducedData/CarbonReduced/C_13CcutFull_Data.root");
	
	/*
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/decodeNew_TDL_5001.root" ); //change the path to the files and the file name accordingly
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/decodeNew_TDL_5002.root" );
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/decodeNew_TDL_5003.root" );
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/decodeNew_TDL_5004.root" );
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/decodeNew_TDL_5006.root" );
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/decodeNew_TDL_5007.root" );
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/decodeNew_TDL_5009.root" );
	*/
    
    /*
    //Carbon data with target
	chain->Add ( "/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/decodeNewC_SdrNewInBeamCal_Target_TDL_YuPedestal5001.root" ); //change the path to the files and the file name accordingly
	chain->Add ( "/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/decodeNewC_SdrNewInBeamCal_Target_TDL_YuPedestal5002.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/decodeNewC_SdrNewInBeamCal_Target_TDL_YuPedestal5003.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/decodeNewC_SdrNewInBeamCal_Target_TDL_YuPedestal5004.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/decodeNewC_SdrNewInBeamCal_Target_TDL_YuPedestal5006.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/decodeNewC_SdrNewInBeamCal_Target_TDL_YuPedestal5007.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/decodeNewC_SdrNewInBeamCal_Target_TDL_YuPedestal5009.root" );
	*/
    
    
    
    
    
    
    
    /*
    chain->Add ( "/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/decodeNewC_Sd1NewAlphaCal_Target_TDL_YuPedestal5001.root" );
    chain->Add ( "/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/decodeNewC_Sd1NewAlphaCal_Target_TDL_YuPedestal5002.root" );
    chain->Add ( "/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/decodeNewC_Sd1NewAlphaCal_Target_TDL_YuPedestal5003.root" );
    chain->Add ( "/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/decodeNewC_Sd1NewAlphaCal_Target_TDL_YuPedestal5004.root" );
    chain->Add ( "/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/decodeNewC_Sd1NewAlphaCal_Target_TDL_YuPedestal5006.root" );
    chain->Add ( "/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/decodeNewC_Sd1NewAlphaCal_Target_TDL_YuPedestal5007.root" );
    chain->Add ( "/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withTarget/decodeNewC_Sd1NewAlphaCal_Target_TDL_YuPedestal5009.root" );
	*/
    
    // Carbon no target data
    
    //chain->Add ( "/home/jerome/12Be_exp/Processed_files/C_notarget/decode4992.root" );
    
    
    //chain->Add ( "/home/jerome/12Be_exp/Analysis/TargetThickness/FromCBeam/withoutTarget/decodeNewC_SdrNewInBeamCal_NoTarget_TDL_YuPedestal4992.root" );
	
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
	TFile *f_out = new TFile(Form("/home/jerome/12Be_exp/Analysis/BeamOffset/carbon/C_pedestal_TRIUMF_DL_BeamOffset_-3_4_TargetDistance%.0fmm_TargetThickness%.2fum.root",0-Yuz,TThickness),"RECREATE");
	
	/*TFile *f_cut = TFile::Open("/home/jerome/12Be_exp/scripts/13Ccut.root"); //PID cut for C
	TCutG *pidcut = (TCutG*) f_cut->Get("CUTG"); // for C
	*/
	
	TFile *f_cut = TFile::Open("/home/jerome/12Be_exp/Analysis/CarbonGain/13CcutFull.root"); //Full PID cut for C
	TCutG *pidcut = (TCutG*) f_cut->Get("CcutCalSd1rSd2rFull"); // for C
	
    
	/*TFile *f_cut = TFile::Open("/home/jerome/12Be_exp/Analysis/Be/BeCut1.root"); //PID cut for Be
	TCutG *pidcut = (TCutG*) f_cut->Get("BeCut_CalSd1rSd2r"); */ // for Be
	
	//definition of histograms
	//calculate the variable bins for the Yu detector to plot the angles
	double Yubins[17]={0};
	
	for (int i=0; i<17; i++){
		Yubins[i]=180-TMath::RadToDeg()*atan2((50+(16-i)*4.94),(0-Yuz));
		//cout << Yubins[i] << endl;
	}
	
	TH2D *hSd1rSd2r = new TH2D("hSd1rSd2r","Sd1r vs Sd2r",500,0,5000,500,0,5000); //non calibrated Sd1r vs Sd2r
    TH2D *hSd1rSd2s = new TH2D("hSd1rSd2s","Sd1r vs Sd2s",500,0,5000,500,0,5000); //non calibrated Sd1r vs Sd2s
	TH2D *hCSd1rSd2r = new TH2D("hCSd1rSd2r","Calibrated Sd1r vs Sd2r",1500,0,150,600,0,25); // calibrated Sd1r vs Sd2r
	TH2D *hCSd1rSd2rIC = new TH2D("hCSd1rSd2rIC","Cal Sd1r vs Sd2r, gated by IC",1500,0,150,600,0,25); // calibrated Sd1r vs Sd2r and gated by IC
	TH2D *hCSd1rSd2rICCut = new TH2D("hCSd1rSd2rICCut","Cal Sd1r vs Sd2r, gated by IC gated by the channel corr",1500,0,150,600,0,25); // calibrated Sd1r vs Sd2r gated by IC and required to be
	TH2D *hCSd1rSd2rYdIC = new TH2D("hCSd1rSd2rYdIC","Cal Sd1r vs Sd2r, gated by Yd & IC",1500,0,150,600,0,25); // calibrated Sd1r vs Sd2r and gated by Yd & IC
	TH2D *hCSd1rSd2rYdICCut = new TH2D("hCSd1rSd2rYdICCut","Cal Sd1r vs Sd2r, gated by Yd & IC & chann corr",1500,0,150,600,0,25); // calibrated Sd1r vs Sd2r gated by Yd, IC and required to be correlated with each other
	
    TH1D *hCSd1rEn = new TH1D("hCSd1rEn","Calibrated SD1r Energy",500,0,50);
    TH1D *hCSd1r_0_1_2_En = new TH1D("hCSd1r_0_1_2_En","Calibrated SD1r Energy (Rings 0, 1, and 2)",500,0,50);
    
    TH1D *hCSd1sEn = new TH1D("hCSd1sEn","Calibrated SD1s Energy",500,0,50);
    TH1D *hCSd2sEn = new TH1D("hCSd2sEn","Calibrated SD2s Energy",500,0,100);
    
    TH1D *hCSd2rEn = new TH1D("hCSd2rEn","Calibrated SD2r Energy",500,0,120);
    TH1D *hCSd2r_0_1_2_En = new TH1D("hCSd2r_0_1_2_En","Calibrated SD2r Energy (Rings 0, 1, and 2)",500,0,120);
    
  	TH1D *hYuMul = new TH1D("hYuMul","hYuMul",10,0,10);
	
	TH1D *hYuEn = new TH1D("hYuEn","hYuEn",2500,0,10); //energy singles in Yu
	TH1D *hYuEnT = new TH1D("hYuEnT","hYuEn with target energy loss",2500,0,10);
	TH1D *hYuEn1 = new TH1D("hYuEn1","hYuEn1",2500,0,10);
	
	TH1D *hYuEnPID = new TH1D("hYuEnPID","hYuEnPID",6000,0,6);//6000,0,6; 2400,0,2.4
	
	TH1D *hYuEnIC = new TH1D("hYuEnIC","hYuEnIC",6000,0,6);
	TH1D *hYuEnIC1 = new TH1D("hYuEnIC1","hYuEnIC1",6000,0,6);
	
	TH2D *hYuAn = new TH2D("hYuAn","YuE vs Angle",16,Yubins,250,0,10); //Energy vs angle in Yu no gates
	TH2D *hYuAnT = new TH2D("hYuAnT","YuE vs Angle with target energy loss",16,Yubins,250,0,10);
	TH2D *hYuAn1 = new TH2D("hYuAn1","YuE vs Angle",16,Yubins,250,0,10);
	
	TH2D *hYuAnIC = new TH2D("hYuAnIC","YuE vs Angle with a gate on the IC",16,Yubins,6000,0,6); //Energy vs angle in Yu with a gate on IC//240,0,2.4
	TH2D *hYuAnICT = new TH2D("hYuAnICT","YuE vs Angle with a gate on the IC and with target energy loss",16,Yubins,6000,0,6);//240,0,2.4
	TH2D *hYuAnIC1 = new TH2D("hYuAnIC1","YuE vs Angle with a gate on the IC",16,Yubins,6000,0,6);//240,0,2.4
	
	TH2D *hYuAnPID = new TH2D("hYuAnPID","YuE vs Angle with an IC and PID gate",16,Yubins,6000,0,6); //Energy vs angle in Yu with a gate on IC and PID//240,0,2.4
	TH2D *hYuAnPIDT = new TH2D("hYuAnPIDT","YuE vs Angle with an IC and PID gate and with target energy loss",16,Yubins,6000,0,6);//240,0,2.4
	TH2D *hYuAnPID1 = new TH2D("hYuAnPID1","YuE vs Angle with an IC and PID gate",16,Yubins,6000,0,6);//240,0,2.4
    
	TH2D *hYuAnPIDLeft = new TH2D("hYuAnPIDLeft","YuE vs Angle with an IC and PID gate, left side of Yu",16,Yubins,6000,0,6);
	TH2D *hYuAnPIDRight = new TH2D("hYuAnPIDRight","YuE vs Angle with an IC and PID gate, right side of Yu",16,Yubins,6000,0,6);
	TH2D *hYuAnPIDTop = new TH2D("hYuAnPIDTop","YuE vs Angle with an IC and PID gate, top part of Yu",16,Yubins,6000,0,6);
	TH2D *hYuAnPIDBottom = new TH2D("hYuAnPIDBottom","YuE vs Angle with an IC and PID gate, bottom part of Yu",16,Yubins,6000,0,6);
	TH2D *hYuAnPID315_45 = new TH2D("hYuAnPID315_45","YuE vs Angle with an IC and PID gate, 315-45",16,Yubins,6000,0,6);
	TH2D *hYuAnPID45_135 = new TH2D("hYuAnPID45_135","YuE vs Angle with an IC and PID gate, 45-135",16,Yubins,6000,0,6);
	TH2D *hYuAnPID135_225 = new TH2D("hYuAnPID135_225","YuE vs Angle with an IC and PID gate, 135-225",16,Yubins,6000,0,6);
	TH2D *hYuAnPID225_315 = new TH2D("hYuAnPID225_315","YuE vs Angle with an IC and PID gate, 225-315",16,Yubins,6000,0,6);
	TH2D *hYuAnPID0_90 = new TH2D("hYuAnPID0_90","YuE vs Angle with an IC and PID gate, 0-90",16,Yubins,6000,0,6);
	TH2D *hYuAnPID90_180 = new TH2D("hYuAnPID90_180","YuE vs Angle with an IC and PID gate, 90-180",16,Yubins,6000,0,6);
	TH2D *hYuAnPID180_270 = new TH2D("hYuAnPID180_270","YuE vs Angle with an IC and PID gate, 180-270",16,Yubins,6000,0,6);
	TH2D *hYuAnPID270_360 = new TH2D("hYuAnPID270_360","YuE vs Angle with an IC and PID gate, 270-360",16,Yubins,6000,0,6);
	TH2D *hYuAnPIDRs0_7 = new TH2D("hYuAnPIDRs0_7","YuE vs Angle with an IC and PID gate, rings 0 to 7",16,Yubins,6000,0,6);
	//TH2D *hYuAnPIDR4_7 = new TH2D("hYuAnPIDR4_7","YuE vs Angle with an IC and PID gate, rings 4 to 7",16,Yubins,6000,0,6);
	//TH2D *hYuAnPIDR8_11 = new TH2D("hYuAnPIDR8_11","YuE vs Angle with an IC and PID gate, rings 8 to 11",16,Yubins,6000,0,6);
	TH2D *hYuAnPIDRs8_15 = new TH2D("hYuAnPIDRs12_15","YuE vs Angle with an IC and PID gate, rings 8 to 15",16,Yubins,6000,0,6);
	
    TH2D *hICSd1rSd2r = new TH2D("hICSd1rSd2r","IC Ch. No. vs Sd1r+Sd2r Energy in MeV",600,0,120,400,0,4000);
    
	//Q values
	TH1D *hQval = new TH1D("hQval","Q values",100,-3.5,4);
	TH1D *hQvalT = new TH1D("hQvalT","Q values with target energy loss",100,-3.5,4);
	TH1D *hQval1 = new TH1D("hQval1","Q values",160,-4,4);
	TH2D *hQvalAn = new TH2D("hQvalAn","Q values vs Angle with energy loss",16,Yubins,160,-4,4);
	TH2D *hQvalAnT = new TH2D("hQvalAnT","Q values vs Angle with target energy loss",16,Yubins,160,-4,4);
	
	TH1D *hQvalLeft = new TH1D("hQvalLeft","Q values, left side of Yu",160,-4,4);
	TH1D *hQvalRight = new TH1D("hQvalRight","Q values, right side of Yu",160,-4,4);
	TH1D *hQvalTop = new TH1D("hQvalTop","Q values, top part of Yu",160,-4,4);
	TH1D *hQvalBottom = new TH1D("hQvalBottom","Q values, bottom part of Yu",160,-4,4);
	TH1D *hQval315_45 = new TH1D("hQval315_45","Q values, 315-45",160,-4,4);
	TH1D *hQval45_135 = new TH1D("hQval45_135","Q values, 45-135",160,-4,4);
	TH1D *hQval135_225 = new TH1D("hQval135_225","Q values, 135-225",160,-4,4);
	TH1D *hQval225_315 = new TH1D("hQval225_315","Q values, 225-315",160,-4,4);
	TH1D *hQval0_90 = new TH1D("hQval0_90","Q values, 0-90",160,-4,4);
	TH1D *hQval90_180 = new TH1D("hQval90_180","Q values, 90-180",160,-4,4);
	TH1D *hQval180_270 = new TH1D("hQval180_270","Q values, 180-270",160,-4,4);
	TH1D *hQval270_360 = new TH1D("hQval270_360","Q values, 270-360",160,-4,4);
	TH1D *hQvalRs0_7 = new TH1D("hQvalRs0_7","Q values, rings 0 to 7",160,-4,4);
    
    TH1D *hQvalR[16];
    for (int i=0;i<16;i++){
        hQvalR[i] = new TH1D(Form("hQvalR%d",i),Form("Q values, ring %d",i),160,-4,4);
    }
    
    TH1D *hQvalS[8];
    for (int i=0;i<8;i++){
        hQvalS[i] = new TH1D(Form("hQvalS%d",i),Form("Q values, sector %d",i),160,-4,4);
    }
    
    TH1D *hQvalRs0_3 = new TH1D("hQvalRs0_3","Q values, rings 0 to 3",160,-4,4);
	TH1D *hQvalRs4_7 = new TH1D("hQvalRs4_7","Q values, rings 4 to 7",160,-4,4);
	TH1D *hQvalRs8_11 = new TH1D("hQvalRs8_11","Q values, rings 8 to 11",160,-4,4);
    TH1D *hQvalRs12_15 = new TH1D("hQvalRs12_15","Q values, rings 12 to 15",160,-4,4);
    
	TH1D *hQvalRs8_15 = new TH1D("hQvalRs8_15","Q values, rings 8 to 15",160,-4,4);
    TH1D *hQvalRs4_12 = new TH1D("hQvalRs4_12","Q values, rings 4 to 12",160,-4,4);
    TH1D *hQvalRs12_4 = new TH1D("hQvalRs12_4","Q values, rings 12 to 4",160,-4,4);
	
	TH2D *hYuEnM = new TH2D("hYuEnM","YuE1 vs YuE2 for multiplicity 2",6000,0,6,6000,0,6);
	TH2D *hYuEnMICPID = new TH2D("hYuEnMICPID","YuE1 vs YuE2 for multiplicity 2 with IC and PID gates",6000,0,6,6000,0,6);//2400,0,2.4
	
	TH2D *hYuAnPIDSec[8];//YuE vs Angle with an IC and PID gate sectorwise
	TH2D *hYuAnPIDRing[16];//YuE vs Angle with an IC and PID gate ringwise
	TH2D *hYuAnPIDR0_7[8];
    TH2D *hYuAnPIDR8_15[8];
    TH1D *hQvalR0_7[8];
    TH1D *hQvalR8_15[8];
    
    
	for(int i=0; i<8; i++){
		string namehYuAnPIDSec = Form("hYuAnPIDSec_%i",i);
		hYuAnPIDSec[i] = new TH2D(namehYuAnPIDSec.c_str(),"YuE vs Angle with an IC and PID gate",16,Yubins,6000,0,6);//240,0,2.4
	}
	
	for(int i=0; i<16; i++){
		string namehYuAnPIDRing = Form("hYuAnPIDRing_%i",i);
		hYuAnPIDRing[i] = new TH2D(namehYuAnPIDRing.c_str(),"YuE vs Sector Number with an IC and PID gate for each ring",8,0,8,6000,0,6);//240,0,2.4
	}
    
    for(int i=0; i<8; i++){
		string namehYuAnPIDR0_7 = Form("hYuAnPIDR0_7_%i",i);
		hYuAnPIDR0_7[i] = new TH2D(namehYuAnPIDR0_7.c_str(),Form("YuE vs Angle with an IC and PID gate for rings 0-7 with ring_%i removed",i),16,Yubins,6000,0,6);//240,0,2.4
        string namehYuAnPIDR8_15 = Form("hYuAnPIDR8_15_%i",i+8);
		hYuAnPIDR8_15[i] = new TH2D(namehYuAnPIDR8_15.c_str(),Form("YuE vs Angle with an IC and PID gate for rings 8-15 with ring_%i removed",i+8),16,Yubins,6000,0,6);//240,0,2.4
        string namehQvalR0_7 = Form("hQvalR0_7_%i",i);
		hQvalR0_7[i] = new TH1D(namehQvalR0_7.c_str(),Form("Q values for rings 0-7 with ring_%i removed",i),160,-4,4);//240,0,2.4
        string namehQvalR8_15 = Form("hQvalR8_15_%i",i+8);
		hQvalR8_15[i] = new TH1D(namehQvalR8_15.c_str(),Form("Q values for rings 8-15 with ring_%i removed",i+8),160,-4,4);//240,0,2.4
	}
	
    
	//start reading the tree
	int ev = chain->GetEntries(); //get the total number of entries
	cout << "Total number of events =" << ev << endl;
	
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
	
	float kBeam; //put the correct value; beam energy 110.22 MeV (49.99 um) at the center of the target, 110.36(43 um); 111.22 MeV at the front of the target
    
	if(TT43){
        kBeam = 110.36;
    }else if(TT45){
        kBeam = 110.32;
    }else if(TT47){
        kBeam = 110.28;
    }else if(TT49_99){
        kBeam = 110.22;
    }
	
	float mbeam = 12 * amu;  //mass of the beam (12Be or 12C) in MeV
	float mrecoil = 13 * amu;  //mass of the recoil (13Be or 13C) in MeV
	float mejec = 1 * amu; //mass of the proton
    
    cout << TThickness << " " << kBeam << endl;
	double Qval, QvalT, Qval1;; //Q value variable
	vector<double> *YuAngle = new vector<double>; //Yu angle variable definition
    
    
	double ringB[16], ringAl[16], ringD[16], a[16], b[16], c[16], d[16], e[16], f[16], g[16], h[16], k[16], loss, loss1, loss2, loss3;

	ifstream Bfit;
	Bfit.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/B_Dead_layer/B_fitparameters_TRIUMF.txt");

	if(Bfit.is_open()){
		Bfit.ignore(256,'\n');
		
		for(int i=0;i<16;i++){
			Bfit >> ringB[i] >> a[i] >> b[i] >> c[i];
		}
	}

	ifstream Alfit;
	Alfit.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Al_Dead_layer/Al_fitparameters_TRIUMF.txt");

	if(Alfit.is_open()){
		Alfit.ignore(256,'\n');
		
		for(int i=0;i<16;i++){
			Alfit >> ringAl[i] >> d[i] >> e[i] >> f[i];
		}
	}
    
    
	ifstream Dfit;
	Dfit.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Target/D2_fitparameters_TRIUMF.txt");
	
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
	
	//double x, y;
    double phi[8], r[16], R[8][16], DetAngle[8][16];
    phi[0]=22.5;
    for (int i=1;i<8;i++){
        phi[i]=phi[0]+45*i;
        //cout << phi[i] << endl;
    }
    
    for (int i=0;i<16;i++){
        r[i]=50+((129.-50.)/16)*(0.5+i);
        //cout << r[i] << endl;
    }
    
    int seed = 123;
	TRandom3* ran = new TRandom3(seed);
	
	for (int i=0;i<8;i++){
        for (int j=0;j<16;j++){
            double x_yu = r[j]*sin(phi[i]*M_PI/180.);
            double y_yu = r[j]*cos(phi[i]*M_PI/180.);
            double x_target = -3;
            double y_target = 4;
            TVector3 vec1(0, 0, 0-Yuz);
            TVector3 vec2(x_yu - x_target, y_yu - y_target, Yuz);
            double angle = vec1.Angle(vec2);
            // std::cout << i << '\t' << j << '\t' << angle*180./M_PI << std::endl;
            //R[i][j][xindex][yindex]=TMath::Sqrt(pow((x[xindex]-A),2)+pow((y[yindex]-B),2));
            //DetAngle[i][j][xindex][yindex]=180-TMath::RadToDeg()*TMath::ATan(R[i][j][xindex][yindex]/(0-Yuz));
            DetAngle[i][j] = angle*180./M_PI;
            //cout << DetAngle[i][j][xindex][yindex] << endl;
        }
    }
	

	double shift=0;//0.0883486;//0.0904;//0.0924;//0.154;
	double YuEnergyLoss, YuEnergyShift;
	////////////////////////////////
	// Event by event starts here //
	////////////////////////////////
	
	for(int ev_num = 0; ev_num < ev; ev_num++) {
		if ( ev_num%50000==0 ) cout << "Current event = " << ev_num << "\r"<< flush;
		chain->GetEntry ( ev_num ); //get the current entry
				
		//Defining a structure YuDetector
		vector<YuDet> YuDetector;
		for(size_t i = 0; i < YuEnergy->size(); i++) {
			if(YuChannel->at(i) == 82 || YuChannel->at(i) == 96 || YuChannel->at(i) == 106 || YuChannel->at(i) == 111) continue;
			YuDet hit = {YuMul, YuChannel->at(i), YuRing->at(i), YuSector->at(i), YuEnergyRaw->at(i), YuEnergy->at(i)};
			YuDetector.push_back(hit);
		}
		
		//Defining a structure Sd1rDetector
		vector<S3Det> Sd1rDetector;
		for(size_t i = 0; i < Sd1rEnergy->size(); i++) {
			S3Det hit = {Sd1rMul, Sd1rChannel->at(i), Sd1rEnergyRaw->at(i), Sd1rEnergy->at(i)};
			Sd1rDetector.push_back(hit);
		}
		
		//Defining a structure Sd2rDetector
		vector<S3Det> Sd2rDetector;
		for(size_t i = 0; i < Sd2rEnergy->size(); i++) {
			S3Det hit = {Sd2rMul, Sd2rChannel->at(i), Sd2rEnergyRaw->at(i), Sd2rEnergy->at(i)};
			Sd2rDetector.push_back(hit);
		}
		
		//Sorting Yu
		if(!YuDetector.empty()) std::sort(YuDetector.begin(), YuDetector.end(), sortByEnergyYY1());
		
		//Sorting Sd1r
		if(!Sd1rDetector.empty()) std::sort(Sd1rDetector.begin(), Sd1rDetector.end(), sortByEnergyS3());
		
		//Sorting Sd2r
		if(!Sd2rDetector.empty()) std::sort(Sd2rDetector.begin(), Sd2rDetector.end(), sortByEnergyS3());
		
		if(Sd1rDetector.empty() || Sd2rDetector.empty()) continue;
		if(YuDetector.empty()) continue;
        if(ICEnergyRaw < 1500 || ICEnergyRaw > 2200) continue;
		if(!pidcut->IsInside(Sd2rDetector[0].energy,Sd1rDetector[0].energy)) continue;
		
		
		
		if(YuDetector.size()>1) {
			hYuEnM->Fill(YuDetector[0].energy,YuDetector[1].energy);
			hYuEnMICPID->Fill(YuDetector[0].energy,YuDetector[1].energy);
		}
		
		
		//Sd1r vs Sd2r calibrated
		if((Sd1rDetector.size()>0 && Sd2rDetector.size()>0) && Sd1rMul==1 && Sd2rMul==1){
            hSd1rSd2r->Fill(Sd2rEnergyRaw->at(0),Sd1rEnergyRaw->at(0)); //Non-calibrated Sd1r vs Sd2r energy - PID plot
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
			
			hCSd1rEn->Fill(Sd1rDetector[0].energy);
            hCSd2rEn->Fill(Sd2rDetector[0].energy);
            
            if(Sd1sEnergy->size()>0 && Sd1sEnergy->at(0)>0){
                hCSd1sEn->Fill(Sd1sEnergy->at(0));
            }
            
            if(Sd2sEnergy->size()>0 && Sd2sEnergy->at(0)>0){
                hCSd2sEn->Fill(Sd2sEnergy->at(0));
            }
            
            if (Sd1rDetector[0].channel<3) {
                hCSd1r_0_1_2_En->Fill(Sd1rDetector[0].energy);
            }
            
            if (Sd2rDetector[0].channel<3) {
                hCSd2r_0_1_2_En->Fill(Sd2rDetector[0].energy);
            }
			
		} //end of Sd1r vs Sd2r calibrated
		
		hICSd1rSd2r->Fill(Sd1rDetector[0].energy+Sd2rDetector[0].energy,ICEnergyRaw);
		
		
		//fill in the Yu detector
		if (!YuDetector.empty() && YuDetector[0].energy>1){
			//calculate the angles for the Yu detector
			yuM = Ydr0 + ( YuDetector[0].ring*YChWidth )+YChMiddle;
			TVector3 YuposM ( 0,yuM,Yuz); //shifting the detecor
			YuposM.SetPhi ( TMath::Pi() /2-YuDetector[0].sector *TMath::Pi() /4 ); //Pi/2 because the center of the sector is at 90degrees, Pi/4 is because there is 8 sectors so it is 2Pi/8
			//YuthetaM = beam.Angle ( YuposM ) *TMath::RadToDeg();
			//YuthetaM=ran->Uniform(Yubins[15-YuDetector[0].ring],Yubins[15-YuDetector[0].ring+1]);
            YuthetaM=DetAngle[YuDetector[0].sector][YuDetector[0].ring];
            
            //Randomization
            /*double phiRand = ran->Uniform(phi[YuDetector[0].sector]-21,phi[YuDetector[0].sector]+21); //Each detector covers 42 degreees
            double rRand = ran->Uniform(r[YuDetector[0].ring]-4.9375/2, r[YuDetector[0].ring]+4.9375/2); //Each ring has a width of 4.9375 mm
            double x_yu = rRand*sin(phiRand);
            double y_yu = rRand*cos(phiRand);
            double x_target = 2;
            double y_target = -2.6;
            TVector3 vec1(0, 0, 0-Yuz);
            TVector3 vec2(x_yu - x_target, y_yu - y_target, Yuz);
            double angle = vec1.Angle(vec2);
            YuthetaM=angle*180./M_PI;*/
                
			YuAngle->push_back(YuthetaM);
            YuEnergyShift=YuDetector[0].energy-shift;
							
			//Adding the energy lost by the protons through the target and the deadlayers of the Yu detector (11/15/2018)  
				
			// loss1=(-b[YuDetector[0].ring]-TMath::Sqrt(pow(b[YuDetector[0].ring],2)-4*(a[YuDetector[0].ring]-YuDetector[0].energy)*c[YuDetector[0].ring]))/(2*(a[YuDetector[0].ring]-YuDetector[0].energy));// Boron dead layer
            // loss2=(-e[YuDetector[0].ring]-TMath::Sqrt((e[YuDetector[0].ring]*e[YuDetector[0].ring])-4*(d[YuDetector[0].ring]-(YuDetector[0].energy+loss1))*f[YuDetector[0].ring]))/(2*(d[YuDetector[0].ring]-(YuDetector[0].energy+loss1)));// Aluminium dead layer
			// loss3=(-h[YuDetector[0].ring]-TMath::Sqrt((h[YuDetector[0].ring]*h[YuDetector[0].ring])-4*(g[YuDetector[0].ring]-(YuDetector[0].energy+loss1+loss2))*k[YuDetector[0].ring]))/(2*(g[YuDetector[0].ring]-(YuDetector[0].energy+loss1+loss2)));// Target
			// loss=loss1+loss2+loss3; //
				
			//if(YuDetector[0].ring == 0) {
				//double newEnergy = protonEL->AddBack(YuEnergy->at(0) + loss1 + loss2, 30./1000./cos(YuAngleMap[YuDetector[0].ring]*M_PI/180.));
				//std::cout << YuEnergy->at(0) + loss1 + loss2 << '\t' << YuEnergy->at(0) + loss1 + loss2 + loss3 << '\t' << newEnergy << std::endl;
			//}
            double DThickness = ran->Uniform(0,49.99);//60.78992
            double YuEnergyLoss = protonELB->AddBack(YuEnergyShift, 5e-5/fabs(cos(YuthetaM*M_PI/180.)));
            YuEnergyLoss = protonELA->AddBack(YuEnergyLoss, 1e-4/fabs(cos(YuthetaM*M_PI/180.)));
			YuEnergyLoss = protonELD->AddBack(YuEnergyLoss, (TThickness/2)/1000./fabs(cos(YuthetaM*M_PI/180.)));//
			//YuEnergyLoss=YuDetector[0].energy+loss;
			//cout << ev_num << "	" << YuDetector[0].ring << "		" << YuEnergy->at ( 0 ) << loss1 << loss2 << loss3 << loss << endl;
			hYuEn->Fill(YuEnergyLoss);
			//hYuEnT->Fill(YuDetector[0].energy+loss3);
			hYuAn->Fill(YuthetaM,YuEnergyLoss);
			//hYuAnT->Fill(YuthetaM,YuDetector[0].energy+loss3);
			hYuAnIC->Fill(YuthetaM,YuEnergyLoss);
			//hYuAnICT->Fill(YuthetaM,YuDetector[0].energy+loss3);
			hYuEnIC->Fill(YuEnergyLoss);
			hYuEnPID->Fill(YuEnergyLoss);
			hYuAnPID->Fill(YuthetaM,YuEnergyLoss);
			//hYuAnPIDT->Fill(YuthetaM,YuDetector[0].energy+loss3);
			hYuAnPIDSec[YuDetector[0].sector]->Fill(YuthetaM,YuEnergyLoss);
			hYuAnPIDRing[YuDetector[0].ring]->Fill(YuDetector[0].sector,YuEnergyLoss);
			Qval = ( 1+mejec/mrecoil ) * ( YuEnergyLoss) - ( 1 - mbeam/mrecoil ) * ( kBeam ) - 2 * TMath::Sqrt ( mbeam*mejec* ( YuEnergyLoss) * ( kBeam ) ) / mrecoil * TMath::Cos ( YuthetaM * TMath::Pi() / 180. );
			hQval->Fill ( Qval );
			hQvalAn->Fill(YuthetaM,Qval);
			//QvalT = ( 1+mejec/mrecoil ) * ( YuDetector[0].energy+loss3) - ( 1 - mbeam/mrecoil ) * ( kBeam ) - 2 * TMath::Sqrt ( mbeam*mejec* ( YuDetector[0].energy+loss3) * ( kBeam ) ) / mrecoil * TMath::Cos ( YuthetaM * TMath::Pi() / 180. );
            //std::cout << YuDetector[0].ring << '\t' << YuthetaM << std::endl;
			//hQvalT->Fill ( QvalT);
			//hQvalAnT->Fill(YuthetaM,QvalT);
						
			//Four sectors; left and right
			if(YuDetector[0].channel>=64 && YuDetector[0].channel<128){
				hYuAnPIDLeft->Fill(YuthetaM,YuEnergyLoss);
				hQvalLeft->Fill ( Qval );
			}
			else if(YuDetector[0].channel<64){
				hYuAnPIDRight->Fill(YuthetaM,YuEnergyLoss);
				hQvalRight->Fill ( Qval );
			}
						
			//Four sectors; top and bottom
			if(YuDetector[0].channel<32 || YuDetector[0].channel>=96){
				hYuAnPIDTop->Fill(YuthetaM,YuEnergyLoss);
				hQvalTop->Fill ( Qval );
			}
			else if(YuDetector[0].channel>=32 && YuDetector[0].channel<96){
				hYuAnPIDBottom->Fill(YuthetaM,YuEnergyLoss);
				hQvalBottom->Fill ( Qval );
			}
						
			//Adjacent sectors; first configuration
			if(YuDetector[0].channel>=112 || YuDetector[0].channel<16){
				hYuAnPID315_45->Fill(YuthetaM,YuEnergyLoss);
				hQval315_45->Fill ( Qval );
			}
			else if(YuDetector[0].channel>=16 && YuDetector[0].channel<48){
				hYuAnPID45_135->Fill(YuthetaM,YuEnergyLoss);
				hQval45_135->Fill ( Qval );
			}
			else if(YuDetector[0].channel>=48 && YuDetector[0].channel<80){
				hYuAnPID135_225->Fill(YuthetaM,YuEnergyLoss);
				hQval135_225->Fill ( Qval );
			}
			else if(YuDetector[0].channel>=80 && YuDetector[0].channel<112){
				hYuAnPID225_315->Fill(YuthetaM,YuEnergyLoss);
				hQval225_315->Fill ( Qval );
			}
						
			//Adjacent sectors; second configuration
			if(YuDetector[0].channel<32){
				hYuAnPID0_90->Fill(YuthetaM,YuEnergyLoss);
				hQval0_90->Fill ( Qval );
			}
			else if(YuDetector[0].channel>=32 && YuDetector[0].channel<64){
				hYuAnPID90_180->Fill(YuthetaM,YuEnergyLoss);
				hQval90_180->Fill ( Qval );
			}
			else if(YuDetector[0].channel>=64 && YuDetector[0].channel<96){
				hYuAnPID180_270->Fill(YuthetaM,YuEnergyLoss);
				hQval180_270->Fill ( Qval );
			}
			else if(YuDetector[0].channel>=96){
				hYuAnPID270_360->Fill(YuthetaM,YuEnergyLoss);
				hQval270_360->Fill ( Qval );
			}
			
			if(YuDetector[0].ring<4){
                //hYuAnPIDRs0_3->Fill(YuthetaM,YuEnergyLoss);
                hQvalRs0_3->Fill ( Qval );
            }else if(YuDetector[0].ring>=4 && YuDetector[0].ring<8){
                //hYuAnPIDRs0_3->Fill(YuthetaM,YuEnergyLoss);
                hQvalRs4_7->Fill ( Qval );
            }else if(YuDetector[0].ring>=8 && YuDetector[0].ring<12){
                hQvalRs8_11->Fill ( Qval );
            }else{
                hQvalRs12_15->Fill ( Qval );
            }
			
        
            hQvalR[YuDetector[0].ring]->Fill ( Qval );
            hQvalS[YuDetector[0].sector]->Fill ( Qval );

                
					
			//Rings in groups of 8
			if(YuDetector[0].ring<8){
                hYuAnPIDRs0_7->Fill(YuthetaM,YuEnergyLoss);
                hQvalRs0_7->Fill ( Qval );
                for(int i=0; i<8; i++){
                    if(YuDetector[0].ring != i){
                        hYuAnPIDR0_7[i]->Fill(YuthetaM,YuEnergyLoss);
                        hQvalR0_7[i]->Fill ( Qval );
                    }
                }
            }else if(YuDetector[0].ring>=8){
                hYuAnPIDRs8_15->Fill(YuthetaM,YuEnergyLoss);
                hQvalRs8_15->Fill ( Qval );
                for(int i=8; i<16; i++){
                    if(YuDetector[0].ring != i){
                        hYuAnPIDR8_15[i-8]->Fill(YuthetaM,YuEnergyLoss);
                        hQvalR8_15[i-8]->Fill ( Qval );
                    }
                }
			}
			
            if(YuDetector[0].ring>=4 && YuDetector[0].ring<12){
                //hYuAnPIDRs4_12->Fill(YuthetaM,YuEnergyLoss);
                hQvalRs4_12->Fill ( Qval );
            }else if(YuDetector[0].ring<4 || YuDetector[0].ring>=12){
               //hYuAnPIDRs12_4->Fill(YuthetaM,YuEnergyLoss);
                hQvalRs12_4->Fill ( Qval );
            }

			// Without the energy loss through the target and the dead layers
			hYuAn1->Fill(YuthetaM,YuDetector[0].energy); //all, no cuts
			hYuAnIC1->Fill(YuthetaM,YuEnergy->at (0));
			hYuEnIC1->Fill(YuEnergy->at (0));
			hYuEn1->Fill(YuDetector[0].energy);
			hYuMul->Fill(YuMul);
	
			hYuAnPID1->Fill ( YuthetaM,YuDetector[0].energy );
			Qval1 = ( 1+mejec/mrecoil ) * ( YuDetector[0].energy ) - ( 1 - mbeam/mrecoil ) * ( kBeam ) - 2 * TMath::Sqrt ( mbeam*mejec* ( YuDetector[0].energy ) * ( kBeam ) ) / mrecoil * TMath::Cos ( YuthetaM * TMath::Pi() / 180. );
			hQval1->Fill ( Qval1 );
		} 
		
	} //end of the main while loop
    
  
	YuAngle->clear();
    
	f_out->cd();
    
	//write the histograms in the output root file
    hSd1rSd2r->Write();
	hCSd1rSd2r->Write();
	hCSd1rSd2rIC->Write(); //calibrated ring vs ring with IC cut
	hCSd1rSd2rICCut->Write(); //calibrated ring vs ring with IC cut and diagonal taken into account
	hCSd1rSd2rYdIC->Write();
	hCSd1rSd2rYdICCut->Write();
    hCSd1rEn->Write();
    hCSd2rEn->Write();
    
    hCSd1sEn->Write();
    hCSd2sEn->Write();
    
    hCSd1r_0_1_2_En->Write();
    hCSd2r_0_1_2_En->Write();
    hICSd1rSd2r->Write();
    
	//Yu
	hYuAn->Write();
	hYuAnT->Write();
	hYuAn1->Write();
	hYuAnIC->Write();
	hYuAnICT->Write();
	hYuAnIC1->Write();
	hYuAnPID->Write();
	hYuAnPIDT->Write();
	hYuAnPID1->Write();
	hYuAnPIDLeft->Write();
	hYuAnPIDRight->Write();
	hYuAnPIDTop->Write();
	hYuAnPIDBottom->Write();
	hYuAnPID315_45->Write();
	hYuAnPID45_135->Write();
	hYuAnPID135_225->Write();
	hYuAnPID225_315->Write();
	hYuAnPID0_90->Write();
	hYuAnPID90_180->Write();
	hYuAnPID180_270->Write();
	hYuAnPID270_360->Write();
	hYuAnPIDRs0_7->Write();
	//hYuAnPIDR4_7->Write();
	//hYuAnPIDR8_11->Write();
	hYuAnPIDRs8_15->Write();
	hYuEn->Write();
	hYuEnT->Write();
	hYuEnIC->Write();
	hYuEnIC1->Write();
	hYuMul->Write();
	hYuEnPID->Write();
	hYuEn1->Write();
	hQval->Write();
	hQvalT->Write();
	hQval1->Write();
	hQvalLeft->Write();
	hQvalRight->Write();
	hQvalTop->Write();
	hQvalBottom->Write();
	hQval315_45->Write();
	hQval45_135->Write();
	hQval135_225->Write();
	hQval225_315->Write();
	hQval0_90->Write();
	hQval90_180->Write();
	hQval180_270->Write();
	hQval270_360->Write();
	hQvalRs0_7->Write();
    
    hQvalRs0_3->Write();
	hQvalRs4_7->Write();
	hQvalRs8_11->Write();
    hQvalRs12_15->Write();
    
    
	hQvalRs8_15->Write();
    hQvalRs4_12->Write();
    hQvalRs12_4->Write();
	hQvalAn->Write();
	hQvalAnT->Write();
	hYuEnM->Write();
	hYuEnMICPID->Write();
	for(int i=0; i<8; i++){hYuAnPIDSec[i]->Write();}
	for(int i=0; i<16; i++){hYuAnPIDRing[i]->Write();}
	for(int i=0; i<8; i++){
        hYuAnPIDR0_7[i]->Write();
        hYuAnPIDR8_15[i]->Write();
        hQvalR0_7[i]->Write();
        hQvalR8_15[i]->Write();
    }
    for(int i=0; i<16; i++){hQvalR[i]->Write();}
    for(int i=0; i<8; i++){hQvalS[i]->Write();}
	
	f_out->Close();
	
	return -Yuz;
	
} //end of the program
