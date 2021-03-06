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
	
struct sortByEnergyYY1 {
	inline bool operator() (const YuDet& EnYY1_1,
							const YuDet& EnYY1_2){
		return (EnYY1_1.energy>EnYY1_2.energy);
	}
};

struct sortByEnergyS3 {
	inline bool operator() (const Sd1rDet& EnS3_1,
							const Sd1rDet& EnS3_2){
		return (EnS3_1.energy>EnS3_2.energy);
	}
};


double AnalysisBeamOffset(){
	//open the output file
	TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/C_pedestal_TRIUMF_DL_BeamOffset_2.root","RECREATE"); //change the path and name accordingly
  
	//Open the input files
	TChain *chain = new TChain ( "AutoTree" ); //chain the desired input files
    
    EnergyLoss* protonELB = new EnergyLoss("Proton_Boron.dat");
    EnergyLoss* protonELA = new EnergyLoss("Proton_Aluminum.dat");
	EnergyLoss* protonELD = new EnergyLoss("Proton_DeuteriumTarget.dat");
    
    protonELB->AddBackHigh(15.);
    protonELA->AddBackHigh(15.);
    protonELD->AddBackHigh(15.);
	
	//Carbon data with target
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/decodeNew_TDL_5001.root" ); //change the path to the files and the file name accordingly
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/decodeNew_TDL_5002.root" );
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/decodeNew_TDL_5003.root" );
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/decodeNew_TDL_5004.root" );
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/decodeNew_TDL_5006.root" );
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/decodeNew_TDL_5007.root" );
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/decodeNew_TDL_5009.root" );
	
	/*chain->Add ( "/home/jerome/12Be_exp/Analysis/with_pedestal/decode_Yupedestal_5001.root" ); //change the path to the files and the file name accordingly
	chain->Add ( "/home/jerome/12Be_exp/Analysis/with_pedestal/decode_Yupedestal_5002.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/with_pedestal/decode_Yupedestal_5003.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/with_pedestal/decode_Yupedestal_5004.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/with_pedestal/decode_Yupedestal_5006.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/with_pedestal/decode_Yupedestal_5007.root" );
	chain->Add ( "/home/jerome/12Be_exp/Analysis/with_pedestal/decode_Yupedestal_5009.root" );*/
	
	
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
	
	TFile *f_cut = TFile::Open("/home/jerome/12Be_exp/Analysis/CarbonGain/13CcutFull.root"); //Full PID cut for C
	TCutG *pidcut = (TCutG*) f_cut->Get("CcutCalSd1rSd2rFull"); // for C
	
	/*TFile *f_cut = TFile::Open("/home/jerome/12Be_exp/Analysis/Be/BeCut1.root"); //PID cut for Be
	TCutG *pidcut = (TCutG*) f_cut->Get("BeCut_CalSd1rSd2r"); */ // for Be
	
	//definition of histograms
	//calculate the variable bins for the Yu detector to plot the angles
	double Yubins[17]={0};
	
	for (int i=0; i<17; i++){
		Yubins[i]=180-TMath::RadToDeg()*TMath::ATan((50+(16-i)*4.94)/(0-Yuz));
		//cout << Yubins[i] << endl;
	}
	
	//TH2D *hYuAnPID = new TH2D("hYuAnPID","YuE vs Angle with an IC and PID gate",16,Yubins,6000,0,6); //Energy vs angle in Yu with a gate on IC and PID//240,0,2.4
	//TH1D *hQval = new TH1D("hQval","Q values",400,-10,10);
	
	//start reading the tree
	int ev = chain->GetEntries(); //get the total number of entries
	cout << "Total number of events =" << ev << endl;
	int ev_num=0;
	
	//variable calculation for the YY1 detectors used in angle calculations
	float YChWidth = ( Ydr1 - Ydr0 ) /16.;
	float YChMiddle = YChWidth/2;
	float yuM=0, ydM=0;
	double YuthetaM, YuthetaR, YdthetaM; //angle for Yu/Yd
	TVector3 beam(0,0,1); //beam vector
    
	//variable calculation for the S3 detectors used in angle calculations
	float SdChWidth = ( Sdr1 - Sdr0 ) /24.;
	float SdChMiddle = SdChWidth/2;
	float SdM=0;
	double SdthetaM; //angle for Sd1
	
	//variable definition for the Q value calculations
	float amu = 931.5; // atomic mass unit in MeV
	float massEjec = 938.28; //mass of the proton in MeV/c2 
	float kBeam = 110; //put the correct value; beam energy
	float mbeam = 12 * amu;  //mass of the beam (12Be or 12C) in MeV
	float mrecoil = 13 * amu;  //mass of the recoil (13Be or 13C) in MeV
	float mejec = 1 * amu; //mass of the proton
    

	double Qval; //Q value variable
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
	
	double A, B;

	double OffsetInitialx=-1;
    double OffsetFinalx=1;
    double OffsetInitialy=0;
    double OffsetFinaly=-1;
    double interval=0.2;
    int OffsetLengthx=1+(OffsetFinalx-OffsetInitialx)/interval;
    int OffsetLengthy=1+fabs(OffsetFinaly-OffsetInitialy)/interval;
    int Length=OffsetLengthx*OffsetLengthy;
    
    double x[OffsetLengthx], y[OffsetLengthy];
    
    for (int i=0;i<OffsetLengthx;i++){
        x[i]=OffsetInitialx+i*interval;
    }
    
    for (int i=0;i<OffsetLengthy;i++){
        y[i]=OffsetInitialy+i*interval;
    }
    
    double phi[8], r[16], R[8][16][OffsetLengthx][OffsetLengthy], DetAngle[8][16][OffsetLengthx][OffsetLengthy];
    phi[0]=22.5;
    for (int i=1;i<8;i++){
        phi[i]=phi[0]+45*i;
        //cout << phi[i] << endl;
    }
    
    for (int i=0;i<16;i++){
        r[i]=50+((129.-50.)/16)*(0.5+i);
        //cout << r[i] << endl;
    }
    
    //Offsetting the beam
        for (int xindex=0;xindex<OffsetLengthx;xindex++){
            for (int yindex=0;yindex<OffsetLengthy;yindex++){
                //Calculating the angles
                for (int i=0;i<8;i++){
                    for (int j=0;j<16;j++){
                        A=r[j]*TMath::Sin(phi[i]);
                        B=r[j]*TMath::Cos(phi[i]);
                        R[i][j][xindex][yindex]=TMath::Sqrt(pow((x[xindex]-A),2)+pow((y[yindex]-B),2));
                        DetAngle[i][j][xindex][yindex]=180-TMath::RadToDeg()*TMath::ATan(R[i][j][xindex][yindex]/85);
                        //cout << DetAngle[i][j][xindex][yindex] << endl;
                    }
                }
            }
        }

        

	TH2D *hYuAnPID[Length];
	TH1D *hQval[Length];
	for(int i=0;i<OffsetLengthx;i++){
        for(int j=0;j<OffsetLengthy;j++){
		//string namehYuAnPID = Form("hYuAnPID_%i",i);
		//hYuAnPID[i] = new TH2D(namehYuAnPID.c_str(),"YuE vs Angle with an IC and PID gate",16,Yubins,6000,0,6);//240,0,2.4
		string namehQval = Form("hQval_%i",i*OffsetLengthy+j);
		hQval[i*OffsetLengthy+j] = new TH1D(namehQval.c_str(),"Q values",160,-4,4);//240,0,2.4
        }
	}
	
	double YuEnergyLoss, YuEnergyShift;
	//TF1* fit_func1 = new TF1 ( "fit_func1","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",-10,10);
	int seed = 123;
	TRandom3* ran = new TRandom3(seed);

	double shift=0.154;//0.0924;
	
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
		
		if(Sd1rDetector.empty() || Sd2rDetector.empty()) continue;
		if(YuDetector.empty()) continue;
		if(!pidcut->IsInside(Sd2rDetector[0].energy,Sd1rDetector[0].energy) || ICEnergyRaw<1500 || ICEnergyRaw>2200) continue;
		
		//Sorting Yu
		if(!YuDetector.empty()){
			std::sort(YuDetector.begin(),YuDetector.end(),
			sortByEnergyYY1());
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
		
		
		//fill in the Yu detector
        if (YuDetector[0].energy>1){
            //calculate the angles for the Yu detector
            //yuM = Ydr0+ ( YuDetector[0].ring*YChWidth )+YChMiddle;
            //TVector3 YuposM ( 0,yuM,Yuz); //shifting the detecor
            //YuposM.SetPhi ( TMath::Pi() /2-YuDetector[0].sector *TMath::Pi() /4 ); //Pi/2 because the center of the sector is at 90degrees, Pi/4 is because there is 8 sectors so it is 2Pi/8
            //YuthetaM = beam.Angle ( YuposM ) *TMath::RadToDeg();
            //YuthetaM=ran->Uniform(Yubins[15-YuDetector[0].ring],Yubins[15-YuDetector[0].ring+1]);
        
            YuEnergyShift=YuDetector[0].energy-shift;
            //Adding the energy lost by the protons through the target and the deadlayers of the Yu detector (11/15/2018)  
            /*loss1=(-b[YuDetector[0].ring]-TMath::Sqrt(pow(b[YuDetector[0].ring],2)-4*(a[YuDetector[0].ring]-YuEnergyShift)*c[YuDetector[0].ring]))/(2*(a[YuDetector[0].ring]-YuEnergyShift));// Boron dead layer
            loss2=(-e[YuDetector[0].ring]-TMath::Sqrt((e[YuDetector[0].ring]*e[YuDetector[0].ring])-4*(d[YuDetector[0].ring]-(YuEnergyShift+loss1))*f[YuDetector[0].ring]))/(2*(d[YuDetector[0].ring]-(YuEnergyShift+loss1)));// Aluminium dead layer
            loss3=(-h[YuDetector[0].ring]-TMath::Sqrt((h[YuDetector[0].ring]*h[YuDetector[0].ring])-4*(g[YuDetector[0].ring]-(YuEnergyShift+loss1+loss2))*k[YuDetector[0].ring]))/(2*(g[YuDetector[0].ring]-(YuEnergyShift+loss1+loss2)));// Target
            loss=loss1+loss2+loss3; //
            YuEnergyLoss=YuEnergyShift+loss;*/
                

            //Offsetting the beam
            for (int xindex=0;xindex<OffsetLengthx;xindex++){
                for (int yindex=0;yindex<OffsetLengthy;yindex++){
                    //Calculating the angles
                    YuthetaM=DetAngle[YuDetector[0].sector][YuDetector[0].ring][xindex][yindex];
                    //cout << ev_num << "		" << YuDetector[0].ring << "		" << YuthetaM << "		" << YuthetaM <<endl;

                    //if(pidcut->IsInside(Sd2rDetector[0].energy,Sd1rDetector[0].energy) && ICEnergyRaw>620 && ICEnergyRaw<1100 ){
                    double YuEnergyLoss = protonELB->AddBack(YuEnergyShift, 5e-5/fabs(cos(YuthetaM*M_PI/180.)));
                    YuEnergyLoss = protonELA->AddBack(YuEnergyLoss, 1e-4/fabs(cos(YuthetaM*M_PI/180.)));
                    YuEnergyLoss = protonELD->AddBack(YuEnergyLoss, 30./1000./fabs(cos(YuthetaM*M_PI/180.)));
                    Qval = ( 1+mejec/mrecoil ) * ( YuEnergyLoss) - ( 1 - mbeam/mrecoil ) * ( kBeam ) - 2 * TMath::Sqrt ( mbeam*mejec* ( YuEnergyLoss) * ( kBeam ) ) / mrecoil * TMath::Cos ( YuthetaM * TMath::Pi() / 180. );
                    hQval[xindex*OffsetLengthy+yindex]->Fill ( Qval );
                    //hQvalAn->Fill(YuthetaM,Qval);
                }
            }
            YuAngle->push_back(YuthetaM);
        }
	} //end of the main for loop
    
  
	YuAngle->clear();
    
	f_out->cd();
    
	//write the histograms in the output root file
	
	for(int i=0; i<Length; i++){
		//hYuAnPID[i]->Write();
		hQval[i]->Write();
	}
	
	f_out->Close();
	
	return 0;
	
} //end of the program
