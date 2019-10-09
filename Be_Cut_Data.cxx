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


void Be_Cut_Data(){

    

    TChain *chain = new TChain ( "AutoTree" );

    
    //Be data with target
		
	//First Half	
	
	for(int run_num=5021;run_num<5096;run_num++){
		if(run_num==5026||run_num==5040||run_num==5043||run_num==5046||run_num==5047||run_num==5059||run_num==5062||run_num==5063||run_num==5093||run_num==5099||run_num==5101){
			continue;
		}else{
			//string f_name=Form("/home/jerome/12Be_exp/Analysis/Be/AlphaOnly/decode_Yupedestal_%i.root",run_num);
			//string f_name=Form("/home/jerome/12Be_exp/Analysis/TRIUMF_DL/decodeBe_Yupedestal_TDL_%i.root",run_num);
			string f_name=Form("/home/jerome/12Be_exp/Analysis/Be_newdecode/decodeBe_Yupedestal_%i.root",run_num);// new decode
			//string f_name=Form("/home/jerome/12Be_exp/Analysis/TargetThickness/FromBeBeam/decodeNewBe_TDL_YuPedestal%i.root",run_num);// new decode
			chain->Add(f_name.c_str());
		}
	}
	
	//Second Half
	
    for(int run_num=5096;run_num<5113;run_num++){
        if(run_num==5093||run_num==5099||run_num==5101){
			continue;
		}else{
            //string f_name=Form("/home/jerome/12Be_exp/Analysis/Be/AlphaOnly/decode_Yupedestal_%i.root",run_num);
            //string f_name=Form("/home/jerome/12Be_exp/Analysis/TRIUMF_DL/decodeBe_Yupedestal_TDL_%i.root",run_num);
            string f_name=Form("/home/jerome/12Be_exp/Analysis/Be_newdecode/decodeBe_Yupedestal_%i.root",run_num);// new decode
            //string f_name=Form("/home/jerome/12Be_exp/Analysis/TargetThickness/FromBeBeam/decodeNewBe_TDL_YuPedestal%i.root",run_num);// new decode
            chain->Add(f_name.c_str());
        }
	}
    
	for(int run_num=5115;run_num<5133;run_num++){
		//string f_name=Form("/home/jerome/12Be_exp/Analysis/Be/AlphaOnly/decode_Yupedestal_%i.root",run_num);
		//string f_name=Form("/home/jerome/12Be_exp/Analysis/TRIUMF_DL/decodeBe_Yupedestal_TDL_%i.root",run_num);
		string f_name=Form("/home/jerome/12Be_exp/Analysis/Be_newdecode/decodeBe_Yupedestal_%i.root",run_num);// new decode
        //string f_name=Form("/home/jerome/12Be_exp/Analysis/TargetThickness/FromBeBeam/decodeNewBe_TDL_YuPedestal%i.root",run_num);// new decode
		chain->Add(f_name.c_str());
	}
	
	
	for(int run_num=5180;run_num<5223;run_num++){
		if(run_num==5187||run_num==5200){
			continue;
		}else{
			//string f_name=Form("/home/jerome/12Be_exp/Analysis/Be/AlphaOnly/decode_Yupedestal_%i.root",run_num);
			//string f_name=Form("/home/jerome/12Be_exp/Analysis/TRIUMF_DL/decodeBe_Yupedestal_TDL_%i.root",run_num);
			string f_name=Form("/home/jerome/12Be_exp/Analysis/Be_newdecode/decodeBe_Yupedestal_%i.root",run_num);// new decode
            //string f_name=Form("/home/jerome/12Be_exp/Analysis/TargetThickness/FromBeBeam/decodeNewBe_TDL_YuPedestal%i.root",run_num);// new decode
			chain->Add(f_name.c_str());
		}
	}

    //define the input variables
	//YYD detector
	int YdMul_in;
	vector<int>* YdChannel_in = new vector<int>();
	vector<double>* YdEnergyRaw_in= new vector<double>();
	vector<double>* YdEnergy_in= new vector<double>();
	vector<int>* YdRing_in= new vector<int>();
	vector<int>* YdSector_in= new vector<int>();
	
	//YYU detector
	int YuMul_in;
	vector<int>* YuChannel_in= new vector<int>();
	vector<double>* YuEnergyRaw_in= new vector<double>();
	vector<double>* YuEnergy_in= new vector<double>();
	vector<int>* YuRing_in= new vector<int>();
	vector<int>* YuSector_in= new vector<int>();
	
	//Purpose of two readouts is to increase the range of operation (multiply it by 2 in total)
	//CsI 1 detector
	int CsI1Mul_in;
	vector<int>* CsI1Channel_in= new vector<int>(); //Channel corresponds to one christal - 1 means that the gain is set to one value
	vector<double>* CsI1EnergyRaw_in= new vector<double>();
	
	//CsI 2 detector
	int CsI2Mul_in;
	vector<int>* CsI2Channel_in= new vector<int>();//Channel corresponds to one christal - 2 means that the gain is set to a different value than 1
	vector<double>* CsI2EnergyRaw_in= new vector<double>();
	
	//IC chamber
	int ICChannel_in;
	double ICEnergyRaw_in;
	
	//Sdd1 detector
	int Sd1rMul_in;
	vector<int>* Sd1rChannel_in= new vector<int>(); //Ring Channels
	vector<double>* Sd1rEnergyRaw_in= new vector<double>();
	vector<double>* Sd1rEnergy_in= new vector<double>();
 
	//Sd1s detector
	int Sd1sMul_in;
	vector<int>* Sd1sChannel_in = new vector<int>(); //Sector Channels
	vector<double>* Sd1sEnergyRaw_in= new vector<double>();
	vector<double>* Sd1sEnergy_in= new vector<double>();
	
	//Sd2r detector
	int Sd2rMul_in;
	vector<int>* Sd2rChannel_in= new vector<int>(); //Ring Channels
	vector<double>* Sd2rEnergyRaw_in= new vector<double>();
	vector<double>* Sd2rEnergy_in= new vector<double>();
	
	//Sd2s detector
	int Sd2sMul_in;
	vector<int>* Sd2sChannel_in= new vector<int>(); //Sector Channels
	vector<double>* Sd2sEnergyRaw_in= new vector<double>();
	vector<double>* Sd2sEnergy_in= new vector<double>();
	
	//Sur
	int SurMul_in;
	vector<int>* SurChannel_in= new vector<int>(); //Ring Channels
	vector<double>* SurEnergyRaw_in= new vector<double>();
	
	//Sus
	int SusMul_in;
	vector<int>* SusChannel_in= new vector<int>(); //Sector Channels
	vector<double>* SusEnergyRaw_in= new vector<double>();
	
	//reading the input tree
	chain->SetBranchAddress ( "YdMul",&YdMul_in );
	chain->SetBranchAddress ( "YdChannel",&YdChannel_in );
	chain->SetBranchAddress ( "YdEnergyRaw",&YdEnergyRaw_in );
	chain->SetBranchAddress ( "YdEnergy",&YdEnergy_in );
	chain->SetBranchAddress ( "YdRing",&YdRing_in );
	chain->SetBranchAddress ( "YdSector",&YdSector_in );
	
	chain->SetBranchAddress ( "YuMul",&YuMul_in );
	chain->SetBranchAddress ( "YuChannel",&YuChannel_in );
	chain->SetBranchAddress ( "YuEnergyRaw",&YuEnergyRaw_in );
	chain->SetBranchAddress ( "YuEnergy",&YuEnergy_in );
	chain->SetBranchAddress ( "YuRing",&YuRing_in );
	chain->SetBranchAddress ( "YuSector",&YuSector_in );
	
	chain->SetBranchAddress ( "CsI1Mul",&CsI1Mul_in );
	chain->SetBranchAddress ( "CsI1Channel",&CsI1Channel_in );
	chain->SetBranchAddress ( "CsI1EnergyRaw",&CsI1EnergyRaw_in );
	
	chain->SetBranchAddress ( "CsI2Mul",&CsI2Mul_in );
	chain->SetBranchAddress ( "CsI2Channel",&CsI2Channel_in );
	chain->SetBranchAddress ( "CsI2EnergyRaw",&CsI2EnergyRaw_in );
	
	chain->SetBranchAddress ( "ICChannel",&ICChannel_in );
	chain->SetBranchAddress ( "ICEnergyRaw",&ICEnergyRaw_in );
	
	chain->SetBranchAddress ( "Sd1rMul",&Sd1rMul_in );
	chain->SetBranchAddress ( "Sd1rChannel",&Sd1rChannel_in );
	chain->SetBranchAddress ( "Sd1rEnergyRaw",&Sd1rEnergyRaw_in );
	chain->SetBranchAddress ( "Sd1rEnergy",&Sd1rEnergy_in );
	
	chain->SetBranchAddress ( "Sd1sMul",&Sd1sMul_in );
	chain->SetBranchAddress ( "Sd1sChannel",&Sd1sChannel_in );
	chain->SetBranchAddress ( "Sd1sEnergyRaw",&Sd1sEnergyRaw_in );
	chain->SetBranchAddress ( "Sd1sEnergy",&Sd1sEnergy_in );
	
	chain->SetBranchAddress ( "Sd2rMul",&Sd2rMul_in );
	chain->SetBranchAddress ( "Sd2rChannel",&Sd2rChannel_in );
	chain->SetBranchAddress ( "Sd2rEnergyRaw",&Sd2rEnergyRaw_in );
	chain->SetBranchAddress ( "Sd2rEnergy",&Sd2rEnergy_in );
	
	chain->SetBranchAddress ( "Sd2sMul",&Sd2sMul_in );
	chain->SetBranchAddress ( "Sd2sChannel",&Sd2sChannel_in );
	chain->SetBranchAddress ( "Sd2sEnergyRaw",&Sd2sEnergyRaw_in );
	chain->SetBranchAddress ( "Sd2sEnergy",&Sd2sEnergy_in );
	
	chain->SetBranchAddress ( "SurMul",&SurMul_in );
	chain->SetBranchAddress ( "SurChannel",&SurChannel_in );
	chain->SetBranchAddress ( "SurEnergyRaw",&SurEnergyRaw_in );
	
	chain->SetBranchAddress ( "SusMul",&SusMul_in );
	chain->SetBranchAddress ( "SusChannel",&SusChannel_in );
	chain->SetBranchAddress ( "SusEnergyRaw",&SusEnergyRaw_in );
    
    
    
    //output variables
	//YYD detector
	int YdMulO;
	vector<int> YdChannelO;
	vector<double> YdEnergyRawO;
	vector<double> YdEnergyO;
	vector<int> YdRingO;
	vector<int> YdSectorO;

	//YYU detector
	int YuMulO;
	vector<int> YuChannelO;
	vector<double> YuEnergyRawO;
	vector<double> YuEnergyO;
	vector<int> YuRingO;
	vector<int> YuSectorO;
	
	//Purpose of two readouts is to increase the range of operation (multiply it by 2 in total)
	//CsI 1 detector
	int CsI1MulO;
	vector<int> CsI1ChannelO; //Channel corresponds to one christal - 1 means that the gain is set to one value
	vector<double> CsI1EnergyRawO;

	//CsI 2 detector
	int CsI2MulO;
	vector<int> CsI2ChannelO;//Channel corresponds to one christal - 2 means that the gain is set to a different value than 1
	vector<double> CsI2EnergyRawO;

	//IC chamber
	int ICChannelO;
	double ICEnergyRawO;
	
	//Sdd1 detector
	int Sd1rMulO;
	vector<int> Sd1rChannelO; //Ring Channels
	vector<double> Sd1rEnergyRawO;
	vector<double> Sd1rEnergyO;

	//Sd1s detector
	int Sd1sMulO;
	vector<int> Sd1sChannelO ; //Sector Channels
	vector<double> Sd1sEnergyRawO;
	vector<double> Sd1sEnergyO;
	
	//Sd2r detector
	int Sd2rMulO;
	vector<int> Sd2rChannelO; //Ring Channels
	vector<double> Sd2rEnergyRawO;
	vector<double> Sd2rEnergyO;
	
	//Sd2s detector
	int Sd2sMulO;
	vector<int> Sd2sChannelO; //Sector Channels
	vector<double> Sd2sEnergyRawO;
	vector<double> Sd2sEnergyO;

	//Sur
	int SurMulO;
	vector<int> SurChannelO; //Ring Channels
	vector<double> SurEnergyRawO;

	//Sus
	int SusMulO;
	vector<int> SusChannelO; //Sector Channels
	vector<double> SusEnergyRawO;

	//open the output file
	
	//string fOut_name = Form ( "/home/jerome/12Be_exp/Analysis/Su_calibration/decode_%i_%d.root",run_num,index); //Open the output file; put the correct path and file name
	//string fOut_name = Form ( "/home/jerome/12Be_exp/Analysis/TargetThickness/FromBeBeam/decodeNewBe_TDL_YuPedestal%i.root",run_num ); //Open the output file; put the correct path and file name
    //string fOut_name = Form ( "/home/jerome/12Be_exp/Analysis/TargetThickness/FromBeBeam/Be_noTarget/decodeNewBe_noTarget_TDL_YuPedestal%i.root",run_num ); //Be no target
    //std::cout << fOut_name << std::endl;
	TFile* f_out = new TFile ( "/home/jerome/12Be_exp/Analysis/BeCutData.root","RECREATE" );
	TTree* tr_out = new TTree ( "AutoTree","AutoTree" ); //create the new root tree
	//Branches created in the new file; the leafs are removed from the three as well as arrays of vectors or whatever that was
	tr_out->Branch ( "YdMulO",&YdMulO );
	tr_out->Branch ( "YdChannelO",&YdChannelO );
	tr_out->Branch ( "YdEnergyRawO",&YdEnergyRawO );
	tr_out->Branch ( "YdEnergyO",&YdEnergyO );
	tr_out->Branch ( "YdRingO",&YdRingO );
	tr_out->Branch ( "YdSectorO",&YdSectorO );

	tr_out->Branch ( "YuMulO",&YuMulO );
	tr_out->Branch ( "YuChannelO",&YuChannelO );
	tr_out->Branch ( "YuEnergyRawO",&YuEnergyRawO );
	tr_out->Branch ( "YuEnergyO",&YuEnergyO );
	tr_out->Branch ( "YuRingO",&YuRingO );
	tr_out->Branch ( "YuSectorO",&YuSectorO );
	
	tr_out->Branch ( "CsI1MulO",&CsI1MulO );
	tr_out->Branch ( "CsI1ChannelO",&CsI1ChannelO );
	tr_out->Branch ( "CsI1EnergyRawO",&CsI1EnergyRawO );

	tr_out->Branch ( "CsI2MulO",&CsI2MulO );
	tr_out->Branch ( "CsI2ChannelO",&CsI2ChannelO );
	tr_out->Branch ( "CsI2EnergyRawO",&CsI2EnergyRawO );
	
	tr_out->Branch ( "ICChannelO",&ICChannelO );
	tr_out->Branch ( "ICEnergyRawO",&ICEnergyRawO );

	tr_out->Branch ( "Sd1rMulO",&Sd1rMulO );
	tr_out->Branch ( "Sd1rChannelO",&Sd1rChannelO );
	tr_out->Branch ( "Sd1rEnergyRawO",&Sd1rEnergyRawO );
	tr_out->Branch ( "Sd1rEnergyO",&Sd1rEnergyO );

	tr_out->Branch ( "Sd1sMulO",&Sd1sMulO );
	tr_out->Branch ( "Sd1sChannelO",&Sd1sChannelO );
	tr_out->Branch ( "Sd1sEnergyRawO",&Sd1sEnergyRawO );
	tr_out->Branch ( "Sd1sEnergyO",&Sd1sEnergyO );

	tr_out->Branch ( "Sd2rMulO",&Sd2rMulO );
	tr_out->Branch ( "Sd2rChannelO",&Sd2rChannelO );
	tr_out->Branch ( "Sd2rEnergyRawO",&Sd2rEnergyRawO );
	tr_out->Branch ( "Sd2rEnergyO",&Sd2rEnergyO );

	tr_out->Branch ( "Sd2sMulO",&Sd2sMulO );
	tr_out->Branch ( "Sd2sChannelO",&Sd2sChannelO );
	tr_out->Branch ( "Sd2sEnergyRawO",&Sd2sEnergyRawO );
	tr_out->Branch ( "Sd2sEnergyO",&Sd2sEnergyO );
	
	tr_out->Branch ( "SurMulO",&SurMulO );
	tr_out->Branch ( "SurChannelO",&SurChannelO );
	tr_out->Branch ( "SurEnergyRawO",&SurEnergyRawO );
	
	tr_out->Branch ( "SusMulO",&SusMulO );
	tr_out->Branch ( "SusChannelO",&SusChannelO );
	tr_out->Branch ( "SusEnergyRawO",&SusEnergyRawO );
    
    TFile *f_cut = TFile::Open("/home/jerome/12Be_exp/Analysis/Be/cuts/BeCutIC.root"); //PID cut for Be
	TCutG *pidcut = (TCutG*) f_cut->Get("BeCut_CalSd1rSd2rIC");
    
    long long ev = chain->GetEntries(); //get the total number of events in the tree
	cout << "Total number of events =" << ev << endl;
	long long ev_num=0; 
    
    for(long long ev_num = 0; ev_num < ev; ev_num++) {
    //for(long long ev_num = 0; ev_num < 200000; ev_num++) {
    //while( ev_num<=ev ){ //entering the main while loop and looping over all the events in the Iris tree
		//if ( ev_num%1000==0 ) cout << "Current event = " << ev_num << "\r"<< endl; // flush;
        if ( ev_num%100000==0 ) cout << "Current event = " << ev_num << "\r"<< flush;
		chain->GetEntry ( ev_num ); //pulling up the current event from the Iris tree
        
        //clear the vectors so you do not multiply the data
        YdMulO = 0;
		YdChannelO.clear();
		YdEnergyRawO.clear();
		YdEnergyO.clear();
		YdRingO.clear();
		YdSectorO.clear();
		
        YuMulO = 0;
		YuChannelO.clear();
		YuEnergyRawO.clear();
		YuEnergyO.clear();
		YuRingO.clear();
		YuSectorO.clear();
		
        CsI1MulO = 0;
		CsI1ChannelO.clear();
		CsI1EnergyRawO.clear();
        CsI2MulO = 0;
		CsI2ChannelO.clear();
		CsI2EnergyRawO.clear();
		
        Sd1rMulO = 0;
		Sd1rChannelO.clear();
		Sd1rEnergyRawO.clear();
		Sd1rEnergyO.clear();
        Sd1sMulO = 0;
		Sd1sChannelO.clear();
		Sd1sEnergyRawO.clear();
		Sd1sEnergyO.clear();

        Sd2rMulO = 0;
		Sd2rChannelO.clear();
		Sd2rEnergyRawO.clear();
		Sd2rEnergyO.clear();
        Sd2sMulO = 0;
		Sd2sChannelO.clear();
		Sd2sEnergyRawO.clear();
		Sd2sEnergyO.clear();

        SurMulO = 0;
		SurChannelO.clear();
		SurEnergyRawO.clear();
        SusMulO = 0;
		SusChannelO.clear();
		SusEnergyRawO.clear();

		ICChannelO = -10;
		ICEnergyRawO = -10;
        
        
        //Defining a structure YuDetector
		vector<YuDet> YuDetector;
		for(size_t i = 0; i < YuEnergy_in->size(); i++) {
			if(YuChannel_in->at(i) == 82 || YuChannel_in->at(i) == 96 || YuChannel_in->at(i) == 106 || YuChannel_in->at(i) == 111) continue;
			YuDet hit = {YuMul_in, YuChannel_in->at(i), YuRing_in->at(i), YuSector_in->at(i), YuEnergyRaw_in->at(i), YuEnergy_in->at(i)};
			YuDetector.push_back(hit);
		}
		
		//Defining a structure Sd1rDetector
		vector<S3Det> Sd1rDetector;
		for(size_t i = 0; i < Sd1rEnergy_in->size(); i++) {
			S3Det hit = {Sd1rMul_in, Sd1rChannel_in->at(i), Sd1rEnergyRaw_in->at(i), Sd1rEnergy_in->at(i)};
			Sd1rDetector.push_back(hit);
		}
		
		//Defining a structure Sd2rDetector
		vector<S3Det> Sd2rDetector;
		for(size_t i = 0; i < Sd2rEnergy_in->size(); i++) {
			S3Det hit = {Sd2rMul_in, Sd2rChannel_in->at(i), Sd2rEnergyRaw_in->at(i), Sd2rEnergy_in->at(i)};
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
        if(ICEnergyRaw_in < 620 || ICEnergyRaw_in > 1100) continue;
		if(!pidcut->IsInside(Sd2rDetector[0].energy,Sd1rDetector[0].energy)) continue;
		
		
		//////////////////////////////////////////////////////////
		//						Yd DETECTOR						//
		//////////////////////////////////////////////////////////
        //cout << YdMul_in << endl;
		//treating the Yd detector
        
		if(YdChannel_in->size()>0){
			YdMulO=YdMul_in; //set the multiplicity to 0
			//cout << YdMulO << endl;
			//  cout << " Entered the Yd loop " << endl;
			//cout << "Size " << det->TYdADC.size() << endl;
            for ( unsigned int i=0; i < YdChannel_in->size(); i++ ){ //loop over the size of the YdADC vector, this is always = 128
                YdChannelO.push_back (YdChannel_in->at(i));
                YdEnergyRawO.push_back (YdEnergyRaw_in->at(i));
                YdEnergyO.push_back (YdEnergy_in->at(i));
                YdRingO.push_back (YdRing_in->at(i));
                YdSectorO.push_back (YdSector_in->at(i));
            }
		}
		
		//treating the CsI detector, the same way as it was done for the Yd detector
		//The signals from th CsI detector were split in two and send to two different gain settings, so the channels for the CsI are multiplied
		//////////////////////////////////////////////////////////
		//						CsI1 DETECTOR					//
		//////////////////////////////////////////////////////////

		if ( CsI1Channel_in->size() >0 ){
			//E layerDet
			CsI1MulO=CsI1Mul_in;
			for ( unsigned int i=0; i < CsI1Channel_in->size(); i++ ){
                CsI1ChannelO.push_back (CsI1Channel_in->at(i));
                CsI1EnergyRawO.push_back (CsI1EnergyRaw_in->at(i));
			}
				
		}

		//////////////////////////////////////////////////////////
		//						CsI2 DETECTOR					//
		//////////////////////////////////////////////////////////
		
		if ( CsI2Channel_in->size() >0 ){
			//E layerDet
			CsI2MulO=CsI2Mul_in;
			for ( unsigned int i=0; i < CsI2Channel_in->size(); i++ ){
                CsI2ChannelO.push_back (CsI2Channel_in->at(i));
                CsI2EnergyRawO.push_back (CsI2EnergyRaw_in->at(i));
			}
				
		}
		

		//////////////////////////////////////////////////////////
		//						YU DETECTOR						//
		//Treat the Yu detector, same as all the previous ones	//
		//////////////////////////////////////////////////////////
		
		if(YuChannel_in->size()>0){
			YuMulO=YuMul_in; //set the multiplicity to 0
			//  cout << " Entered the Yd loop " << endl;
			//cout << "Size " << det->TYdADC.size() << endl;
            for ( unsigned int i=0; i < YuMul_in; i++ ){ //loop over the size of the YdADC vector, this is always = 128
                YuChannelO.push_back (YuChannel_in->at(i));
                YuEnergyRawO.push_back (YuEnergyRaw_in->at(i));
                YuEnergyO.push_back (YuEnergy_in->at(i));
                YuRingO.push_back (YuRing_in->at(i));
                YuSectorO.push_back (YuSector_in->at(i));
            }
		}
		
			
		//////////////////////////////////////////////////////////
		//						IC CHAMBER						//
		//////////////////////////////////////////////////////////
		//The IC chamber has only one channel, it is plugged in channel 15
        ICChannelO = ICChannel_in;
        ICEnergyRawO = ICEnergyRaw_in;
		
		
		//////////////////////////////////////////////////////////
		//						Sd1r DETECTOR					//
		//////////////////////////////////////////////////////////
		//Treating the Sd1r part of the Sd1 detector. There are 2 Sd detectors. Sd1 is the one making the dE layer in the dE-E telechope. It has 24 rings and 32 sectors. The rings are oriented toward the beam. 
		//The data is gathered in a way that the sectors and rings are registered separately which gives the rings as Sd1r. The Sd2 detector is treated the same as Sd1. 
        
		if (Sd1rChannel_in->size() >0 ){
			// cout << " size of the vector " << det->TSd1rADC.size() << endl;
			Sd1rMulO=Sd1rMul_in;
			for ( unsigned int i=0; i < Sd1rChannel_in->size(); i++ ){
				Sd1rChannelO.push_back (Sd1rChannel_in->at(i));
                Sd1rEnergyRawO.push_back (Sd1rEnergyRaw_in->at(i));
                Sd1rEnergyO.push_back (Sd1rEnergy_in->at(i));
			}
		}
		
		//////////////////////////////////////////////////////////
		//						Sd1s DETECTOR					//
		//////////////////////////////////////////////////////////
		//the sector side of the Sd3 detector not facing the beam - dE layer
		
		if (Sd1sChannel_in->size() >0 ){
			// cout << " size of the vector " << det->TSd1rADC.size() << endl;
			Sd1sMulO=Sd1sMul_in;
			for ( unsigned int i=0; i < Sd1sChannel_in->size(); i++ ){
				Sd1sChannelO.push_back (Sd1sChannel_in->at(i));
                Sd1sEnergyRawO.push_back (Sd1sEnergyRaw_in->at(i));
                Sd1sEnergyO.push_back (Sd1sEnergy_in->at(i));
			}
		}

		//////////////////////////////////////////////////////////
		//						Sd2r DETECTOR					//
		//////////////////////////////////////////////////////////
		//The Sd2 detector faces the beem with its sector side; This is the treatment of the ring side of the detector
		
		if (Sd2rChannel_in->size() >0 ){
			// cout << " size of the vector " << det->TSd1rADC.size() << endl;
			Sd2rMulO=Sd2rMul_in;
			for ( unsigned int i=0; i < Sd2rChannel_in->size(); i++ ){
				Sd2rChannelO.push_back (Sd2rChannel_in->at(i));
                Sd2rEnergyRawO.push_back (Sd2rEnergyRaw_in->at(i));
                Sd2rEnergyO.push_back (Sd2rEnergy_in->at(i));
			}
		}
		
		//////////////////////////////////////////////////////////
		//						Sd2s DETECTOR					//
		//////////////////////////////////////////////////////////
		//The sector side of the S3 detector facing the beam - E layer
		if (Sd2sChannel_in->size() >0 ){
			// cout << " size of the vector " << det->TSd1rADC.size() << endl;
			Sd2sMulO=Sd2sMul_in;
			for ( unsigned int i=0; i < Sd2sChannel_in->size(); i++ ){
				Sd2sChannelO.push_back (Sd2sChannel_in->at(i));
                Sd2sEnergyRawO.push_back (Sd2sEnergyRaw_in->at(i));
                Sd2sEnergyO.push_back (Sd2sEnergy_in->at(i));
			}
		}
		
		//////////////////////////////////////////////////////////
		//						Sur DETECTOR	`				//
		//////////////////////////////////////////////////////////
		//another S3 detector positioned upstream in front of the Yu detector
		//The ring side of this detector is facing the target and its NOT facing the beam
		if (SurChannel_in->size() >0 ){
			// cout << " size of the vector " << det->TSd1rADC.size() << endl;
			SurMulO=SurMul_in;
			for ( unsigned int i=0; i < SurChannel_in->size(); i++ ){
				SurChannelO.push_back (SurChannel_in->at(i));
                SurEnergyRawO.push_back (SurEnergyRaw_in->at(i));
			}
		}
		
		//////////////////////////////////////////////////////////
		//						Sus DETECTOR	`				//
		//////////////////////////////////////////////////////////
		//sector side of the upstream S3 detector is facing the beam
		if (SusChannel_in->size() >0 ){
			// cout << " size of the vector " << det->TSd1rADC.size() << endl;
			SusMulO=SusMul_in;
			for ( unsigned int i=0; i < SusChannel_in->size(); i++ ){
				SusChannelO.push_back (SusChannel_in->at(i));
                SusEnergyRawO.push_back (SusEnergyRaw_in->at(i));
			}
		}


		tr_out->Fill(); //fill out the tree

		//ev_num++; //looping over the while loop
	} //end of the main while loop

    
    f_out->cd(); //access the output file
	tr_out->Write(); //write the output file
	f_out->Close(); //close the output file

}
