#define Reduce_cxx
#include "Reduce.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <string.h>
#include <vector>
#include <iostream>

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

int main() {
    TChain *chain = new TChain("AutoTree");
    
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
    
    Reduce t(chain);
    t.Loop();
    
    return 0;
}


void Reduce::Loop() {

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    
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

    Long64_t nbytes = 0, nb = 0;
    for(Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if(ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        
        if ( jentry%1000==0 ) cout << "Current event = " << jentry << "\r"<< endl; // flush;
        
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
        
        if(YuMul > 0) {
            cout << YuMul << '\t' << YuChannel->size() << '\t' << YuChannel[0][0] << endl;
        }
        
//         //Defining a structure YuDetector
// 		vector<YuDet> YuDetector;
// 		for(size_t i = 0; i < YuMul; i++) {
// 			if(YuChannel[i] == 82 || YuChannel[i] == 96 || YuChannel[i] == 106 || YuChannel[i] == 111) continue;
// 			YuDet hit = {YuMul, YuChannel[i], YuRing[i], YuSector[i], YuEnergyRaw[i], YuEnergy[i]};
// 			YuDetector.push_back(hit);
// 		}
// 		
// 		//Defining a structure Sd1rDetector
// 		vector<S3Det> Sd1rDetector;
// 		for(size_t i = 0; i < Sd1rMul; i++) {
// 			S3Det hit = {Sd1rMul, Sd1rChannel[i], Sd1rEnergyRaw[i], Sd1rEnergy[i]};
// 			Sd1rDetector.push_back(hit);
// 		}
// 		
// 		//Defining a structure Sd2rDetector
// 		vector<S3Det> Sd2rDetector;
// 		for(size_t i = 0; i < Sd2rMul; i++) {
// 			S3Det hit = {Sd2rMul, Sd2rChannel[i], Sd2rEnergyRaw[i], Sd2rEnergy[i]};
// 			Sd2rDetector.push_back(hit);
// 		}
//         
//         //Sorting Yu
// 		if(!YuDetector.empty()) std::sort(YuDetector.begin(), YuDetector.end(), sortByEnergyYY1());
// 		
// 		//Sorting Sd1r
// 		if(!Sd1rDetector.empty()) std::sort(Sd1rDetector.begin(), Sd1rDetector.end(), sortByEnergyS3());
// 		
// 		//Sorting Sd2r
// 		if(!Sd2rDetector.empty()) std::sort(Sd2rDetector.begin(), Sd2rDetector.end(), sortByEnergyS3());
//         
//         
//         if(Sd1rDetector.empty() || Sd2rDetector.empty()) continue;
// 		if(YuDetector.empty()) continue;
//         if(ICEnergyRaw < 620 || ICEnergyRaw > 1100) continue;
// 		if(!pidcut->IsInside(Sd2rDetector[0].energy,Sd1rDetector[0].energy)) continue;
        
        
    }
}
