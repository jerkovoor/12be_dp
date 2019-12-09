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


void SingleRunAnalysis(){
	//open the output file
	//TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/Be/CarbonGain/Be_pedestal_TRIUMF_DL_CarbonGain_QvalMinGS_Random_NewDecode_Cut2_Yu.root","RECREATE"); //change the path and name accordingly
    //TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n2.root","RECREATE"); //change the path and name accordingly
    //TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2_CutData_TargetFront.root","RECREATE");
    
    TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/TargetThickness/FromBeBeam/Be_target/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2_target_thickness_Sd1rAlphaCal_Sd2rNewInBeamCal_SingleRun.root","RECREATE");
    
    //TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/TargetThickness/FromBeBeam/Be_noTarget/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n1_2_noTarget_target_thickness_Sd1rAlphaCal_Sd2rNewInBeamCal.root","RECREATE");
    
	//Open the input files
	
    
    //First Half	
	TH1D *hCSd1rEn[140];//(5223-5021)-13-(5179-5133+1)-2=140 runs
    TH1D *hCSd1r_0_1_2_En[140];
    TH1D *hCSd1rEn_YuGated[140];
    TH1D *hCSd2rEn[140];
    TH1D *hCSd2r_0_1_2_En[140];
    TH1D *hCSd2rEn_YuGated[140];
    int num=0;
	for(int run_num=5021;run_num<5223;run_num++){//5113
		if(run_num==5026||run_num==5040||run_num==5043||run_num==5046||run_num==5047||run_num==5059||run_num==5062||run_num==5063||run_num==5093||run_num==5099||run_num==5101||run_num==5113||run_num==5114||(5133<=run_num && run_num<5180)||run_num==5187||run_num==5200){
			continue;
		}else{
			//string f_name=Form("/home/jerome/12Be_exp/Analysis/Be/AlphaOnly/decode_Yupedestal_%i.root",run_num);
			//string f_name=Form("/home/jerome/12Be_exp/Analysis/TRIUMF_DL/decodeBe_Yupedestal_TDL_%i.root",run_num);
			//string f_name=Form("/home/jerome/12Be_exp/Analysis/Be_newdecode/decodeBe_Yupedestal_%i.root",run_num);// new decode
			string f_name=Form("/home/jerome/12Be_exp/Analysis/TargetThickness/FromBeBeam/Be_target/decodeNewBe_Sd1rAlphaCal_Sd2rNewInBeamCal_TDL_YuPedestal%i.root",run_num);// new decode
            TFile* f_in = new TFile (f_name.c_str());
			TTree* tr_out = (TTree*)f_in->Get("AutoTree");
            
            //YYU detector
            int YuMul;
            vector<int>* YuChannel= new vector<int>();
            vector<double>* YuEnergyRaw= new vector<double>();
            vector<double>* YuEnergy= new vector<double>();
            vector<int>* YuRing= new vector<int>();
            vector<int>* YuSector= new vector<int>();
            
            //Sd1r detector
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
            
            tr_out->SetBranchAddress ( "YuMul",&YuMul );
            tr_out->SetBranchAddress ( "YuChannel",&YuChannel );
            tr_out->SetBranchAddress ( "YuEnergyRaw",&YuEnergyRaw );
            tr_out->SetBranchAddress ( "YuEnergy",&YuEnergy );
            tr_out->SetBranchAddress ( "YuRing",&YuRing );
            tr_out->SetBranchAddress ( "YuSector",&YuSector );
            
            tr_out->SetBranchAddress ( "Sd1rMul",&Sd1rMul );
            tr_out->SetBranchAddress ( "Sd1rChannel",&Sd1rChannel );
            tr_out->SetBranchAddress ( "Sd1rEnergyRaw",&Sd1rEnergyRaw );
            tr_out->SetBranchAddress ( "Sd1rEnergy",&Sd1rEnergy );
            
            tr_out->SetBranchAddress ( "Sd1sMul",&Sd1sMul );
            tr_out->SetBranchAddress ( "Sd1sChannel",&Sd1sChannel );
            tr_out->SetBranchAddress ( "Sd1sEnergyRaw",&Sd1sEnergyRaw );
            tr_out->SetBranchAddress ( "Sd1sEnergy",&Sd1sEnergy );
            
            tr_out->SetBranchAddress ( "Sd2rMul",&Sd2rMul );
            tr_out->SetBranchAddress ( "Sd2rChannel",&Sd2rChannel );
            tr_out->SetBranchAddress ( "Sd2rEnergyRaw",&Sd2rEnergyRaw );
            tr_out->SetBranchAddress ( "Sd2rEnergy",&Sd2rEnergy );
            
            tr_out->SetBranchAddress ( "Sd2sMul",&Sd2sMul );
            tr_out->SetBranchAddress ( "Sd2sChannel",&Sd2sChannel );
            tr_out->SetBranchAddress ( "Sd2sEnergyRaw",&Sd2sEnergyRaw );
            tr_out->SetBranchAddress ( "Sd2sEnergy",&Sd2sEnergy );
            
            
            cout << run_num << "\t" << num << endl;
            hCSd1rEn[num] = new TH1D(Form("hCSd1rEn_%i",num),Form("Calibrated SD1r Energy_%i",run_num),500,0,50);
            hCSd1rEn_YuGated[num] = new TH1D(Form("hCSd1rEn_YuGated_%i",num),Form("Calibrated SD1r Energy_%i anti-gated on Yu",run_num),500,0,50);
            
            
            hCSd1r_0_1_2_En[num] = new TH1D(Form("hCSd1r_0_1_2_En_%i",num),Form("Calibrated SD1r Energy (Rings 0, 1, and 2)_%i",run_num),500,0,50);
            
            //TH1D *hCSd1sEn = new TH1D("hCSd1sEn","Calibrated SD1s Energy",500,0,50);
            
            hCSd2rEn[num] = new TH1D(Form("hCSd2rEn_%i",num),Form("Calibrated SD2r Energy_%i",run_num),500,0,120);
            hCSd2rEn_YuGated[num] = new TH1D(Form("hCSd2rEn_YuGated_%i",num),Form("Calibrated SD2r Energy_%i anti-gated on Yu",run_num),500,0,120);
            
            
            hCSd2r_0_1_2_En[num] = new TH1D(Form("hCSd2r_0_1_2_En_%i",num),Form("Calibrated SD2r Energy (Rings 0, 1, and 2)_%i",run_num),500,0,120);
            
            //TH1D *hCSd2sEn = new TH1D("hCSd2sEn","Calibrated SD2s Energy",500,0,120);
            //start reading the tree
            int ev = tr_out->GetEntries(); //get the total number of entries
            cout << "Total number of events =" << ev << endl;
            int ev_num=0;
            
            ////////////////////////////////
            // Event by event starts here //
            ////////////////////////////////
            
                    
            for(int ev_num = 0; ev_num < ev; ev_num++) {
                if ( ev_num%50000==0 ) cout << "Current event = " << ev_num << "\r"<< flush;
                tr_out->GetEntry ( ev_num ); //get the current entry
                
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
                //if(YuDetector.empty()) continue;
                //if(ICEnergyRaw < 620 || ICEnergyRaw > 1100) continue;
                //if(!pidcut->IsInside(Sd2rDetector[0].energy,Sd1rDetector[0].energy)) continue;

                hCSd1rEn[num]->Fill(Sd1rDetector[0].energy);
                hCSd2rEn[num]->Fill(Sd2rDetector[0].energy);
                
                if(YuDetector.empty()){//Anti-gating on Yu
                    hCSd1rEn_YuGated[num]->Fill(Sd1rDetector[0].energy);
                    hCSd2rEn_YuGated[num]->Fill(Sd2rDetector[0].energy);
                }
                
                if (Sd1rDetector[0].channel<3) {
                    hCSd1r_0_1_2_En[num]->Fill(Sd1rDetector[0].energy);
                }
                
                if (Sd2rDetector[0].channel<3) {
                    hCSd2r_0_1_2_En[num]->Fill(Sd2rDetector[0].energy);
                }
            } //end of the main while loop
		}
		num=num+1;
	}
	
	f_out->cd();
    
    
    for(int i=0; i<140; i++){
        hCSd1rEn[i]->Write();
        hCSd2rEn[i]->Write();
        hCSd1r_0_1_2_En[i]->Write();
        hCSd2r_0_1_2_En[i]->Write();
        hCSd1rEn_YuGated[i]->Write();
        hCSd2rEn_YuGated[i]->Write();
    }
    
    f_out->Close();
	
	
	
	
    /*
	for(int run_num=5115;run_num<5133;run_num++){
		//string f_name=Form("/home/jerome/12Be_exp/Analysis/Be/AlphaOnly/decode_Yupedestal_%i.root",run_num);
		//string f_name=Form("/home/jerome/12Be_exp/Analysis/TRIUMF_DL/decodeBe_Yupedestal_TDL_%i.root",run_num);
		//string f_name=Form("/home/jerome/12Be_exp/Analysis/Be_newdecode/decodeBe_Yupedestal_%i.root",run_num);// new decode
        string f_name=Form("/home/jerome/12Be_exp/Analysis/TargetThickness/FromBeBeam/Be_target/decodeNewBe_Sd1rAlphaCal_Sd2rNewInBeamCal_TDL_YuPedestal%i.root",run_num);// new decode
		tr_out->Add(f_name.c_str());
	}
	
	
	for(int run_num=5180;run_num<5223;run_num++){
		if(run_num==5187||run_num==5200){
			continue;
		}else{
			//string f_name=Form("/home/jerome/12Be_exp/Analysis/Be/AlphaOnly/decode_Yupedestal_%i.root",run_num);
			//string f_name=Form("/home/jerome/12Be_exp/Analysis/TRIUMF_DL/decodeBe_Yupedestal_TDL_%i.root",run_num);
			//string f_name=Form("/home/jerome/12Be_exp/Analysis/Be_newdecode/decodeBe_Yupedestal_%i.root",run_num);// new decode
            string f_name=Form("/home/jerome/12Be_exp/Analysis/TargetThickness/FromBeBeam/Be_target/decodeNewBe_Sd1rAlphaCal_Sd2rNewInBeamCal_TDL_YuPedestal%i.root",run_num);// new decode
			tr_out->Add(f_name.c_str());
		}
	}
*/

}
