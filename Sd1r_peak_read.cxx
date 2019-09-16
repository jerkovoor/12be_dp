//12C in-beam calibration of Sd detectors

using namespace std;

#include <TFile.h>
#include "TString.h"
#include "string"
#include "TTree.h"
#include "TBranch.h"
#include "iostream"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include <vector>
#include <TObject.h>
#include <TClass.h>
#include "TSpectrum.h"
#include <algorithm>
#include "TF1.h"
#include "fstream"

double Sd1r_peak_read()
{
  
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
    vector<int>* ICChannel= new vector<int>;
    vector<double>* ICEnergyRaw= new vector<double>;

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
 
    //Sd2s detector
    int Sd2sMul;
    vector<int>* Sd2sChannel= new vector<int>; //Sector Channels
    vector<double>* Sd2sEnergyRaw= new vector<double>;

    //Sur
    int SurMul;
    vector<int>* SurChannel= new vector<int>; //Ring Channels
    vector<double>* SurEnergyRaw= new vector<double>;

    //Sus
    int SusMul;
    vector<int>* SusChannel= new vector<int>; //Sector Channels
    vector<double>* SusEnergyRaw= new vector<double>;
    
  
    TChain* chain = new TChain ( "AutoTree" );
    
    //Alpha data
    chain->Add ( "/home/jerome/12Be_exp/Processed_files/rootindpndt_4972.root" ); //change the path to the files and the file name accordingly
    chain->Add ( "/home/jerome/12Be_exp/Processed_files/rootindpndt_4973.root" );
    chain->Add ( "/home/jerome/12Be_exp/Processed_files/rootindpndt_4974.root" );
    chain->Add ( "/home/jerome/12Be_exp/Processed_files/rootindpndt_4975.root" );
    chain->Add ( "/home/jerome/12Be_exp/Processed_files/rootindpndt_4977.root" );
    
    /*
    //Carbon data with target
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/target_thickness/12C/decode5001.root" ); //change the path to the files and the file name accordingly
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/target_thickness/12C/decode5002.root" );
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/target_thickness/12C/decode5003.root" );
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/target_thickness/12C/decode5004.root" );
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/target_thickness/12C/decode5006.root" );
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/target_thickness/12C/decode5007.root" );
	chain->Add ( "/home/jerome/12Be_exp/Processed_files/target_thickness/12C/decode5009.root" );
    */
    
    
    //Be Data
    /*
    
    for(int run_num=5021;run_num<5113;run_num++){
		if(run_num==5026||run_num==5040||run_num==5043||run_num==5046||run_num==5047||run_num==5059||run_num==5062||run_num==5063||run_num==5093||run_num==5099||run_num==5101||run_num==5093||run_num==5099||run_num==5101){
			continue;
		}else{

			string f_name=Form("/home/jerome/12Be_exp/Raw_root/12Be/tree%i.root",run_num);// new decode

			chain->Add(f_name.c_str());
		}
	}
    
	for(int run_num=5115;run_num<5133;run_num++){

		string f_name=Form("/home/jerome/12Be_exp/Raw_root/12Be/tree%i.root",run_num);// new decode

		chain->Add(f_name.c_str());
	}
	
	
	for(int run_num=5180;run_num<5223;run_num++){
		if(run_num==5187||run_num==5200){
			continue;
		}else{

			string f_name=Form("/home/jerome/12Be_exp/Raw_root/12Be/tree%i.root",run_num);// new decode

			chain->Add(f_name.c_str());
		}
	}
	
	*/
    
    TFile *f_out = new TFile ( "/home/jerome/12Be_exp/Analysis/Sd_Calibration/Sd1_AlphaPeaks.root","RECREATE" );
   
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

    chain->SetBranchAddress ( "Sd2sMul",&Sd2sMul );
    chain->SetBranchAddress ( "Sd2sChannel",&Sd2sChannel );
    chain->SetBranchAddress ( "Sd2sEnergyRaw",&Sd2sEnergyRaw );
 
    chain->SetBranchAddress ( "SurMul",&SurMul );
    chain->SetBranchAddress ( "SurChannel",&SurChannel );
    chain->SetBranchAddress ( "SurEnergyRaw",&SurEnergyRaw );
  
    chain->SetBranchAddress ( "SusMul",&SusMul );
    chain->SetBranchAddress ( "SusChannel",&SusChannel );
    chain->SetBranchAddress ( "SusEnergyRaw",&SusEnergyRaw );
    
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
    
    for ( int i=0; i<128; i++ )
    {
        string Yd_name = Form ( "Yd_Ch%i",i );
        hYd[i] =  new TH1D ( Yd_name.c_str(),"Yd_Ch",1024,0,4096 );

        string Yu_name = Form ( "Yu_Ch%i",i );
        hYu[i] =  new TH1D ( Yu_name.c_str(),"Yu_Ch",1024,0,4096); 

    }
    
    TH1D *hYuSectorRing[8][16];
    
    for(int i=0; i<8; i++)
    {
     for(int j=0; j<16; j++)
     {
            string nameSR = Form("Yu_Sec%i_Ring%i", i,j);
            hYuSectorRing[i][j] = new TH1D ( nameSR.c_str(),"Yu", 350,2400,3800 );
        }
    }
    
    for(int i=0; i<16; i++)
    {
     string CsI1_name = Form("CsI1_%i",i); 
     hCsI1[i] = new TH1D (CsI1_name.c_str(),"CsI1",1024,0,4096);
     
     string CsI2_name = Form("CsI2_%i",i); 
     hCsI2[i] = new TH1D (CsI2_name.c_str(),"CsI2",1024,0,4096);
      
    }
    
    for(int i=0; i<24; i++)
    {
        string Sd1r_name = Form ( "Sd1r_Ch%i",i );
        hSd1r[i] =  new TH1D ( Sd1r_name.c_str(),"Sd1r_Ch",500,0,2000 );

        string Sd2r_name = Form ( "Sd2r_Ch%i",i );
        hSd2r[i] =  new TH1D ( Sd2r_name.c_str(),"Sd2r_Ch",100,100,200 );

        string Sur_name = Form ( "Sur_Ch%i",i );
        hSur[i] =  new TH1D ( Sur_name.c_str(),"Sur_Ch",1024,0,4096 );
    }
    
    for(int i=0; i<32; i++)
    {
        string Sd1s_name = Form ( "Sd1s_Ch%i",i );
        hSd1s[i] =  new TH1D ( Sd1s_name.c_str(),"Sd1s_Ch",500,0,2000 );

        string Sd2s_name = Form ( "Sd2s_Ch%i",i );
        hSd2s[i] =  new TH1D ( Sd2s_name.c_str(),"Sd2s_Ch",100,100,200 );

        string Sus_name = Form ( "Sus_Ch%i",i );
        hSus[i] =  new TH1D ( Sus_name.c_str(),"Sus_Ch",1024,0,4096 );
    }
    
    
    
    int ev = chain->GetEntries();
    cout << "Total number of events =" << ev << endl;
    int ev_num=0;
    
    //Event by event starts
    
    while ( ev_num<=ev ){
        if ( ev_num%10000==0 ) cout << "Current event = " << ev_num << "\r"<< flush;
        chain->GetEntry ( ev_num );

        if ( YdEnergyRaw->size() >0 && YdMul==1 ){if ( YdEnergyRaw->at(0)>0 ){hYd[YdChannel->at(0)]->Fill ( YdEnergyRaw->at(0) );}}
        //cout << "I am at Yd" << endl;
        if(YuChannel->size()>0 && YuEnergyRaw->size()>0 && YuChannel->size()==YuEnergyRaw->size()){
            for(unsigned int i=0; i<YuEnergyRaw->size(); i++){ //cout << "I am inside for loop" << endl;
                if(YuEnergyRaw->at(i)>0){//cout << "Yu if condition" << endl;
                    hYu[YuChannel->at ( i )]->Fill ( YuEnergyRaw->at ( i ) );
                    hYuSectorRing[YuSector->at(i)][YuRing->at(i)]->Fill(YuEnergyRaw->at ( i ) ); 
                }
            }
        }
      
       
    
        if ( CsI1EnergyRaw->size() >0 && CsI1Mul==1 ){if ( CsI1EnergyRaw->at ( 0 ) >0 ){hCsI1[CsI1Channel->at ( 0 )]->Fill ( CsI1EnergyRaw->at ( 0 ) );}}
        if ( CsI2EnergyRaw->size() >0 && CsI2Mul==1 ){if ( CsI2EnergyRaw->at ( 0 ) >0 ){hCsI2[CsI2Channel->at ( 0 )]->Fill ( CsI2EnergyRaw->at ( 0 ) );}}
        if ( Sd1rEnergyRaw->size() >0 && Sd1rMul==1 ){if ( Sd1rEnergyRaw->at ( 0 ) >500 ){hSd1r[Sd1rChannel->at ( 0 )]->Fill ( Sd1rEnergyRaw->at ( 0 ) );}}
        if ( Sd1sEnergyRaw->size() >0 && Sd1sMul==1 ){if ( Sd1sEnergyRaw->at ( 0 ) >500 ){hSd1s[Sd1sChannel->at ( 0 )]->Fill ( Sd1sEnergyRaw->at ( 0 ) );}}
        if ( Sd2rEnergyRaw->size() >0 && Sd2rMul==1 ){if ( Sd2rEnergyRaw->at ( 0 ) >0 ){hSd2r[Sd2rChannel->at ( 0 )]->Fill ( Sd2rEnergyRaw->at ( 0 ) );}}
        if ( Sd2sEnergyRaw->size() >0 && Sd2sMul==1 ){if ( Sd2sEnergyRaw->at ( 0 ) >0 ){hSd2s[Sd2sChannel->at ( 0 )]->Fill ( Sd2sEnergyRaw->at ( 0 ) );}}
        if ( SurEnergyRaw->size() >0 && SurMul==1 ){if ( SurEnergyRaw->at ( 0 ) >0 ){hSur[SurChannel->at ( 0 )]->Fill ( SurEnergyRaw->at ( 0 ) );}}
        if ( SusEnergyRaw->size() >0 && SusMul==1 ){if ( SusEnergyRaw->at ( 0 ) >0 ){hSus[SusChannel->at ( 0 )]->Fill ( SusEnergyRaw->at ( 0 ) );}}

        ev_num++;
    }//end of the while loop
        
    //Event by event ends    
    f_out->cd();

    for ( int i=0; i<24; i++ ){hSd1r[i]->Write();}
    for ( int i=0; i<32; i++ ){hSd1s[i]->Write();}
    
    f_out->Close();

    return 0;
}
