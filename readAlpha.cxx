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

double readAlpha()
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

    // Yd calibration before the experiment
  /*  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/Calibration_data/Alpha4972.root" );
    chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/Calibration_data/Alpha4973.root" );
    chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/Calibration_data/Alpha4974.root" );
    chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/Calibration_data/Alpha4975.root" );
    chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/Calibration_data/Alpha4977.root" );
    chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/Calibration_data/Alpha4978.root" );*/
   // chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/Calibration_data/Alpha4980.root" );

  
   chain->Add ( "/home/jerome/12Be_exp/Processed_files/rootindpndt_5226_1.root" );
   chain->Add ( "/home/jerome/12Be_exp/Processed_files/rootindpndt_5226_2.root" );
   
   //chain->Add ( "/home/jerome/12Be_exp/Processed_files/pedestalindpndt_4815.root" );
   
   /*chain->Add ( "/home/jerome/12Be_exp/Processed_files/rootindpndt_4972.root" );
   chain->Add ( "/home/jerome/12Be_exp/Processed_files/rootindpndt_4973.root" );
   chain->Add ( "/home/jerome/12Be_exp/Processed_files/rootindpndt_4974.root" );
   chain->Add ( "/home/jerome/12Be_exp/Processed_files/rootindpndt_4975.root" );
   chain->Add ( "/home/jerome/12Be_exp/Processed_files/rootindpndt_4977.root" );*/
   
   
  
  //after unpacking for the 2nd time
  /*chain->Add("/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/AlphaTest5225_1.root");
  chain->Add("/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/AlphaTest5225_2.root");
  chain->Add("/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/AlphaTest5225_3.root");*/

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
    
    TFile *f_out = new TFile ( "HisAlpha_Yd_calib.root","RECREATE" );
    
    int ev = chain->GetEntries();
    cout << "Total number of events =" << ev << endl;
    int ev_num=0;
    
    while ( ev_num<=ev )
    {
        if ( ev_num%10000==0 ) cout << "Current event = " << ev_num << "\r"<< flush;
        chain->GetEntry ( ev_num );

       if ( YdEnergyRaw->size() >0 && YdMul==1 ){if ( YdEnergyRaw->at(0)>500 ){hYd[YdChannel->at(0)]->Fill ( YdEnergyRaw->at(0) );}}
       //cout << "I am at Yd" << endl;
       if(YuChannel->size()>0 && YuEnergyRaw->size()>0 && YuChannel->size()==YuEnergyRaw->size()){
	for(unsigned int i=0; i<YuEnergyRaw->size(); i++)
	{ //cout << "I am inside for loop" << endl;
	  if(YuEnergyRaw->at(i)>500)
	 {//cout << "Yu if condition" << endl;
	    hYu[YuChannel->at ( i )]->Fill ( YuEnergyRaw->at ( i ) );
	    hYuSectorRing[YuSector->at(i)][YuRing->at(i)]->Fill(YuEnergyRaw->at ( i ) ); 
	 }}}
      
       
    
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
    
    TF1* fit_func6 = new TF1 ( "fit_func6","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp(-pow((x-[10]),2)/(2*[2]*[2]))+[11]*TMath::Exp(-pow((x-[12]),2)/(2*[2]*[2]))",500,5000 ); //removed the liner part for comparison
    TF1* fit_func3 = new TF1 ( "fit_func3","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))",0,5000 );
    TF1* fit_func9 = new TF1 ( "fit_func9","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp(-pow((x-[10]),2)/(2*[2]*[2]))+[11]*TMath::Exp(-pow((x-[12]),2)/(2*[2]*[2]))+[13]*TMath::Exp(-pow((x-[14]),2)/(2*[2]*[2]))+[15]*TMath::Exp(-pow((x-[16]),2)/(2*[2]*[2]))+[17]*TMath::Exp(-pow((x-[18]),2)/(2*[2]*[2]))+[19]+[20]*x",2000,5000 );
    
    //Search for peaks in Yd strips
    int npeaks=3;
    TSpectrum* s = new TSpectrum ( 2*npeaks );
    int nfound=0;
    double* txp = s->GetPositionX();
    double xp[3];
    double a[128],b[128],c[128],ped[128];

    for ( int j=0; j<128; j++ )
    {
        nfound = s->Search ( hYd[j],2,"",0.1 );
       
        for ( int p=0; p<nfound; p++ ){xp[p] = txp[p];}
        
        /*ped[j] = min ( min ( xp[0], xp[1] ),min ( xp[2],xp[3] ) );
        c[j] = max ( max ( xp[0], xp[1] ),max ( xp[2],xp[3] ) );

       if((xp[0]!=ped[j] && xp[0]!=c[j]) && (xp[1]!=ped[j] && xp[1]!=c[j])){a[j]=min(xp[0],xp[1]);b[j]=max(xp[0],xp[1]);}
       else if ((xp[0]!=ped[j] && xp[0]!=c[j]) && (xp[2]!=ped[j] && xp[2]!=c[j])){a[j]=min(xp[0],xp[2]);b[j]=max(xp[0],xp[2]);}
       else if ((xp[0]!=ped[j] && xp[0]!=c[j]) && (xp[3]!=ped[j] && xp[3]!=c[j])){a[j]=min(xp[0],xp[3]);b[j]=max(xp[0],xp[3]);}
       else if ((xp[1]!=ped[j] && xp[1]!=c[j]) && (xp[2]!=ped[j] && xp[2]!=c[j])){a[j]=min(xp[1],xp[2]);b[j]=max(xp[1],xp[2]);}
       else if ((xp[1]!=ped[j] && xp[1]!=c[j]) && (xp[3]!=ped[j] && xp[3]!=c[j])){a[j]=min(xp[1],xp[3]);b[j]=max(xp[1],xp[3]);}
       else if ((xp[2]!=ped[j] && xp[2]!=c[j]) && (xp[3]!=ped[j] && xp[3]!=c[j])){a[j]=min(xp[2],xp[3]);b[j]=max(xp[2],xp[3]);}
       
      // cout << j << "pedestal " << ped[j] << " peak 1 " << a[j] << " peak 2 " << b[j] << " peak 3 " << c[j] << endl;
               
    }*/

   /* for ( int j=0; j<128; j++ )
    {
        nfound = s->Search ( hYu[j],2,"",0.1 );

        for ( int p=0; p<nfound; p++ ){ xp[p] = txp[p];}
                
       // ped[j] = min ( min ( xp[0], xp[1] ),min ( xp[2],xp[3] ) );
       // c[j] = max ( max ( xp[0], xp[1] ),max ( xp[2],xp[3] ) );

      // if((xp[0]!=ped[j] && xp[0]!=c[j]) && (xp[1]!=ped[j] && xp[1]!=c[j])){a[j]=min(xp[0],xp[1]);b[j]=max(xp[0],xp[1]);}
      // else if ((xp[0]!=ped[j] && xp[0]!=c[j]) && (xp[2]!=ped[j] && xp[2]!=c[j])){a[j]=min(xp[0],xp[2]);b[j]=max(xp[0],xp[2]);}
     //  else if ((xp[0]!=ped[j] && xp[0]!=c[j]) && (xp[3]!=ped[j] && xp[3]!=c[j])){a[j]=min(xp[0],xp[3]);b[j]=max(xp[0],xp[3]);}
     //  else if ((xp[1]!=ped[j] && xp[1]!=c[j]) && (xp[2]!=ped[j] && xp[2]!=c[j])){a[j]=min(xp[1],xp[2]);b[j]=max(xp[1],xp[2]);}
     //  else if ((xp[1]!=ped[j] && xp[1]!=c[j]) && (xp[3]!=ped[j] && xp[3]!=c[j])){a[j]=min(xp[1],xp[3]);b[j]=max(xp[1],xp[3]);}
     //  else if ((xp[2]!=ped[j] && xp[2]!=c[j]) && (xp[3]!=ped[j] && xp[3]!=c[j])){a[j]=min(xp[2],xp[3]);b[j]=max(xp[2],xp[3]);}*/
       
        a[j] = min ( min (xp[0],xp[1]), xp[2]);
        c[j] = max ( max (xp[0],xp[1]), xp[2]);
    
    if(xp[0]!=a[j] && xp[0]!=c[j]) {b[j] = xp[0];}
    else if(xp[1]!=a[j] && xp[1]!=c[j]) {b[j] = xp[1];}
    else if(xp[2]!=a[j] && xp[2]!=c[j]) {b[j] = xp[2];}
       
      // cout << j << "pedestal " << ped[j] << " peak 1 " << a[j] << " peak 2 " << b[j] << " peak 3 " << c[j] << endl;
      cout << j << " peak 1 " << a[j] << " peak 2 " << b[j] << " peak 3 " << c[j] << endl;
               
    }
  
    //open a txt file to store the positions of the peaks
    ofstream Yd_alpha;
    Yd_alpha.open ( "Yd_AlphaPeaks.txt" );
    Yd_alpha << " strip / 1 / 2 / 3  \n";

  
    ofstream Yd_alphaCounts;
    Yd_alphaCounts.open ( "YdCounts.txt" );
    Yd_alphaCounts << " strip / Counts 1 / 2 / 3  \n";
      
    for ( int i=0; i<128; i++ )
    {
        cout << i << " /a " << a[i] << " /b " << b[i] << " /c " << c[i] << endl;
        //if ( i<55 )
        //{
           /* fit_func9->SetParLimits ( 1,a[i]-100,a[i]-50 );
            fit_func9->SetParLimits ( 4,a[i]-30,a[i]-10 );
            fit_func9->SetParLimits ( 6,a[i]-10,a[i]+10 );
            fit_func9->SetParLimits ( 8,b[i]-100,b[i]-50 );
            fit_func9->SetParLimits ( 10,b[i]-30,b[i]-10 );
            fit_func9->SetParLimits ( 12,b[i]-10,b[i]+10 );
            fit_func9->SetParLimits ( 14,c[i]-100,c[i]-50 );
            fit_func9->SetParLimits ( 16,c[i]-30,c[i]-10 );
            fit_func9->SetParLimits ( 18,c[i]-10,c[i]+10 );
            fit_func9->SetParLimits ( 2,0,25 );

            hYd[i]->Fit ( "fit_func9","","",a[i]-100,c[i]+100 );

            //Yd_alpha << i << " " << ped[i] << " "  << fit_func9->GetParameter ( 6 ) << " " << fit_func9->GetParameter ( 12 ) << " " << fit_func9->GetParameter ( 18 ) << endl;
	    Yd_alpha << i << " "  << fit_func9->GetParameter ( 6 ) << " " << fit_func9->GetParameter ( 12 ) << " " << fit_func9->GetParameter ( 18 ) << endl;
            Yu_alphaCounts << i << " " << fit_func9->GetParameter ( 5 ) << " " <<  fit_func9->GetParameter ( 11 ) << " " <<  fit_func9->GetParameter ( 17 ) << endl;*/
        //}
        //else if ( i>55 )
       // {
            fit_func6->SetParLimits ( 1,a[i]-30,a[i]-10 );
            fit_func6->SetParLimits ( 4,a[i]-10,a[i]+10 );
            fit_func6->SetParLimits ( 6,b[i]-30,b[i]-10 );
            fit_func6->SetParLimits ( 8,b[i]-10,b[i]+10 );
            fit_func6->SetParLimits ( 10,c[i]-30,c[i]-10 );
            fit_func6->SetParLimits ( 12,c[i]-10,c[i]+10 );
            fit_func6->SetParLimits ( 2,5,15 );
            
            //if(i==44||i==46||i==60||i==106){hYd[i]->Fit ( "fit_func6","","",2500,3500 );}

            hYd[i]->Fit ( "fit_func6","","",a[i]-100,c[i]+100 );

           // Yd_alpha << i << " " << ped[i] << " "  << fit_func6->GetParameter ( 4 ) << " " << fit_func6->GetParameter ( 8 ) << " " << fit_func6->GetParameter ( 12 ) << endl;
	    Yd_alpha << i << " "  << fit_func6->GetParameter ( 4 ) << " " << fit_func6->GetParameter ( 8 ) << " " << fit_func6->GetParameter ( 12 ) << endl;
            Yd_alphaCounts << i << " " << fit_func6->GetParameter ( 3 ) << " " <<  fit_func6->GetParameter ( 7 ) << " " <<  fit_func6->GetParameter ( 11 ) << endl;
     //   }
           
    }
  //end of the fit on the Yd detector
  
  //find peaks and fir the Sd1r detector
  /*  double Sdrxp[3];
    double Sdra[24],Sdrb[24],Sdrc[24],Sdrped[24];

    for ( int j=0; j<24; j++ )
    {
        nfound = s->Search ( hSd2r[j],2,"",0.1 );
       
        for ( int p=0; p<nfound; p++ ){Sdrxp[p] = txp[p];}
        
       // Sdrped[j] = min ( min ( Sdrxp[0], Sdrxp[1] ),min ( Sdrxp[2],Sdrxp[3] ) );
       // Sdrc[j] = max ( max ( Sdrxp[0], Sdrxp[1] ),max ( Sdrxp[2],Sdrxp[3] ) );

      // if((Sdrxp[0]!=Sdrped[j] && Sdrxp[0]!=Sdrc[j]) && (Sdrxp[1]!=Sdrped[j] && Sdrxp[1]!=Sdrc[j])){Sdra[j]=min(Sdrxp[0],Sdrxp[1]);Sdrb[j]=max(Sdrxp[0],Sdrxp[1]);}
      // else if ((Sdrxp[0]!=Sdrped[j] && Sdrxp[0]!=Sdrc[j]) && (Sdrxp[2]!=Sdrped[j] && Sdrxp[2]!=Sdrc[j])){Sdra[j]=min(Sdrxp[0],Sdrxp[2]);Sdrb[j]=max(Sdrxp[0],Sdrxp[2]);}
      // else if ((Sdrxp[0]!=Sdrped[j] && Sdrxp[0]!=Sdrc[j]) && (Sdrxp[3]!=Sdrped[j] && Sdrxp[3]!=Sdrc[j])){Sdra[j]=min(Sdrxp[0],Sdrxp[3]);Sdrb[j]=max(Sdrxp[0],Sdrxp[3]);}
      // else if ((Sdrxp[1]!=Sdrped[j] && Sdrxp[1]!=Sdrc[j]) && (Sdrxp[2]!=Sdrped[j] && Sdrxp[2]!=Sdrc[j])){Sdra[j]=min(Sdrxp[1],Sdrxp[2]);Sdrb[j]=max(Sdrxp[1],Sdrxp[2]);}
      // else if ((Sdrxp[1]!=Sdrped[j] && Sdrxp[1]!=Sdrc[j]) && (Sdrxp[3]!=Sdrped[j] && Sdrxp[3]!=Sdrc[j])){Sdra[j]=min(Sdrxp[1],Sdrxp[3]);Sdrb[j]=max(Sdrxp[1],Sdrxp[3]);}
    //   else if ((Sdrxp[2]!=Sdrped[j] && Sdrxp[2]!=Sdrc[j]) && (Sdrxp[3]!=Sdrped[j] && Sdrxp[3]!=Sdrc[j])){Sdra[j]=min(Sdrxp[2],Sdrxp[3]);Sdrb[j]=max(Sdrxp[2],Sdrxp[3]);}
       
        Sdra[j] = min ( min ( Sdrxp[0],Sdrxp[1] ), Sdrxp[2] );
        Sdrc[j] = max ( max ( Sdrxp[0],Sdrxp[1] ), Sdrxp[2] );

    if(Sdrxp[0]!=Sdra[j] && Sdrxp[0]!=Sdrc[j]) {Sdrb[j] = Sdrxp[0];}
    else if(Sdrxp[1]!=Sdra[j] && Sdrxp[1]!=Sdrc[j]) {Sdrb[j] = Sdrxp[1];}
    else if(Sdrxp[2]!=Sdra[j] && Sdrxp[2]!=Sdrc[j]) {Sdrb[j] = Sdrxp[2];}            
    
    }
  
    //open a txt file to store the positions of the peaks
    ofstream Sd1r_alpha;
    Sd1r_alpha.open ( "Sdr2_AlphaPeaks.txt" );
    Sd1r_alpha << " strip / 1 / 2 / 3  \n";
  
    ofstream Sd1r_alphaCounts;
    Sd1r_alphaCounts.open ( "Sd2rCounts.txt" );
    Sd1r_alphaCounts << " strip / Counts 1 / 2 / 3  \n";
    for ( int i=0; i<24; i++ )
    {
      //cout << i << " a " << Sdra[i] << " b " << Sdrb[i] << " c " << Sdrc[i] << endl;
       if(i<6)
       {fit_func3->SetParLimits ( 1,Sdra[i]-5,Sdra[i]+5 );
        fit_func3->SetParLimits ( 4,Sdrb[i]-5,Sdrb[i]+5 );
        fit_func3->SetParLimits ( 6,Sdrc[i]-5,Sdrc[i]+5 );
        fit_func3->SetParLimits ( 2,0,10 );
       }
       else if(i>6)
       {fit_func3->SetParLimits ( 1,Sdra[i]-5,Sdra[i]+5 );
        fit_func3->SetParLimits ( 4,Sdrb[i]-5,Sdrb[i]+5 );
        fit_func3->SetParLimits ( 6,Sdrc[i]-5,Sdrc[i]+5 );
        fit_func3->SetParLimits ( 2,0,6 );
        }
        hSd2r[i]->Fit ( "fit_func3","","",Sdra[i]-10,Sdrc[i]+10 );

        Sd1r_alpha << i << " "  << fit_func3->GetParameter ( 1 ) << " " << fit_func3->GetParameter ( 4 ) << " " << fit_func3->GetParameter ( 6 ) << endl;
        Sd1r_alphaCounts << i << " " << fit_func3->GetParameter ( 0 ) << " " <<  fit_func3->GetParameter ( 3 ) << " " <<  fit_func3->GetParameter ( 5 ) << endl;
    }//end of fit on Sd1r detector
    
    //find peaks and fir the Sd1s detector
   //Search for peaks in Yd strips
    double Sdsxp[3];
    double Sdsa[32],Sdsb[32],Sdsc[32],Sdsped[32];

    for ( int j=0; j<32; j++ )
    {
        nfound = s->Search ( hSd2s[j],2,"",0.1 );
       
        for ( int p=0; p<nfound; p++ ){Sdsxp[p] = txp[p];}
        
      //  Sdsped[j] = min ( min ( Sdsxp[0], Sdsxp[1] ),min ( Sdsxp[2],Sdsxp[3] ) );
      //  Sdsc[j] = max ( max ( Sdsxp[0], Sdsxp[1] ),max ( Sdsxp[2],Sdsxp[3] ) );

     //  if((Sdsxp[0]!=Sdsped[j] && Sdsxp[0]!=Sdsc[j]) && (Sdsxp[1]!=Sdsped[j] && Sdsxp[1]!=Sdsc[j])){Sdsa[j]=min(Sdsxp[0],Sdsxp[1]);Sdsb[j]=max(Sdsxp[0],Sdsxp[1]);}
     //  else if ((Sdsxp[0]!=Sdsped[j] && Sdsxp[0]!=Sdsc[j]) && (Sdsxp[2]!=Sdsped[j] && Sdsxp[2]!=Sdsc[j])){Sdsa[j]=min(Sdsxp[0],Sdsxp[2]);Sdsb[j]=max(Sdsxp[0],Sdsxp[2]);}
     //  else if ((Sdsxp[0]!=Sdsped[j] && Sdsxp[0]!=Sdsc[j]) && (Sdsxp[3]!=Sdsped[j] && Sdsxp[3]!=Sdsc[j])){Sdsa[j]=min(Sdsxp[0],Sdsxp[3]);Sdsb[j]=max(Sdsxp[0],Sdsxp[3]);}
     //  else if ((Sdsxp[1]!=Sdsped[j] && Sdsxp[1]!=Sdsc[j]) && (Sdsxp[2]!=Sdsped[j] && Sdsxp[2]!=Sdsc[j])){Sdsa[j]=min(Sdsxp[1],Sdsxp[2]);Sdsb[j]=max(Sdsxp[1],Sdsxp[2]);}
     //  else if ((Sdsxp[1]!=Sdsped[j] && Sdsxp[1]!=Sdsc[j]) && (Sdsxp[3]!=Sdsped[j] && Sdsxp[3]!=Sdsc[j])){Sdsa[j]=min(Sdsxp[1],Sdsxp[3]);Sdsb[j]=max(Sdsxp[1],Sdsxp[3]);}
     //  else if ((Sdsxp[2]!=Sdsped[j] && Sdsxp[2]!=Sdsc[j]) && (Sdsxp[3]!=Sdsped[j] && Sdsxp[3]!=Sdsc[j])){Sdsa[j]=min(Sdsxp[2],Sdsxp[3]);Sdsb[j]=max(Sdsxp[2],Sdsxp[3]);}
      
        Sdsa[j] = min ( min ( Sdsxp[0],Sdsxp[1] ), Sdsxp[2] );
        Sdsc[j] = max ( max ( Sdsxp[0],Sdsxp[1] ), Sdsxp[2] );

    if(Sdsxp[0]!=Sdsa[j] && Sdsxp[0]!=Sdsc[j]) {Sdsb[j] = Sdsxp[0];}
    else if(Sdsxp[1]!=Sdsa[j] && Sdsxp[1]!=Sdsc[j]) {Sdsb[j] = Sdsxp[1];}
    else if(Sdsxp[2]!=Sdsa[j] && Sdsxp[2]!=Sdsc[j]) {Sdsb[j] = Sdsxp[2];}       
                           
    }
  
    //open a txt file to store the positions of the peaks
    ofstream Sd1s_alpha;
    Sd1s_alpha.open ( "Sds2_AlphaPeaks.txt" );
    Sd1s_alpha << " strip / 1 / 2 / 3  \n";
  
    ofstream Sd1s_alphaCounts;
    Sd1s_alphaCounts.open ( "Sd2sCounts.txt" );
    Sd1s_alphaCounts << " strip / Counts 1 / 2 / 3  \n";
    for ( int i=0; i<32; i++ )
    {
            fit_func3->SetParLimits ( 1,Sdsa[i]-5,Sdsa[i]+5 );
            fit_func3->SetParLimits ( 4,Sdsb[i]-5,Sdsb[i]+5 );
            fit_func3->SetParLimits ( 6,Sdsc[i]-5,Sdsc[i]+5 );
            fit_func3->SetParLimits ( 2,0,10 );

        
        hSd2s[i]->Fit ( "fit_func3","","",Sdsa[i]-10,Sdsc[i]+10 );

        Sd1s_alpha << i << " "  << fit_func3->GetParameter ( 1 ) << " " << fit_func3->GetParameter ( 4 ) << " " << fit_func3->GetParameter ( 6 ) << endl;
        Sd1s_alphaCounts << i << " " << fit_func3->GetParameter ( 0 ) << " " <<  fit_func3->GetParameter ( 3 ) << " " <<  fit_func3->GetParameter ( 5 ) << endl;
    }//end of fit on Sd1s detector
        */
  
    /*double xp[3];
    double a[128],b[128],c[128],ped[128];

    for ( int j=0; j<128; j++ )
    {
        nfound = s->Search ( hYu[j],2,"",0.1 );*/
//cout << "nfound for loop" << endl;
  
/*  for(int i=0; i<8; i++)
  {
    for(int j=0; j<16; j++)
    {
      nfound = s->Search ( hYuSectorSing[i][j],2,"",0.1 );*/
      
            //for ( int p=0; p<nfound; p++ ){xp[p] = txp[p];}
                
       // ped[j] = min ( min ( xp[0], xp[1] ),min ( xp[2],xp[3] ) );
       // c[j] = max ( max ( xp[0], xp[1] ),max ( xp[2],xp[3] ) );

      // if((xp[0]!=ped[j] && xp[0]!=c[j]) && (xp[1]!=ped[j] && xp[1]!=c[j])){a[j]=min(xp[0],xp[1]);b[j]=max(xp[0],xp[1]);}
      // else if ((xp[0]!=ped[j] && xp[0]!=c[j]) && (xp[2]!=ped[j] && xp[2]!=c[j])){a[j]=min(xp[0],xp[2]);b[j]=max(xp[0],xp[2]);}
     //  else if ((xp[0]!=ped[j] && xp[0]!=c[j]) && (xp[3]!=ped[j] && xp[3]!=c[j])){a[j]=min(xp[0],xp[3]);b[j]=max(xp[0],xp[3]);}
     //  else if ((xp[1]!=ped[j] && xp[1]!=c[j]) && (xp[2]!=ped[j] && xp[2]!=c[j])){a[j]=min(xp[1],xp[2]);b[j]=max(xp[1],xp[2]);}
     //  else if ((xp[1]!=ped[j] && xp[1]!=c[j]) && (xp[3]!=ped[j] && xp[3]!=c[j])){a[j]=min(xp[1],xp[3]);b[j]=max(xp[1],xp[3]);}
     //  else if ((xp[2]!=ped[j] && xp[2]!=c[j]) && (xp[3]!=ped[j] && xp[3]!=c[j])){a[j]=min(xp[2],xp[3]);b[j]=max(xp[2],xp[3]);}
       
 /*       a[(i*16)+j] = min ( min (xp[0],xp[1]), xp[2]);
        c[(i*16)+j] = max ( max (xp[0],xp[1]), xp[2]);
    
    if(xp[0]!=a[(i*16)+j] && xp[0]!=c[(i*16)+j]) {b[(i*16)+j] = xp[0];}
    else if(xp[1]!=a[(i*16)+j] && xp[1]!=c[(i*16)+j]) {b[(i*16)+j] = xp[1];}
    else if(xp[2]!=a[(i*16)+j] && xp[2]!=c[(i*16)+j]) {b[(i*16)+j] = xp[2];}*/
    
    /* a[j] = min ( min (xp[0],xp[1]), xp[2]);
     c[j] = max ( max (xp[0],xp[1]), xp[2]);
    
    if(xp[0]!=a[j] && xp[0]!=c[j]) {b[j] = xp[0];}
    else if(xp[1]!=a[j] && xp[1]!=c[j]) {b[j] = xp[1];}
    else if(xp[2]!=a[j] && xp[2]!=c[j]) {b[j] = xp[2];}
       
          //cout << j << "pedestal " << ped[j] << " peak 1 " << a[j] << " peak 2 " << b[j] << " peak 3 " << c[j] << endl;
          cout << j << " peak 1 " << a[j] << " peak 2 " << b[j] << " peak 3 " << c[j] << endl;*/
               
   // }
//  }
    //open a txt file to store the positions of the peaks
   /* ofstream Yu_alpha;
    Yu_alpha.open ( "Yu_AlphaPeaks5226.txt" );
    Yu_alpha << " strip / 1 / 2 / 3  \n";

  
    ofstream Yu_alphaCounts;
    Yu_alphaCounts.open ( "YuCounts5226.txt" );
    Yu_alphaCounts << " strip / Counts 1 / 2 / 3  \n";
      
    for ( int i=0; i<128; i++ )
    {
      cout << i << " /a " << a[i] << " /b " << b[i] << " /c " << c[i] << endl;
      
       
        fit_func6->SetParLimits ( 1,a[i]-30,a[i]-10 );
        fit_func6->SetParLimits ( 4,a[i]-10,a[i]+10 );
        fit_func6->SetParLimits ( 6,b[i]-30,b[i]-10 );
        fit_func6->SetParLimits ( 8,b[i]-10,b[i]+10 );
        fit_func6->SetParLimits ( 10,c[i]-30,c[i]-10 );
        fit_func6->SetParLimits ( 12,c[i]-10,c[i]+10 );
        fit_func6->SetParLimits ( 2,5,15 );*/
        
      /*if (40<i<127)  
          {   
        fit_func6->SetParLimits ( 1,a[i]-30,a[i]-10 );
        fit_func6->SetParLimits ( 4,a[i]-10,a[i]+10 );
        fit_func6->SetParLimits ( 6,b[i]-30,b[i]-10 );
        fit_func6->SetParLimits ( 8,b[i]-10,b[i]+10 );
        fit_func6->SetParLimits ( 10,c[i]-30,c[i]-10 );
        fit_func6->SetParLimits ( 12,c[i]-10,c[i]+10 );
        fit_func6->SetParLimits ( 2,0,15 );}*/
      
   //   fit_func3->SetParLimits(1,a[i]-30,a[i]+30);
   //   fit_func3->SetParLimits(4,b[i]-30,b[i]+30);
  //    fit_func3->SetParLimits(6,c[i]-30,c[i]+30);
  //    fit_func3->SetParLimits(2,0,30);
	
	/*fit_func9->SetParLimits ( 1,a[i]-45,a[i]-30 );
        fit_func9->SetParLimits ( 4,a[i]-30,a[i]-15 );
        fit_func9->SetParLimits ( 6,a[i]-15,a[i]+15 );
        fit_func9->SetParLimits ( 8,b[i]-45,b[i]-30 );
        fit_func9->SetParLimits ( 10,b[i]-30,b[i]-15 );
        fit_func9->SetParLimits ( 12,b[i]-15,b[i]+15 );
	fit_func9->SetParLimits ( 8,c[i]-45,c[i]-30 );
        fit_func9->SetParLimits ( 10,c[i]-30,c[i]-15 );
        fit_func9->SetParLimits ( 12,c[i]-15,c[i]+15 );
        fit_func9->SetParLimits ( 2,0,30 );*/
      
	
	    /*cout << " 1 " << a[i] << " 2 " << b[i] << " 3 " << c[i] << endl;
        hYu[i]->Fit ( "fit_func6","","",a[i]-100,c[i]+100 );
	
	//Yu_alpha << i  << " "  << fit_func9->GetParameter ( 6 ) << " " << fit_func9->GetParameter ( 12 ) << " " << fit_func9->GetParameter ( 18 ) << endl;
        //Yu_alphaCounts << i << " " << fit_func9->GetParameter (5 ) << " " <<  fit_func9->GetParameter ( 11 ) << " " <<  fit_func9->GetParameter ( 17 ) << endl;

        Yu_alpha << i  << " "  << fit_func6->GetParameter ( 4 ) << " " << fit_func6->GetParameter ( 8 ) << " " << fit_func6->GetParameter ( 12 ) << endl;
        Yu_alphaCounts << i << " " << fit_func6->GetParameter (3 ) << " " <<  fit_func6->GetParameter ( 7 ) << " " <<  fit_func6->GetParameter ( 11 ) << endl;*/
       // Yu_alpha << i  << " "  << fit_func3->GetParameter ( 1 ) << " " << fit_func3->GetParameter ( 4 ) << " " << fit_func3->GetParameter ( 6 ) << endl;
      //  Yu_alphaCounts << i << " " << fit_func3->GetParameter (0 ) << " " <<  fit_func3->GetParameter ( 3 ) << " " <<  fit_func3->GetParameter ( 5 ) << endl;
  //  }
  
  
   /* for ( int i=0; i<8; i++ )
    {
        for ( int j=0; j<16; j++ )
        {
   
        fit_func6->SetParLimits ( 1,a[(i*16)+j]-30,a[(i*16)+j]-10 );
        fit_func6->SetParLimits ( 4,a[(i*16)+j]-10,a[(i*16)+j]+10 );
        fit_func6->SetParLimits ( 6,b[(i*16)+j]-30,b[(i*16)+j]-10 );
        fit_func6->SetParLimits ( 8,b[(i*16)+j]-10,b[(i*16)+j]+10 );
        fit_func6->SetParLimits ( 10,c[(i*16)+j]-30,c[(i*16)+j]-10 );
        fit_func6->SetParLimits ( 12,c[(i*16)+j]-10,c[(i*16)+j]+10 );
        fit_func6->SetParLimits ( 2,0,15 );

            hYuSectorSing[i][j]->Fit ( "fit_func6","","",a[ ( i*16 ) +j]-100,c[ ( i*16 ) +j]+100 );
            Yu_alpha << (i*16)+j  << " "  << fit_func6->GetParameter ( 4 ) << " " << fit_func6->GetParameter ( 8 ) << " " << fit_func6->GetParameter ( 12 ) << endl;
            Yu_alphaCounts << (i*16)+j << " " << fit_func6->GetParameter ( 3 ) << " " <<  fit_func6->GetParameter ( 7 ) << " " <<  fit_func6->GetParameter ( 11 ) << endl;

        }
    }
    */
  //end of the fit on the Yd detector
  
  //find peaks and fir the Sd1r detector
   //Search for peaks in Yd strips
 /*   double Sdrxp[4];
    double Sdra[24],Sdrb[24],Sdrc[24],Sdrped[24];

    for ( int j=0; j<24; j++ )
    {
        nfound = s->Search ( hSur[j],2,"",0.1 );
       
        for ( int p=0; p<nfound; p++ ){Sdrxp[p] = txp[p];}
        
        Sdrped[j] = min ( min ( Sdrxp[0], Sdrxp[1] ),min ( Sdrxp[2],Sdrxp[3] ) );
        Sdrc[j] = max ( max ( Sdrxp[0], Sdrxp[1] ),max ( Sdrxp[2],Sdrxp[3] ) );

       if((Sdrxp[0]!=Sdrped[j] && Sdrxp[0]!=Sdrc[j]) && (Sdrxp[1]!=Sdrped[j] && Sdrxp[1]!=Sdrc[j])){Sdra[j]=min(Sdrxp[0],Sdrxp[1]);Sdrb[j]=max(Sdrxp[0],Sdrxp[1]);}
       else if ((Sdrxp[0]!=Sdrped[j] && Sdrxp[0]!=Sdrc[j]) && (Sdrxp[2]!=Sdrped[j] && Sdrxp[2]!=Sdrc[j])){Sdra[j]=min(Sdrxp[0],Sdrxp[2]);Sdrb[j]=max(Sdrxp[0],Sdrxp[2]);}
       else if ((Sdrxp[0]!=Sdrped[j] && Sdrxp[0]!=Sdrc[j]) && (Sdrxp[3]!=Sdrped[j] && Sdrxp[3]!=Sdrc[j])){Sdra[j]=min(Sdrxp[0],Sdrxp[3]);Sdrb[j]=max(Sdrxp[0],Sdrxp[3]);}
       else if ((Sdrxp[1]!=Sdrped[j] && Sdrxp[1]!=Sdrc[j]) && (Sdrxp[2]!=Sdrped[j] && Sdrxp[2]!=Sdrc[j])){Sdra[j]=min(Sdrxp[1],Sdrxp[2]);Sdrb[j]=max(Sdrxp[1],Sdrxp[2]);}
       else if ((Sdrxp[1]!=Sdrped[j] && Sdrxp[1]!=Sdrc[j]) && (Sdrxp[3]!=Sdrped[j] && Sdrxp[3]!=Sdrc[j])){Sdra[j]=min(Sdrxp[1],Sdrxp[3]);Sdrb[j]=max(Sdrxp[1],Sdrxp[3]);}
       else if ((Sdrxp[2]!=Sdrped[j] && Sdrxp[2]!=Sdrc[j]) && (Sdrxp[3]!=Sdrped[j] && Sdrxp[3]!=Sdrc[j])){Sdra[j]=min(Sdrxp[2],Sdrxp[3]);Sdrb[j]=max(Sdrxp[2],Sdrxp[3]);}
                           
    }
  
    //open a txt file to store the positions of the peaks
    ofstream Sur_alpha;
    Sur_alpha.open ( "Sur_AlphaPeaks225.txt" );
    Sur_alpha << " strip / pedestal / 1 / 2 / 3  \n";
  
    ofstream Sur_alphaCounts;
    Sur_alphaCounts.open ( "SurCounts225.txt" );
    Sur_alphaCounts << " strip / Counts 1 / 2 / 3  \n";
    for ( int i=0; i<24; i++ )
    {
        fit_func3->SetParLimits ( 1,Sdra[i]-20,Sdra[i]+20 );
            fit_func3->SetParLimits ( 4,Sdrb[i]-20,Sdrb[i]+20 );
            fit_func3->SetParLimits ( 6,Sdrc[i]-20,Sdrc[i]+20 );
            fit_func3->SetParLimits ( 2,0,25 );
        
        hSur[i]->Fit ( "fit_func3","","",Sdra[i]-50,Sdrc[i]+50 );

        Sur_alpha << i << " " << Sdrped[i] << " "  << fit_func3->GetParameter ( 1 ) << " " << fit_func3->GetParameter ( 4 ) << " " << fit_func3->GetParameter ( 6 ) << endl;
        Sur_alphaCounts << i << " " << fit_func3->GetParameter ( 0 ) << " " <<  fit_func3->GetParameter ( 3 ) << " " <<  fit_func3->GetParameter ( 5 ) << endl;
    }//end of fit on Sd1r detector
    
    //find peaks and fir the Sd1s detector
   //Search for peaks in Yd strips
    double Sdsxp[4];
    double Sdsa[32],Sdsb[32],Sdsc[32],Sdsped[32];

    for ( int j=0; j<32; j++ )
    {
        nfound = s->Search ( hSus[j],2,"",0.1 );
       
        for ( int p=0; p<nfound; p++ ){Sdsxp[p] = txp[p];}
        
        Sdsped[j] = min ( min ( Sdsxp[0], Sdsxp[1] ),min ( Sdsxp[2],Sdsxp[3] ) );
        Sdsc[j] = max ( max ( Sdsxp[0], Sdsxp[1] ),max ( Sdsxp[2],Sdsxp[3] ) );

       if((Sdsxp[0]!=Sdsped[j] && Sdsxp[0]!=Sdsc[j]) && (Sdsxp[1]!=Sdsped[j] && Sdsxp[1]!=Sdsc[j])){Sdsa[j]=min(Sdsxp[0],Sdsxp[1]);Sdsb[j]=max(Sdsxp[0],Sdsxp[1]);}
       else if ((Sdsxp[0]!=Sdsped[j] && Sdsxp[0]!=Sdsc[j]) && (Sdsxp[2]!=Sdsped[j] && Sdsxp[2]!=Sdsc[j])){Sdsa[j]=min(Sdsxp[0],Sdsxp[2]);Sdsb[j]=max(Sdsxp[0],Sdsxp[2]);}
       else if ((Sdsxp[0]!=Sdsped[j] && Sdsxp[0]!=Sdsc[j]) && (Sdsxp[3]!=Sdsped[j] && Sdsxp[3]!=Sdsc[j])){Sdsa[j]=min(Sdsxp[0],Sdsxp[3]);Sdsb[j]=max(Sdsxp[0],Sdsxp[3]);}
       else if ((Sdsxp[1]!=Sdsped[j] && Sdsxp[1]!=Sdsc[j]) && (Sdsxp[2]!=Sdsped[j] && Sdsxp[2]!=Sdsc[j])){Sdsa[j]=min(Sdsxp[1],Sdsxp[2]);Sdsb[j]=max(Sdsxp[1],Sdsxp[2]);}
       else if ((Sdsxp[1]!=Sdsped[j] && Sdsxp[1]!=Sdsc[j]) && (Sdsxp[3]!=Sdsped[j] && Sdsxp[3]!=Sdsc[j])){Sdsa[j]=min(Sdsxp[1],Sdsxp[3]);Sdsb[j]=max(Sdsxp[1],Sdsxp[3]);}
       else if ((Sdsxp[2]!=Sdsped[j] && Sdsxp[2]!=Sdsc[j]) && (Sdsxp[3]!=Sdsped[j] && Sdsxp[3]!=Sdsc[j])){Sdsa[j]=min(Sdsxp[2],Sdsxp[3]);Sdsb[j]=max(Sdsxp[2],Sdsxp[3]);}
                           
    }
  
    //open a txt file to store the positions of the peaks
    ofstream Sus_alpha;
    Sus_alpha.open ( "Sus_AlphaPeaks225.txt" );
    Sus_alpha << " strip / pedestal / 1 / 2 / 3  \n";
  
    ofstream Sus_alphaCounts;
    Sus_alphaCounts.open ( "SusCounts225.txt" );
    Sus_alphaCounts << " strip / Counts 1 / 2 / 3  \n";
    for ( int i=0; i<32; i++ )
    {
            fit_func3->SetParLimits ( 1,Sdsa[i]-10,Sdsa[i]+10 );
            fit_func3->SetParLimits ( 4,Sdsb[i]-10,Sdsb[i]+10 );
            fit_func3->SetParLimits ( 6,Sdsc[i]-10,Sdsc[i]+10 );
            fit_func3->SetParLimits ( 2,0,20 );

        
        hSus[i]->Fit ( "fit_func3","","",Sdsa[i]-50,Sdsc[i]+50 );

        Sus_alpha << i << " " << Sdsped[i] << " "  << fit_func3->GetParameter ( 1 ) << " " << fit_func3->GetParameter ( 4 ) << " " << fit_func3->GetParameter ( 6 ) << endl;
        Sus_alphaCounts << i << " " << fit_func3->GetParameter ( 0 ) << " " <<  fit_func3->GetParameter ( 3 ) << " " <<  fit_func3->GetParameter ( 5 ) << endl;
    }//end of fit on Sus detector
        */
    f_out->cd();

    //for ( int i=0; i<128; i++ ){hYd[i]->Write();}
    for ( int i=0; i<128; i++ )
    {
    //hYd[i]->GetXaxis()->SetRange(0,30);
    hYd[i]->Write();
    }
   // for(int i=0; i<8; i++){for(int j=0; j<16; j++){hYuSectorRing[i][j]->Write();}}
    
 //   for ( int i=0; i<16; i++ ){hCsI1[i]->Write();}
 //   for ( int i=0; i<16; i++ ){hCsI2[i]->Write();}
    
    //for ( int i=0; i<24; i++ ){hSd1r[i]->Write();}
  //  for ( int i=0; i<24; i++ ){hSd2r[i]->Write();}
    
   // for ( int i=0; i<32; i++ ){hSd1s[i]->Write();}
  //  for ( int i=0; i<32; i++ ){hSd2s[i]->Write();}
    
   // for ( int i=0; i<24; i++ ){hSur[i]->Write();}
   // for ( int i=0; i<32; i++ ){hSus[i]->Write();}
        

    f_out->Close();

    return 0;
}
