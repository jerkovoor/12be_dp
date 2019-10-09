//M. Vostinar
//The code takes the files which were already treated with decode.cxx
//It is used for calibration purposes
//Fils in the corresponding histograms and uses TSpectrum to search for peaks
//After TSpectrum has found the peaks, depending on the judgment of person fiting one can choose a 3,6 or 9 gauss fit to use on the histogrames
//The histogrames with the fit are saved in a .root file and 2 txt. files contain the position of the tree interesting peaks and the amplitude of these peaks
//The code can be used on all detectors

using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TSpectrum.h"
#include "iostream"
#include "fstream"
#include "TH1.h"
#include "TF1.h"


double readAlpha()
{
  //definition of variables necessary for reading the input root file
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
    int ICChannel;
    double ICEnergyRaw;

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
    vector<double>* Sd2rEnergy= new vector<double>;
 
    //Sd2s detector
    int Sd2sMul;
    vector<int>* Sd2sChannel= new vector<int>; //Sector Channels
    vector<double>* Sd2sEnergyRaw= new vector<double>;
    vector<double>* Sd2sEnergy= new vector<double>;

    //Sur
    int SurMul;
    vector<int>* SurChannel= new vector<int>; //Ring Channels
    vector<double>* SurEnergyRaw= new vector<double>;

    //Sus
    int SusMul;
    vector<int>* SusChannel= new vector<int>; //Sector Channels
    vector<double>* SusEnergyRaw= new vector<double>;
  
    TChain* chain = new TChain ( "AutoTree" );

    //Open the calibration files 
    // Yd calibration before the experiment
    chain->Add ( "/home/jerome/12Be_exp/Processed_files/C_notarget/decode4992.root" ); //change the path to the files and the file name accordingly
    /*chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Alpha/Decode_4973.root" );
    chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Alpha/Decode_4974.root" );
    chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Alpha/Decode_4975.root" );
    chain->Add ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Alpha/Decode_4977.root" );*/
  
  //the 5225 Yu calibration run
 /* chain->Add("/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Alpha/Decode_5225_1.root");
  chain->Add("/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Alpha/Decode_5225_2.root");
  chain->Add("/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Alpha/Decode_5225_3.root");*/

  //read out the input tree
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
    
    //definition of various histogrames
    //one histograme is created for each strip to be able to perform the peak search for calibration
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
        hYu[i] =  new TH1D ( Yu_name.c_str(),"Yu_Ch",175,2400,3800); //modified to see if the calib is better 

    }
    
    TH1D *hYuSectorRing[8][16]; //These histogrames indicate the sector and the ring of the Yu detector. The purpose is to investigate trends conected to sectors/rings
    
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
        hSd1r[i] =  new TH1D ( Sd1r_name.c_str(),"Sd1r_Ch",1024,0,4096);

        string Sd2r_name = Form ( "Sd2r_Ch%i",i );
        hSd2r[i] =  new TH1D ( Sd2r_name.c_str(),"Sd2r_Ch",1024,0,4096 );

        string Sur_name = Form ( "Sur_Ch%i",i );
        hSur[i] =  new TH1D ( Sur_name.c_str(),"Sur_Ch",1024,0,4096 );
    }
    
    for(int i=0; i<32; i++)
    {
        string Sd1s_name = Form ( "Sd1s_Ch%i",i );
        hSd1s[i] =  new TH1D ( Sd1s_name.c_str(),"Sd1s_Ch",1024,0,4096 );

        string Sd2s_name = Form ( "Sd2s_Ch%i",i );
        hSd2s[i] =  new TH1D ( Sd2s_name.c_str(),"Sd2s_Ch",1024,0,4096 );

        string Sus_name = Form ( "Sus_Ch%i",i );
        hSus[i] =  new TH1D ( Sus_name.c_str(),"Sus_Ch",1024,0,4096 );
    }
    
    TFile *f_out = new TFile ( "/home/jerome/12Be_exp/notarget_peak_position/4992_notarget.root","RECREATE" ); //Creating the output .root file put the corect path and/or name
    
    int ev = chain->GetEntries(); //get the total number of entries
    cout << "Total number of events =" << ev << endl;
    int ev_num=0;
    
    while ( ev_num<=ev ) //the main while loop, loops over all the events in the tree and does whatever needs to be done 
    {
        if ( ev_num%10000==0 ) cout << "Current event = " << ev_num << "\r"<< flush;
        chain->GetEntry ( ev_num ); //get the current event

	//the filing of the histograms is performed in the same way for all different detectors, requireing that a channel should or shouldn't be >500 to evade confusion of the TSpectrum 
       if ( YdEnergyRaw->size() >0){if ( YdEnergyRaw->at(0)>500 ){hYd[YdChannel->at(0)]->Fill ( YdEnergyRaw->at(0) );}} //requires that the vector is acctually filled and that the channel is >500 so not to take into acount the pedestal, and estimate only the position of the alpha peaks
     
       if ( YuEnergyRaw->size() >0){if(YuEnergyRaw->at(0)>0)
	 {
	    hYu[YuChannel->at ( 0 )]->Fill ( YuEnergyRaw->at ( 0 ) );
	    hYuSectorRing[YuSector->at(0)][YuRing->at(0)]->Fill(YuEnergyRaw->at ( 0 ) ); 
	 }}
         
     //   if ( CsI1EnergyRaw->size() >0 && CsI1Mul==1 ){if ( CsI1EnergyRaw->at ( 0 ) >0 ){hCsI1[CsI1Channel->at ( 0 )]->Fill ( CsI1EnergyRaw->at ( 0 ) );}} //CsI is not calibrated using the Alpha runs
     //   if ( CsI2EnergyRaw->size() >0 && CsI2Mul==1 ){if ( CsI2EnergyRaw->at ( 0 ) >0 ){hCsI2[CsI2Channel->at ( 0 )]->Fill ( CsI2EnergyRaw->at ( 0 ) );}}
        if ( Sd1rEnergyRaw->size() >0){if ( Sd1rEnergyRaw->at ( 0 ) >0 ){hSd1r[Sd1rChannel->at ( 0 )]->Fill ( Sd1rEnergyRaw->at ( 0 ) );}}  
        if ( Sd1sEnergyRaw->size() >0){if ( Sd1sEnergyRaw->at ( 0 ) >0 ){hSd1s[Sd1sChannel->at ( 0 )]->Fill ( Sd1sEnergyRaw->at ( 0 ) );}}
        if ( Sd2rEnergyRaw->size() >0){if ( Sd2rEnergyRaw->at ( 0 ) >0 ){hSd2r[Sd2rChannel->at ( 0 )]->Fill ( Sd2rEnergyRaw->at ( 0 ) );}}
        if ( Sd2sEnergyRaw->size() >0){if ( Sd2sEnergyRaw->at ( 0 ) >0 ){hSd2s[Sd2sChannel->at ( 0 )]->Fill ( Sd2sEnergyRaw->at ( 0 ) );}}
        if ( SurEnergyRaw->size() >0){if ( SurEnergyRaw->at ( 0 ) >0 ){hSur[SurChannel->at ( 0 )]->Fill ( SurEnergyRaw->at ( 0 ) );}}
        if ( SusEnergyRaw->size() >0){if ( SusEnergyRaw->at ( 0 ) >0 ){hSus[SusChannel->at ( 0 )]->Fill ( SusEnergyRaw->at ( 0 ) );}}

        ev_num++;
    }//end of the while loop
    
    //definition of convoluted gaussian functions
    //1 gaussian
    TF1* fit_func1 = new TF1 ( "fit_func1","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",0,5000);
    TF1* fit_func2 = new TF1 ( "fit_func2","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",0,5000);
    //6 gaussians 
    TF1* fit_func6 = new TF1 ( "fit_func6","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp(-pow((x-[10]),2)/(2*[2]*[2]))+[11]*TMath::Exp(-pow((x-[12]),2)/(2*[2]*[2]))",500,5000 ); //removed the liner part for comparison
    //3 gaussians
    TF1* fit_func3 = new TF1 ( "fit_func3","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*x+[8]",0,5000 );
    //9 gaussians with linear background
    TF1* fit_func9 = new TF1 ( "fit_func9","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp(-pow((x-[10]),2)/(2*[2]*[2]))+[11]*TMath::Exp(-pow((x-[12]),2)/(2*[2]*[2]))+[13]*TMath::Exp(-pow((x-[14]),2)/(2*[2]*[2]))+[15]*TMath::Exp(-pow((x-[16]),2)/(2*[2]*[2]))+[17]*TMath::Exp(-pow((x-[18]),2)/(2*[2]*[2]))+[19]+[20]*x",2000,5000 );
    
    //Use TSpectrum to find the position of the peaks and limit the fit parameters to extract the exact positions of max for each peak
    //Search for peaks in Yd strips
    int npeaks=1; //number of peaks expected to see
    TSpectrum* s = new TSpectrum ( 2*npeaks );
    int nfound=0;
    double* txp = s->GetPositionX(); //returns the channel in which the peak has a max value
    
    /*double xp[3]; //this will store the position of the max   //remove here
    double a[128],b[128],c[128], ped[128]; //these will store the position of the max in correct order

   //search in the YY detectors
    for ( int j=0; j<128; j++ ) //loop over all the strips of the YY detectors
    {
        nfound = s->Search ( hYd[j],2,"",0.1 ); //change the name of the histogram acordingly

        for ( int p=0; p<nfound; p++ ){ xp[p] = txp[p];} //loop over the found peaks and fill in the xp variable
        
        //Treatment in case you are retriving the position of the pedestals at the same time
       // ped[j] = min ( min ( xp[0], xp[1] ),min ( xp[2],xp[3] ) ); //find the min value and assign it to the pedestal
       // c[j] = max ( max ( xp[0], xp[1] ),max ( xp[2],xp[3] ) ); //find the max value and asign it to c

       //search for the parameters a and b which represent the positions of the two peaks between the pedestal and the last alpha peak
      // if((xp[0]!=ped[j] && xp[0]!=c[j]) && (xp[1]!=ped[j] && xp[1]!=c[j])){a[j]=min(xp[0],xp[1]);b[j]=max(xp[0],xp[1]);}
      // else if ((xp[0]!=ped[j] && xp[0]!=c[j]) && (xp[2]!=ped[j] && xp[2]!=c[j])){a[j]=min(xp[0],xp[2]);b[j]=max(xp[0],xp[2]);}
     //  else if ((xp[0]!=ped[j] && xp[0]!=c[j]) && (xp[3]!=ped[j] && xp[3]!=c[j])){a[j]=min(xp[0],xp[3]);b[j]=max(xp[0],xp[3]);}
     //  else if ((xp[1]!=ped[j] && xp[1]!=c[j]) && (xp[2]!=ped[j] && xp[2]!=c[j])){a[j]=min(xp[1],xp[2]);b[j]=max(xp[1],xp[2]);}
     //  else if ((xp[1]!=ped[j] && xp[1]!=c[j]) && (xp[3]!=ped[j] && xp[3]!=c[j])){a[j]=min(xp[1],xp[3]);b[j]=max(xp[1],xp[3]);}
     //  else if ((xp[2]!=ped[j] && xp[2]!=c[j]) && (xp[3]!=ped[j] && xp[3]!=c[j])){a[j]=min(xp[2],xp[3]);b[j]=max(xp[2],xp[3]);}
       
       //search without the pedestal search
        a[j] = min ( min (xp[0],xp[1]), xp[2]); //find the min value
        c[j] = max ( max (xp[0],xp[1]), xp[2]); //find the maximum value
    
    //find the middle peak
    if(xp[0]!=a[j] && xp[0]!=c[j]) {b[j] = xp[0];}
    else if(xp[1]!=a[j] && xp[1]!=c[j]) {b[j] = xp[1];}
    else if(xp[2]!=a[j] && xp[2]!=c[j]) {b[j] = xp[2];}
       
      // cout << j << "pedestal " << ped[j] << " peak 1 " << a[j] << " peak 2 " << b[j] << " peak 3 " << c[j] << endl;
               
    }
  
    //open a txt file to store the positions of the peaks
    ofstream YY_alpha;
    YY_alpha.open ( "Yd_AlphaPeaks.txt" ); //change the path and the name of the file accordingly; these are for the Yd and Yu detectors
    YY_alpha << " Strip / Peak 1 / 2 / 3  \n";

  //open a file to store the amplitude of the peaks
    ofstream YY_alphaCounts;
    YY_alphaCounts.open ( "YdCounts.txt" ); //change the path and the name of the file acordingly; these are for Yd and Yu detectors
    YY_alphaCounts << " Strip / Amplitude 1 / 2 / 3  \n";
      
    for ( int i=0; i<128; i++ ) //loop over all the Yu/Yd strips to fit the data
    { //chosee the desired function to fit with and set the parameters
      
        //One starts from the maximum of the interesting peak and setts the positions of the small peaks on the left hand side of that peak
           /* fit_func9->SetParLimits ( 1,a[i]-100,a[i]-50 );
            fit_func9->SetParLimits ( 4,a[i]-30,a[i]-10 );
            fit_func9->SetParLimits ( 6,a[i]-10,a[i]+10 );
            fit_func9->SetParLimits ( 8,b[i]-100,b[i]-50 );
            fit_func9->SetParLimits ( 10,b[i]-30,b[i]-10 );
            fit_func9->SetParLimits ( 12,b[i]-10,b[i]+10 );
            fit_func9->SetParLimits ( 14,c[i]-100,c[i]-50 );
            fit_func9->SetParLimits ( 16,c[i]-30,c[i]-10 );
            fit_func9->SetParLimits ( 18,c[i]-10,c[i]+10 );
            fit_func9->SetParLimits ( 2,0,25 );*/
	   
	   /*if(i<25)
	   {
	    fit_func6->SetParLimits ( 1,a[i]-30,a[i]-10 );
            fit_func6->SetParLimits ( 4,a[i]-10,a[i]+10 );
            fit_func6->SetParLimits ( 6,b[i]-30,b[i]-10 );
            fit_func6->SetParLimits ( 8,b[i]-10,b[i]+10 );
            fit_func6->SetParLimits ( 10,c[i]-30,c[i]-10 );
            fit_func6->SetParLimits ( 12,c[i]-10,c[i]+10 );
            fit_func6->SetParLimits ( 2,0,10 );
	   }
	   
	   else if(i>24)
	   {*/
            /*fit_func6->SetParLimits ( 1,a[i]-30,a[i]-10 );   //remove here
            fit_func6->SetParLimits ( 4,a[i]-10,a[i]+10 );
            fit_func6->SetParLimits ( 6,b[i]-30,b[i]-10 );
            fit_func6->SetParLimits ( 8,b[i]-10,b[i]+10 );
            fit_func6->SetParLimits ( 10,c[i]-30,c[i]-10 );
            fit_func6->SetParLimits ( 12,c[i]-10,c[i]+10 );
            fit_func6->SetParLimits ( 2,0,25 );
	  // }
	   
            hYd[i]->Fit ( "fit_func6","","",a[i]-100,c[i]+100 ); //fit the desired histograme; change the name to Yu or Yd depending on interest

            //fill the .txt files
         //   YY_alpha << i << " "  << fit_func9->GetParameter ( 6 ) << " " << fit_func9->GetParameter ( 12 ) << " " << fit_func9->GetParameter ( 18 ) << endl;
         //   YY_alphaCounts << i << " " << fit_func9->GetParameter ( 5 ) << " " <<  fit_func9->GetParameter ( 11 ) << " " <<  fit_func9->GetParameter ( 17 ) << endl;
	 
	    YY_alpha << i << " "  << fit_func6->GetParameter ( 4 ) << " " << fit_func6->GetParameter ( 8 ) << " " << fit_func6->GetParameter ( 12 ) << endl;
            YY_alphaCounts << i << " " << fit_func6->GetParameter ( 3 ) << " " <<  fit_func6->GetParameter ( 7 ) << " " <<  fit_func6->GetParameter ( 11 ) << endl;
           
    } */ //end of the fit on the Yd/Yu  detectors  //remove here
 
   
   
   
   
   
   
   
   
   
  //find peaks and for the Sdr detectors
    double Sdrxp[1];
    //double Sdra[24],Sdrb[24],Sdrc[24];
    double Sdra[24];

    for ( int j=0; j<24; j++ )
    {
        nfound = s->Search ( hSd1r[j],2,"",0.1 ); //search for peaks in the hSdr detectors change the name of the histograme according to which peaks you want to search for Sd1r, Sd2r or Sur
       
        for ( int p=0; p<nfound; p++ ){Sdrxp[p] = txp[p];} 
        
        /*Sdra[j] = min ( min ( Sdrxp[0],Sdrxp[1] ), Sdrxp[2] );
        Sdrc[j] = max ( max ( Sdrxp[0],Sdrxp[1] ), Sdrxp[2] );

    if(Sdrxp[0]!=Sdra[j] && Sdrxp[0]!=Sdrc[j]) {Sdrb[j] = Sdrxp[0];} 
    else if(Sdrxp[1]!=Sdra[j] && Sdrxp[1]!=Sdrc[j]) {Sdrb[j] = Sdrxp[1];}
    else if(Sdrxp[2]!=Sdra[j] && Sdrxp[2]!=Sdrc[j]) {Sdrb[j] = Sdrxp[2];}  */          
    
    }
  
    //open a txt file to store the positions of the peaks
    ofstream Sdr_alpha;
    Sdr_alpha.open ( "/home/jerome/12Be_exp/notarget_peak_position/Sd1r_notarget_peak.txt" ); //change the path and the name of the file acordingly
    //Sdr_alpha << " Strip / Peak 1 / 2 / 3  \n";
    Sdr_alpha << " Strip / Peak 1 \n";
  
    ofstream Sdr_alphaCounts;
    Sdr_alphaCounts.open ( "/home/jerome/12Be_exp/notarget_peak_position/Sd1r_notarget_counts.txt" ); //change the path and the name of the file acordingly
    //Sdr_alphaCounts << " Strip / Amplitude 1 / 2 / 3  \n";
    Sdr_alphaCounts << " Strip / Amplitude 1 \n";
    
    for ( int i=0; i<24; i++ ) //loop over all the rings of S3 detectors to set the parameters and fit
    {
     //  if(i<6)
      // {
	/*fit_func3->SetParLimits ( 1,Sdra[i]-5,Sdra[i]+5 ); //remove here
        fit_func3->SetParLimits ( 4,Sdrb[i]-5,Sdrb[i]+5 );
        fit_func3->SetParLimits ( 6,Sdrc[i]-5,Sdrc[i]+5 );
        fit_func3->SetParLimits ( 2,0,10 );*/			//remove here
        
        fit_func1->SetParLimits ( 1,Sdrxp[i]-200,Sdrxp[i]+200 ); //2000,2500
        fit_func1->SetParLimits ( 2,0,200 ); //0,200
     //  }
     //  else if(i>5)
     //  {fit_func3->SetParLimits ( 1,Sdra[i]-3,Sdra[i]+3 );
     //   fit_func3->SetParLimits ( 4,Sdrb[i]-3,Sdrb[i]+3 );
     //   fit_func3->SetParLimits ( 6,Sdrc[i]-3,Sdrc[i]+3 );
     //   fit_func3->SetParLimits ( 2,0,3 );
     //   }
        hSd1r[i]->Fit ( "fit_func1","","",Sdrxp[i]-500,Sdrxp[i]+500 ); //change the histograme accordingly // 1000,3000

        //fill in the .txt files
        //Sdr_alpha << i << " "  << fit_func3->GetParameter ( 1 ) << " " << fit_func3->GetParameter ( 4 ) << " " << fit_func3->GetParameter ( 6 ) << endl;
	Sdr_alpha << i << " "  << fit_func1->GetParameter ( 1 ) << endl;
       // Sdr_alphaCounts << i << " " << fit_func3->GetParameter ( 0 ) << " " <<  fit_func3->GetParameter ( 3 ) << " " <<  fit_func3->GetParameter ( 5 ) << endl;
	Sdr_alphaCounts << i << " " << fit_func1->GetParameter ( 0 ) << endl;
    }//end of fit on Sd1r detector
    
    
    
    
    
    
    
    
    //find peaks and for the Sds detectors
    double Sdsxp[1];
    //double Sdsa[32],Sdsb[32],Sdsc[32],Sdsped[32];

    for ( int j=0; j<32; j++ ) //loop over the sectors of the S3 detectors
    {
        nfound = s->Search ( hSd1s[j],2,"",0.1 ); //change the name of the hitogram according to which detector you want to look at: Sd1s, Sd2s or Sus
       
        for ( int p=0; p<nfound; p++ ){Sdsxp[p] = txp[p];}
           
       /* Sdsa[j] = min ( min ( Sdsxp[0],Sdsxp[1] ), Sdsxp[2] );
        Sdsc[j] = max ( max ( Sdsxp[0],Sdsxp[1] ), Sdsxp[2] );

    if(Sdsxp[0]!=Sdsa[j] && Sdsxp[0]!=Sdsc[j]) {Sdsb[j] = Sdsxp[0];}
    else if(Sdsxp[1]!=Sdsa[j] && Sdsxp[1]!=Sdsc[j]) {Sdsb[j] = Sdsxp[1];}
    else if(Sdsxp[2]!=Sdsa[j] && Sdsxp[2]!=Sdsc[j]) {Sdsb[j] = Sdsxp[2];}     */  
                           
    }
  
    //open a txt file to store the positions of the peaks
    ofstream Sds_alpha;
    Sds_alpha.open ( "/home/jerome/12Be_exp/notarget_peak_position/Sd1s_notarget_peak.txt" );//change the path and the name of the file acordingly
    //Sds_alpha << " Strip / Peak 1 / 2 / 3  \n";
    Sds_alpha << " Strip / Peak 1 \n";
  
    ofstream Sds_alphaCounts;
    Sds_alphaCounts.open ( "/home/jerome/12Be_exp/notarget_peak_position/Sd1s_notarget_counts.txt" );//change the path and the name of the file acordingly
    //Sds_alphaCounts << " Strip / Amplitude 1 / 2 / 3  \n";
    Sds_alphaCounts <<  " Strip / Amplitude 1 \n";
    
    for ( int i=0; i<32; i++ ) //loop over the sectors of S3 detectors to see the parameters and fit
    {
     // if(i<16)
     // {   
	/*fit_func3->SetParLimits ( 1,Sdsa[i]-5,Sdsa[i]+5 );
	fit_func3->SetParLimits ( 4,Sdsb[i]-5,Sdsb[i]+5 );
        fit_func3->SetParLimits ( 6,Sdsc[i]-5,Sdsc[i]+5 );
        fit_func3->SetParLimits ( 2,0,10 );*/
   //   }
        fit_func2->SetParLimits ( 1,Sdsxp[i]-200,Sdsxp[i]+200 );
        fit_func2->SetParLimits ( 2,0,200 ); 
     //       else if(i>15)
	//    {
	  //   fit_func3->SetParLimits ( 1,Sdsa[i]-3,Sdsa[i]+3 );
          //   fit_func3->SetParLimits ( 4,Sdsb[i]-3,Sdsb[i]+3 );
          //   fit_func3->SetParLimits ( 6,Sdsc[i]-3,Sdsc[i]+3 );
          //   fit_func3->SetParLimits ( 2,0,3 );
	  //  }

        
        hSd1s[i]->Fit ( "fit_func2","","",Sdsxp[i]-500,Sdsxp[i]+500 ); //change the histograme name acordingly

        //fill in the .txt files
        //Sds_alpha << i << " "  << fit_func3->GetParameter ( 1 ) << " " << fit_func3->GetParameter ( 4 ) << " " << fit_func3->GetParameter ( 6 ) << endl;
        Sds_alpha << i << " "  << fit_func2->GetParameter ( 1 ) << endl;
	//Sds_alphaCounts << i << " " << fit_func3->GetParameter ( 0 ) << " " <<  fit_func3->GetParameter ( 3 ) << " " <<  fit_func3->GetParameter ( 5 ) << endl;
	Sds_alphaCounts << i << " " << fit_func2->GetParameter ( 0 ) <<  endl;
    }//end of fit on Sd1s detector
        
  
   f_out->cd();

   //Turn on the histogrames you want to write in the .root file
    //for ( int i=0; i<128; i++ ){hYd[i]->Write();}
 //   for ( int i=0; i<128; i++ ){hYu[i]->Write();}
   // for(int i=0; i<8; i++){for(int j=0; j<16; j++){hYuSectorRing[i][j]->Write();}}
    
 //   for ( int i=0; i<16; i++ ){hCsI1[i]->Write();}
 //   for ( int i=0; i<16; i++ ){hCsI2[i]->Write();}
    
    for ( int i=0; i<24; i++ ){hSd1r[i]->Write();}
    //for ( int i=0; i<24; i++ ){hSd2r[i]->Write();}
    
    for ( int i=0; i<32; i++ ){hSd1s[i]->Write();}
   //for ( int i=0; i<32; i++ ){hSd2s[i]->Write();}
    
   // for ( int i=0; i<24; i++ ){hSur[i]->Write();}
   // for ( int i=0; i<32; i++ ){hSus[i]->Write();}
        

    f_out->Close();

    return 0;
} //end of readAlpha.cxx