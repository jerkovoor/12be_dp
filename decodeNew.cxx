using namespace std;

#include <TFile.h>
#include "TString.h"
#include "string"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "iostream"
#include "TMath.h"
#include "/home/jerome/Software/TRIUMF/treeIris/include/IDet.h"
#include <vector>
#include <TObject.h>
#include <TClass.h>
#include "sstream"
#include "fstream"
#include "cstdlib"
#include <algorithm> 

//double decode ( int run_num, int index)
double decode ( int run_num)
{
      
    //string f_name = Form ( "/home/jerome/12Be_exp/Raw_root/tree%i_%d.root",run_num,index); 
    string f_name = Form ( "/home/jerome/12Be_exp/Raw_root/12C/tree%i.root",run_num); 
    TFile* f_in = TFile::Open ( f_name.c_str(),"READ" );
    TTree* tr_in = ( TTree* ) f_in->Get ( "Iris" );

    IDet* det = new IDet();

    Int_t* TYdMul = & ( det->TYdMul );
    std::vector<Int_t>* TYdChannel = & ( det->TYdChannel );
    std::vector<Double_t>* TYdEnergy= & ( det->TYdEnergy );
    std::vector<Int_t>* TYdADC= & ( det->TYdADC );
    std::vector<Int_t>* TYdNo= & ( det->TYdNo );
    std::vector<Int_t>* TYdRing= & ( det->TYdRing );
    std::vector<Double_t>* TYdTheta= & ( det->TYdTheta ); // Yd theta angle

    Int_t* TCsI1Mul= & ( det->TCsI1Mul );
    std::vector<Int_t>* TCsI1Channel= & ( det->TCsI1Channel );
    std::vector<Double_t>* TCsI1Energy= & ( det->TCsI1Energy );
    std::vector<Double_t>* TCsI1Phi= & ( det->TCsI1Phi );
    std::vector<Double_t>* TCsI1ADC= & ( det->TCsI1ADC );

    Int_t* TCsI2Mul= & ( det->TCsI2Mul );
    std::vector<Int_t>* TCsI2Channel= & ( det->TCsI2Channel );
    std::vector<Double_t>* TCsI2Energy= & ( det->TCsI2Energy );
    std::vector<Double_t>* TCsI2Phi= & ( det->TCsI2Phi );
    std::vector<Double_t>* TCsI2ADC= & ( det->TCsI2ADC );

    Double_t* TYdCsI1ETot= & ( det->TYdCsI1ETot );
    Double_t* TYdCsI2ETot= & ( det->TYdCsI2ETot );

    Int_t* TSSBADC= & ( det->TSSBADC );
    Double_t* TSSBEnergy= & ( det->TSSBEnergy );

    Int_t* TScADC= & ( det->TScADC );
    Double_t* TScEnergy= & ( det->TScEnergy );

    std::vector<Int_t>* TTrADC= & ( det->TTrADC );
    std::vector<Double_t>* TTrEnergy= & ( det->TTrEnergy );

    std::vector<Int_t>* TICChannel= & ( det->TICChannel );
    std::vector<Double_t>* TICEnergy= & ( det->TICEnergy );
    std::vector<Double_t>* TICADC= & ( det->TICADC );

    Int_t* TSd1rMul= & ( det->TSd1rMul );
    std::vector<Int_t>* TSd1rChannel= & ( det->TSd1rChannel );
    std::vector<Double_t>* TSd1rEnergy= & ( det->TSd1rEnergy );
    std::vector<Int_t>* TSd1rADC= & ( det->TSd1rADC );

    Int_t* TSd1sMul= & ( det->TSd1sMul );
    std::vector<Int_t>* TSd1sChannel= & ( det->TSd1sChannel );
    std::vector<Double_t>* TSd1sEnergy= & ( det->TSd1sEnergy );
    std::vector<Int_t>* TSd1sADC= & ( det->TSd1sADC );

    Int_t* TSd2rMul= & ( det->TSd2rMul );
    std::vector<Int_t>* TSd2rChannel= & ( det->TSd2rChannel );
    std::vector<Double_t>* TSd2rEnergy= & ( det->TSd2rEnergy );
    std::vector<Int_t>* TSd2rADC= & ( det->TSd2rADC );
    Double_t* TSd2rEnergyCal= & ( det->TSd2rEnergyCal );

    Int_t* TSd2sMul= & ( det->TSd2sMul );
    std::vector<Int_t>* TSd2sChannel= & ( det->TSd2sChannel );
    std::vector<Double_t>* TSd2sEnergy= & ( det->TSd2sEnergy );
    std::vector<Int_t>* TSd2sADC= & ( det->TSd2sADC );

    Double_t* TSdETot= & ( det->TSdETot );
    std::vector<Double_t>* TSd1Theta= & ( det->TSd1Theta );
    std::vector<Double_t>* TSd2Theta= & ( det->TSd2Theta );
    Double_t* TSdThetaCM= & ( det->TSdThetaCM );
    std::vector<Double_t>* TSd1Phi= & ( det->TSd1Phi );
    std::vector<Double_t>* TSd2Phi= & ( det->TSd2Phi );

    Int_t* TYuMul= & ( det->TYuMul );
    std::vector<Int_t>* TYuChannel= & ( det->TYuChannel );
    std::vector<Double_t>* TYuEnergy= & ( det->TYuEnergy );
    std::vector<Int_t>* TYuADC= & ( det->TYuADC );
    std::vector<Int_t>* TYuNo= & ( det->TYuNo );
    std::vector<Int_t>* TYuRing= & ( det->TYuRing );
    std::vector<Double_t>* TYuTheta= & ( det->TYuTheta ); // Yd theta angle

    Int_t* TSurMul= & ( det->TSurMul );
    std::vector<Int_t>* TSurChannel= & ( det->TSurChannel );
    std::vector<Double_t>* TSurEnergy= & ( det->TSurEnergy );
    std::vector<Int_t>* TSurADC= & ( det->TSurADC );
    Double_t* TSurEnergyCal= & ( det->TSurEnergyCal );

    Int_t* TSusMul= & ( det->TSusMul );
    std::vector<Int_t>* TSusChannel= & ( det->TSusChannel );
    std::vector<Double_t>* TSusEnergy= & ( det->TSusEnergy );
    std::vector<Int_t>* TSusADC= & ( det->TSusADC );
    std::vector<Double_t>* TSuTheta= & ( det->TSuTheta );
    std::vector<Double_t>* TSuPhi= & ( det->TSuPhi );

    tr_in->SetBranchAddress ( "det", &det );

    //output variables
    //YYD detector
    int YdMul;
    vector<int> YdChannel;
    vector<double> YdEnergyRaw;
    vector<double> YdEnergy;
    vector<int> YdRing;
    vector<int> YdSector;
  
    //YYU detector
    int YuMul;
    vector<int> YuChannel;
    vector<double> YuEnergyRaw;
    vector<double> YuEnergy;
    vector<int> YuRing;
    vector<int> YuSector;
    
    //Purpose of two readouts is to increase the range of operation (multiply it by 2 in total)
    //CsI 1 detector
    int CsI1Mul;
    vector<int> CsI1Channel; //Channel corresponds to one christal - 1 means that the gain is set to one value
    vector<double> CsI1EnergyRaw;

    //CsI 2 detector
    int CsI2Mul;
    vector<int> CsI2Channel;//Channel corresponds to one christal - 2 means that the gain is set to a different value than 1
    vector<double> CsI2EnergyRaw;

    //IC chamber
    vector<int> ICChannel;
    vector<double> ICEnergyRaw;

    //Sdd1 detector
    int Sd1rMul;
    vector<int> Sd1rChannel; //Ring Channels
    vector<double> Sd1rEnergyRaw;
    vector<double> Sd1rEnergy;
 
    //Sd1s detector
    int Sd1sMul;
    vector<int> Sd1sChannel ; //Sector Channels
    vector<double> Sd1sEnergyRaw;
    vector<double> Sd1sEnergy;
 
    //Sd2r detector
    int Sd2rMul;
    vector<int> Sd2rChannel; //Ring Channels
    vector<double> Sd2rEnergyRaw;
    vector<double> Sd2rEnergy;
 
    //Sd2s detector
    int Sd2sMul;
    vector<int> Sd2sChannel; //Sector Channels
    vector<double> Sd2sEnergyRaw;
    vector<double> Sd2sEnergy;

    //Sur
    int SurMul;
    vector<int> SurChannel; //Ring Channels
    vector<double> SurEnergyRaw;

    //Sus
    int SusMul;
    vector<int> SusChannel; //Sector Channels
    vector<double> SusEnergyRaw;

    //open the output file
 
    //string fOut_name = Form ( "/home/jerome/12Be_exp/Processed_files/Yd_calibration/rootindpndt_loss_%i_%d.root",run_num,index ); 
    string fOut_name = Form ( "/home/jerome/12Be_exp/Processed_files/Yd_calibration/rootindpndt_loss_decodetest_%i.root",run_num); 
    TFile* f_out = new TFile ( fOut_name.c_str(),"RECREATE" );
    TTree* tr_out = new TTree ( "AutoTree","AutoTree" );

    tr_out->Branch ( "YdMul",&YdMul );
    tr_out->Branch ( "YdChannel",&YdChannel );
    tr_out->Branch ( "YdEnergyRaw",&YdEnergyRaw );
    tr_out->Branch ( "YdEnergy",&YdEnergy );
    tr_out->Branch ( "YdRing",&YdRing );
    tr_out->Branch ( "YdSector",&YdSector );

    tr_out->Branch ( "YuMul",&YuMul );
    tr_out->Branch ( "YuChannel",&YuChannel );
    tr_out->Branch ( "YuEnergyRaw",&YuEnergyRaw );
    tr_out->Branch ( "YuEnergy",&YuEnergy );
    tr_out->Branch ( "YuRing",&YuRing );
    tr_out->Branch ( "YuSector",&YuSector ); 
   
    tr_out->Branch ( "CsI1Mul",&CsI1Mul );
    tr_out->Branch ( "CsI1Channel",&CsI1Channel );
    tr_out->Branch ( "CsI1EnergyRaw",&CsI1EnergyRaw );

    tr_out->Branch ( "CsI2Mul",&CsI2Mul );
    tr_out->Branch ( "CsI2Channel",&CsI2Channel );
    tr_out->Branch ( "CsI2EnergyRaw",&CsI2EnergyRaw );
   
    tr_out->Branch ( "ICChannel",&ICChannel );
    tr_out->Branch ( "ICEnergyRaw",&ICEnergyRaw );

    tr_out->Branch ( "Sd1rMul",&Sd1rMul );
    tr_out->Branch ( "Sd1rChannel",&Sd1rChannel );
    tr_out->Branch ( "Sd1rEnergyRaw",&Sd1rEnergyRaw );
    tr_out->Branch ( "Sd1rEnergy",&Sd1rEnergy );

    tr_out->Branch ( "Sd1sMul",&Sd1sMul );
    tr_out->Branch ( "Sd1sChannel",&Sd1sChannel );
    tr_out->Branch ( "Sd1sEnergyRaw",&Sd1sEnergyRaw );
    tr_out->Branch ( "Sd1sEnergy",&Sd1sEnergy );

    tr_out->Branch ( "Sd2rMul",&Sd2rMul );
    tr_out->Branch ( "Sd2rChannel",&Sd2rChannel );
    tr_out->Branch ( "Sd2rEnergyRaw",&Sd2rEnergyRaw );
    tr_out->Branch ( "Sd2rEnergy",&Sd2rEnergy );

    tr_out->Branch ( "Sd2sMul",&Sd2sMul );
    tr_out->Branch ( "Sd2sChannel",&Sd2sChannel );
    tr_out->Branch ( "Sd2sEnergyRaw",&Sd2sEnergyRaw );
    tr_out->Branch ( "Sd2sEnergy",&Sd2sEnergy );
 
    tr_out->Branch ( "SurMul",&SurMul );
    tr_out->Branch ( "SurChannel",&SurChannel );
    tr_out->Branch ( "SurEnergyRaw",&SurEnergyRaw );
  
    tr_out->Branch ( "SusMul",&SusMul );
    tr_out->Branch ( "SusChannel",&SusChannel );
    tr_out->Branch ( "SusEnergyRaw",&SusEnergyRaw );

    int ev = tr_in->GetEntries();
    cout << "Total number of events =" << ev << endl;
    int ev_num=0;

    //open the calibration files
    double stripYd[128], Ydp0[128],Ydp1[128];

    ifstream Alpha_peaksYd;
   // Alpha_peaksYd.open ( "/mnt/i/12Be_S1506/scripts/CalFiles/Yd1_CalPar.txt" );
    Alpha_peaksYd.open ( "/home/jerome/12Be_exp/scripts/YdeCalPar_3peaks.txt" );


    if ( Alpha_peaksYd.is_open() )
    {
        for ( int i=0; i<128; i++ )
        {
            Alpha_peaksYd >> stripYd[i] >> Ydp0[i] >> Ydp1[i] ;
            // cout << " Strip " <<  stripYd[i] << " P0 " << Ydp0[i] << " p1 " << Ydp1[i] << endl;
        }
    }

    else
    {
        cout << "No calib file for Yd ditector "<< endl;
    }
    Alpha_peaksYd.close();

    //Yu detector
    double stripYu[128], Yup0[128],Yup1[128], Yup2[128];

    ifstream Alpha_peaksYu;
    Alpha_peaksYu.open ( "/home/jerome/12Be_exp/scripts/Yu5225_pedestal_loss_AllPar.txt" ); //  YuCalPar5226_3peak_nolin.txt for carbon calibration run 5226
    //YuCalPar5225_3peaks.txt Be calibration

    if ( Alpha_peaksYu.is_open() )
    {
        for ( int i=0; i<128; i++ )
        {
            Alpha_peaksYu >> stripYu[i] >> Yup0[i] >> Yup1[i] ;
            // cout << " Strip " <<  stripYu[i] << " P0 " << Yup0[i] << " p1 " << Yup1[i] << endl;

        }
    }

    else
    {
        cout << " No calib file for Yu detector "<< endl;
    }
    Alpha_peaksYu.close();
    
    //Sd1 detector
    double stripSd1r[24], stripSd1s[32], Sd1rp0[24], Sd1rp1[24], Sd1sp0[32], Sd1sp1[32];
    double stripSd2r[24], stripSd2s[32], Sd2rp0[24], Sd2rp1[24], Sd2sp0[32], Sd2sp1[32];
    
    ifstream Alpha_Sd1r, Alpha_Sd1s;
    ifstream Alpha_Sd2r, Alpha_Sd2s;
    
    Alpha_Sd1r.open("/home/jerome/12Be_exp/scripts/Sd1r_4815Ped.txt");
    if(Alpha_Sd1r.is_open())
    {
     for(int i=0; i<24; i++)
     {
       Alpha_Sd1r >> stripSd1r[i] >> Sd1rp0[i] >> Sd1rp1[i] ;
     }
    }
    else {cout << "No calib file for Sd1r "<< endl;}

    Alpha_Sd1s.open ( "/mnt/i/12Be_S1506/scripts/Calibration/Calib_files/Sd1s_3peaks.txt" );
    if(Alpha_Sd1s.is_open())
    {
     for(int i=0; i<32; i++)
     {
       Alpha_Sd1s >> stripSd1s[i] >> Sd1sp0[i] >> Sd1sp1[i] ;
     }
    }
    else {cout << "No calib file for Sd1s "<< endl;}
    
    Alpha_Sd2r.open("/mnt/i/12Be_S1506/scripts/Calibration/Calib_files/Sd2r_3peaks.txt");
    if(Alpha_Sd2r.is_open())
    {
     for(int i=0; i<24; i++)
     {
       Alpha_Sd2r >> stripSd2r[i] >> Sd2rp0[i] >> Sd2rp1[i] ;
     }
    }
    else {cout << "No calib file for Sd2r "<< endl;}

    Alpha_Sd2s.open ( "/mnt/i/12Be_S1506/scripts/Calibration/Calib_files/Sd2s_3peaks.txt" );
    if(Alpha_Sd2s.is_open())
    {
     for(int i=0; i<32; i++)
     {
       Alpha_Sd2s >> stripSd2s[i] >> Sd2sp0[i] >> Sd2sp1[i] ;
     }
    }
    else {cout << "No calib file for Sd2s "<< endl;}

    while ( ev_num<=ev )
    {
        if ( ev_num%10000==0 ) cout << "Current event = " << ev_num << "\r"<< flush;
        tr_in->GetEntry ( ev_num );

        //reset the variables
        YdMul=-10;
        YuMul=-10;
        CsI1Mul=-10;
        CsI2Mul=-10;
        Sd1rMul=-10;
        Sd1sMul=-10;
        Sd2rMul=-10;
        Sd2sMul=-10;
        SurMul=-10;
        SusMul=-10;

	
	 if(det->TYdADC.size()>0)
	{
	  YdMul=0;
	  
	  //cout << "Size " << det->TYdADC.size() << endl;
	   for ( unsigned int adci=0; adci < det->TYdADC.size(); adci++ )
             {
	       if ( det->TYdADC.at ( adci ) >0 ){YdMul++; }//cout <<  " ALL channel " <<adci << " energy " << det->TYdADC.at ( adci ) << endl;}
	       
	     }
	     
	  double Ydmax = *max_element(det->TYdADC.begin(),det->TYdADC.end());
          int Ydchannel = distance(det->TYdADC.begin(),max_element(det->TYdADC.begin(),det->TYdADC.end()));
	  
	if(Ydmax>0)   
	  
	{YdEnergyRaw.push_back (Ydmax );
           YdChannel.push_back ( Ydchannel );
	   
	 double tempYdE = Ydmax *Ydp0[Ydchannel]+Ydp1[Ydchannel];
         YdEnergy.push_back ( tempYdE/1000 );
	  
	  //cout << "Good channel " << Ydchannel << " Energy " << Ydmax << endl;
	   
	    if ( Ydchannel<16 ){YdSector.push_back ( 0 ); YdRing.push_back ( Ydchannel );} //sector 0
            else if ( Ydchannel>=16 && Ydchannel<32 ){YdSector.push_back ( 1 ); YdRing.push_back ( Ydchannel-16 );} //sector 1
            else if ( Ydchannel>=32 && Ydchannel<48 ){YdSector.push_back ( 2 ); YdRing.push_back ( Ydchannel-32 );} //sector 2
            else if ( Ydchannel>=48 && Ydchannel<64 ){YdSector.push_back ( 3 ); YdRing.push_back ( Ydchannel-48 );} //sector 3
            else if ( Ydchannel>=64 && Ydchannel<80 ){YdSector.push_back ( 4 ); YdRing.push_back ( Ydchannel-64 );} //sector 4
            else if ( Ydchannel>=80 && Ydchannel<96 ){YdSector.push_back ( 5 ); YdRing.push_back ( Ydchannel-80 );} //sector 5
            else if ( Ydchannel>=96 && Ydchannel<112){YdSector.push_back ( 6 ); YdRing.push_back ( Ydchannel-96 );} //sector 6
            else if ( Ydchannel>=112 && Ydchannel<128){YdSector.push_back ( 7 );YdRing.push_back ( Ydchannel-112 );} //sector 7
	}
	}
      
       

        if ( TCsI1Channel->size() >0 )
        {
            //dE layermedDet
             CsI1Mul=0;
            for ( unsigned int adci=0; adci < det->TCsI1ADC.size(); adci++ )
            {
                if ( det->TCsI1ADC.at ( adci ) >0 )
                {
                    CsI1EnergyRaw.push_back ( TCsI1ADC->at ( adci ) );
                    CsI1Channel.push_back ( adci );
                    CsI1Mul++;
                }
            }
        }

        if ( TCsI2Channel->size() >0 )
        {
            //dE layer
            CsI2Mul =0;
            for ( unsigned int adci=0; adci < det->TCsI2ADC.size(); adci++ )
            {
                if ( det->TCsI2ADC.at ( adci ) >0 )
                {
                    CsI2EnergyRaw.push_back ( TCsI2ADC->at ( adci ) );
                    CsI2Channel.push_back ( adci );
                    CsI2Mul++;
                }
            }
        }

        //different approach to Yu fill
        if(det->TYuADC.size()>0)
	{
	  YuMul=0;
	   for ( unsigned int adci=0; adci < det->TYuADC.size(); adci++ )
             {
	       if ( det->TYuADC.at ( adci ) >0 ){YuMul++;}
	     }
	     
	      double max = *max_element(det->TYuADC.begin(),det->TYuADC.end());
          int channel = distance(det->TYuADC.begin(),max_element(det->TYuADC.begin(),det->TYuADC.end()));
	  
	  if(max>0)
	  {
	   YuEnergyRaw.push_back (max );
           YuChannel.push_back ( channel );
	   
	  double tempYuE = max *Yup0[channel]+Yup1[channel];
          YuEnergy.push_back ( tempYuE/1000 );
	   
	    if ( channel<16 ){YuSector.push_back ( 0 ); YuRing.push_back ( channel );} //sector 0
            else if ( channel>=16 && channel<32 ){YuSector.push_back ( 1 ); YuRing.push_back ( channel-16 );} //sector 1
            else if ( channel>=32 && channel<48 ){YuSector.push_back ( 2 ); YuRing.push_back ( channel-32 );} //sector 2
            else if ( channel>=48 && channel<64 ){YuSector.push_back ( 3 ); YuRing.push_back ( channel-48 );} //sector 3
            else if ( channel>=64 && channel<80 ){YuSector.push_back ( 4 ); YuRing.push_back ( channel-64 );} //sector 4
            else if ( channel>=80 && channel<96 ){YuSector.push_back ( 5 ); YuRing.push_back ( channel-80 );} //sector 5
            else if ( channel>=96 && channel<112){YuSector.push_back ( 6 ); YuRing.push_back ( channel-96 );} //sector 6
            else if ( channel>=112 && channel<128){YuSector.push_back ( 7 );YuRing.push_back ( channel-112 );} //sector 7
	  }
	}
    
    

        if ( TICChannel->size() >0 )
        {
            for ( unsigned int adci=0; adci < det->TICADC.size(); adci++ )
            {
                if ( det->TICADC.at ( adci ) >0 )
                {
                    ICEnergyRaw.push_back ( TICADC->at ( adci ) );
                    ICChannel.push_back ( adci );
                }
            }
        }

        if ( TSd1rChannel->size() >0 )
        {
            //The ring side of the S3 detector facing the beam -dE layer
            Sd1rMul=0;
            for ( unsigned int adci=0; adci < det->TSd1rADC.size(); adci++ )
            {
                if ( det->TSd1rADC.at ( adci ) >0 )
                {
                    Sd1rEnergyRaw.push_back ( TSd1rADC->at ( adci ) );
                    Sd1rChannel.push_back ( adci );
                    Sd1rMul++;
    
    if(TSd1rADC->at ( adci )>50)
    {
                    double tempSd1rE = TSd1rADC->at ( adci ) *Sd1rp0[adci]+Sd1rp1[adci];
                    Sd1rEnergy.push_back ( tempSd1rE/1000. );
    }
                }
            }
        }

        //the sector side of the S3 detector not facing the beam - dE layer
        if ( TSd1sChannel->size() >0 )
        {
           Sd1sMul=0;
            for ( unsigned int adci=0; adci < det->TSd1sADC.size(); adci++ )
            {
                if ( det->TSd1sADC.at ( adci ) >0 )
                {
                    Sd1sEnergyRaw.push_back ( TSd1sADC->at ( adci ) );
                    Sd1sChannel.push_back ( adci );
                    Sd1sMul++;
    
     if(TSd1sADC->at ( adci )>0)
    {
                    double tempSd1sE = TSd1sADC->at ( adci ) *Sd1sp0[adci]+Sd1sp1[adci];
                    Sd1sEnergy.push_back ( tempSd1sE/1000. );
    }
                }
            }
        }

        //The ring side of the S3 detector not facing the beam - E layer
        if ( TSd2rChannel->size() >0 )
        {
           Sd2rMul=0;
            for ( unsigned int adci=0; adci < det->TSd2rADC.size(); adci++ )
            {
                if ( det->TSd2rADC.at ( adci ) >0 )
                {
                    Sd2rEnergyRaw.push_back ( TSd2rADC->at ( adci ) );
                    Sd2rChannel.push_back ( adci );
                    Sd2rMul++;
    
     if(TSd2rADC->at ( adci )>0)
    {
                    double tempSd2rE = TSd2rADC->at ( adci ) *Sd2rp0[adci]+Sd2rp1[adci];
                    Sd2rEnergy.push_back ( tempSd2rE/1000. );
    }
                }
            }
        }

        //The sector side of the S3 detector facing the beam - E layer
        if ( TSd2sChannel->size() >0 )
        {
           Sd2sMul=0;
            for ( unsigned int adci=0; adci < det->TSd2sADC.size(); adci++ )
            {
                if ( det->TSd2sADC.at ( adci ) >0 )
                {
                    Sd2sEnergyRaw.push_back ( TSd2sADC->at ( adci ) );
                    Sd2sChannel.push_back ( adci );
                    Sd2sMul++;
    
     if(TSd2sADC->at ( adci )>50)
    {
                    double tempSd2sE = TSd2sADC->at ( adci ) *Sd2sp0[adci]+Sd2sp1[adci];
                    Sd2sEnergy.push_back ( tempSd2sE/1000. );
    }
                }
            }
        }

        if ( TSurChannel->size() >0 )
        {
            //ring side of the upstream S3 detector facing the target
            SurMul=0;
            for ( unsigned int adci=0; adci < det->TSurADC.size(); adci++ )
            {
                if ( det->TSurADC.at ( adci ) >0 )
                {
                    SurEnergyRaw.push_back ( TSurADC->at ( adci ) );
                    SurChannel.push_back ( adci );
                    SurMul++;
                }
            }
        }

        //sector side of the upstream S3 detector not facing the target
        if ( TSusChannel->size() >0 )
        {
          SusMul=0;
            for ( unsigned int adci=0; adci < det->TSusADC.size(); adci++ )
            {
                if ( det->TSusADC.at ( adci ) >0 )
                {
                    SusEnergyRaw.push_back ( TSusADC->at ( adci ) );
                    SusChannel.push_back ( adci );
                    SusMul++;
                }
            }
        }


        tr_out->Fill();

        YdChannel.clear();
        YdEnergyRaw.clear();
        YdEnergy.clear();
        YdRing.clear();
        YdSector.clear();
      
        YuChannel.clear();
        YuEnergyRaw.clear();
        YuEnergy.clear();
        YuRing.clear();
        YuSector.clear();
       
        CsI1Channel.clear();
        CsI1EnergyRaw.clear();
        CsI2Channel.clear();
        CsI2EnergyRaw.clear();
        ICChannel.clear();
        ICEnergyRaw.clear();

        Sd1rChannel.clear();
        Sd1rEnergyRaw.clear();
        Sd1rEnergy.clear();
        Sd1sChannel .clear();
        Sd1sEnergyRaw.clear();
        Sd1sEnergy.clear();

        Sd2rChannel.clear();
        Sd2rEnergyRaw.clear();
	Sd2rEnergy.clear();
        Sd2sChannel.clear();
        Sd2sEnergyRaw.clear();
	Sd2sEnergy.clear();

        SurChannel.clear();
        SurEnergyRaw.clear();
        SusChannel.clear();
        SusEnergyRaw.clear();

        ev_num++;
    } //end of the main while loop

    f_out->cd();
    tr_out->Write();
    f_out->Close();

    return 0;


} //end of the script
