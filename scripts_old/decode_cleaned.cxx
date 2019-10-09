//M. Vostinar
using namespace std;

#include <TFile.h>
#include "TString.h"
#include "string"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "iostream"
#include "TMath.h"
#include "/home/maja/Software/treeIris/include/IDet.h" //the treeIris library dependancy
#include <vector>
#include <TObject.h>
#include <TClass.h>
#include "sstream"
#include "fstream"
#include "cstdlib"
#include <algorithm> 

double decode ( int run_num )
{
      
    string f_name = Form ( "/mnt/f/12Be_S1506/RootFiles/NewFiles/C/tree%i.root",run_num ); //open the raw root file you want to convert to the format used by my codes
    TFile* f_in = TFile::Open ( f_name.c_str(),"READ" );
    TTree* tr_in = ( TTree* ) f_in->Get ( "Iris" ); //reading the Iris tree; 
    //only the ADC leafs of the tree are used in the code. The rest of the data has a vector-array configuration or something like that and presents reliability which introduces unnecessary uncertanties
    
    IDet* det = new IDet();

    Int_t* TYdMul = & ( det->TYdMul );
    std::vector<Int_t>* TYdChannel = & ( det->TYdChannel );
    std::vector<Double_t>* TYdEnergy= & ( det->TYdEnergy );
    std::vector<Int_t>* TYdADC= & ( det->TYdADC ); 
    std::vector<Int_t>* TYdNo= & ( det->TYdNo );
    std::vector<Int_t>* TYdRing= & ( det->TYdRing );
    std::vector<Double_t>* TYdTheta= & ( det->TYdTheta ); 

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
 
    string fOut_name = Form ( "/mnt/f/12Be_S1506/Processed_rootFiles/NewFiles/Carbon/C_YuEnergyLoss_%i.root",run_num ); //Open the output file
    TFile* f_out = new TFile ( fOut_name.c_str(),"RECREATE" );
    TTree* tr_out = new TTree ( "AutoTree","AutoTree" ); //create the new root tree
    //Branches created in the new file; the leafs are removed from the three as well as arrays of vectors or whatever that was
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

    int ev = tr_in->GetEntries(); //get the total number of events in the tree
    cout << "Total number of events =" << ev << endl;
    int ev_num=0; 

    //open the calibration files
    //Opening the Yd calibration file, reading it and assigning the relevant variables to the internaly defined variables
    double stripYd[128], Ydp0[128],Ydp1[128]; //definition of gain and offset used for the Yd calibration

    ifstream Alpha_peaksYd;
    Alpha_peaksYd.open ( "/mnt/f/12Be_S1506/scripts/Calibration/Calib_files/YdCalPar_3peaks.txt" ); //Open the calibration file


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
        cout << "No calib file for Yd ditector "<< endl;;
    }
    Alpha_peaksYd.close();

    //Opening and reading the Yu calibration file, same as for Yd
    double stripYu[128], Yup0[128],Yup1[128], Yup2[128];

    ifstream Alpha_peaksYu;
    Alpha_peaksYu.open ( "/mnt/f/12Be_S1506/scripts/Yu5225_pedestal_loss_AllPar_jerome.txt" ); //  YuCalPar5226_3peak_nolin.txt for carbon calibration run 5226
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
    
    //Openning and readin the calibration files for the Sd1 and Sd2 detectors
    Alpha_Sd1r.open("/home/jerome/12Be_exp/scripts/Sd1r_4815Ped.txt");
    if(Alpha_Sd1r.is_open())
    {
     for(int i=0; i<24; i++)
     {
       Alpha_Sd1r >> stripSd1r[i] >> Sd1rp0[i] >> Sd1rp1[i] ;
     }
    }
    else {cout << "No calib file for Sd1r "<< endl;}

    Alpha_Sd1s.open ( "/mnt/f/12Be_S1506/scripts/Calibration/Calib_files/Sd1s_3peaks.txt" );
    if(Alpha_Sd1s.is_open())
    {
     for(int i=0; i<32; i++)
     {
       Alpha_Sd1s >> stripSd1s[i] >> Sd1sp0[i] >> Sd1sp1[i] ;
     }
    }
    else {cout << "No calib file for Sd1s "<< endl;}
    
    Alpha_Sd2r.open("/mnt/f/12Be_S1506/scripts/Calibration/Calib_files/Sd2r_3peaks.txt");
    if(Alpha_Sd2r.is_open())
    {
     for(int i=0; i<24; i++)
     {
       Alpha_Sd2r >> stripSd2r[i] >> Sd2rp0[i] >> Sd2rp1[i] ;
     }
    }
    else {cout << "No calib file for Sd2r "<< endl;}

    Alpha_Sd2s.open ( "/mnt/f/12Be_S1506/scripts/Calibration/Calib_files/Sd2s_3peaks.txt" );
    if(Alpha_Sd2s.is_open())
    {
     for(int i=0; i<32; i++)
     {
       Alpha_Sd2s >> stripSd2s[i] >> Sd2sp0[i] >> Sd2sp1[i] ;
     }
    }
    else {cout << "No calib file for Sd2s "<< endl;}
    //end of calibration filles

    while ( ev_num<=ev ) //entering the main while loop and looping over all the events in the Iris tree
    {
        if ( ev_num%10000==0 ) cout << "Current event = " << ev_num << "\r"<< flush;
        tr_in->GetEntry ( ev_num ); //pulling up the current event from the Iris tree

        //reset the variables; in this case set the multiplicity to 0, to get the correct multiplicity in the detectors for each event
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

	
	 if(det->TYdADC.size()>0) //treating the Yd detector
	{
	  YdMul=0; //set the multiplicity to 0
	  
	  //cout << "Size " << det->TYdADC.size() << endl;
	   for ( unsigned int adci=0; adci < det->TYdADC.size(); adci++ ) //loop over the size of the YdADC vector, this is always = 128
             {
	       if ( det->TYdADC.at ( adci ) >0 ){YdMul++; }//calculate multiplicity; count every entry >0 for that event
	       
	     }
	     
	  double Ydmax = *max_element(det->TYdADC.begin(),det->TYdADC.end()); //find the maximum value in the YdADC vector; this should correspond to the acctual hit in the strip
          int Ydchannel = distance(det->TYdADC.begin(),max_element(det->TYdADC.begin(),det->TYdADC.end())); //Get the position of the max value. This corresponds to the strip of the detector which fired and is directly related to the position in the vector
	  
	if(Ydmax>0) //if the maximal value in the vector was >0 (IT can be 0) fill in the relevant information 
	  
	{YdEnergyRaw.push_back (Ydmax ); //fill the raw energy in a vector
           YdChannel.push_back ( Ydchannel ); //fill the strip which has firred
	   
	  double tempYdE = Ydmax *Ydp0[Ydchannel]+Ydp1[Ydchannel]; //calibrate the energy using the calibration parameters loaded earlier
          YdEnergy.push_back ( tempYdE/1000 ); //convert the energy from keV to MeV
	  
	  //cout << "Good channel " << Ydchannel << " Energy " << Ydmax << endl;
	   //determine the Sector and the Ring which fired and fill it in the corresponding vector in the tree
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
      
       //treating the CsI detector, the same way as it was done for the Yd detector
       //The signals from th CsI detector were split in two and send to two different gain settings, so the channels for the CsI are multiplied
        if ( TCsI1Channel->size() >0 )
        {
            //E layerDet
             CsI1Mul=0;
            for ( unsigned int adci=0; adci < det->TCsI1ADC.size(); adci++ )
            {
                if ( det->TCsI1ADC.at ( adci ) >0 ){ CsI1Mul++;}
            }
            
            double CsImax = *max_element(det->TCsI1ADC.begin(),det->TCsI1ADC.end());
	    int CsIchannel = distance(det->TCsI1ADC.begin(),max_element(det->TCsI1ADC.begin(),det->TCsI1ADC.end()));
	    
	   if(CsImax>0) { CsI1EnergyRaw.push_back ( CsImax);
            CsI1Channel.push_back (CsIchannel);}
            
        }

        if ( TCsI2Channel->size() >0 )
        {
            //E layer
            CsI2Mul =0;
            for ( unsigned int adci=0; adci < det->TCsI2ADC.size(); adci++ )
            {
                if ( det->TCsI2ADC.at ( adci ) >0 ){CsI2Mul++;}
            }
            
            double CsI2max = *max_element(det->TCsI2ADC.begin(),det->TCsI2ADC.end());
	    int CsI2channel = distance(det->TCsI2ADC.begin(),max_element(det->TCsI2ADC.begin(),det->TCsI2ADC.end()));
	      
	    if(CsI2max>0) { CsI2EnergyRaw.push_back ( CsI2max );
			    CsI2Channel.push_back (CsI2channel );}
        }

        //treat the Yu detector, same as all the previous ones
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
          YuEnergy.push_back ( tempYuE/1000 ); //convert from keV to MeV in the calibrated energy of the Yu detector
	   
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
        
        //treatment of the IC chamber
        //The IC chamber has only one channel, it is plugged in channel 15
        if ( TICChannel->size() >0 )
        {
            for ( unsigned int adci=0; adci < det->TICADC.size(); adci++ )
            {
	      
	       double maxIC = *max_element(det->TICADC.begin(),det->TICADC.end()); // find the maximum value in the Yu vector and asign it to maxYu
	       int channelIC = distance(det->TICADC.begin(),max_element(det->TICADC.begin(),det->TICADC.end()));
	      
                if ( maxIC >0 )
                {
                    ICEnergyRaw.push_back ( maxIC );
                    ICChannel.push_back ( channelIC );
                }
            }
        }
        
        //Treating the Sd1r part of the Sd1 detector. There are 2 Sd detectors. Sd1 is the one making the dE layer in the dE-E telechope. It has 24 rings and 32 sectors. The rings are oriented toward the beam. 
        //The data is gathered in a way that the sectors and rings are registered separately which gives the rings as Sd1r. The Sd2 detector is treated the same as Sd1. 
        if ( det->TSd1rADC.size() >0 )
        {
            Sd1rMul=0;
            for ( unsigned int adci=0; adci < det->TSd1rADC.size(); adci++ )
            {
                if ( det->TSd1rADC.at ( adci ) >0 ){Sd1rMul++;}
                
                    double maxSd1r = *max_element(det->TSd1rADC.begin(),det->TSd1rADC.end());
                    int channelSd1r = distance(det->TSd1rADC.begin(),max_element(det->TSd1rADC.begin(),det->TSd1rADC.end()));
                
		if(max>0){
		  Sd1rEnergyRaw.push_back ( maxSd1r);
		  Sd1rChannel.push_back ( channelSd1r );
		      
		      double tempSd1rE = maxSd1r *Sd1rp0[adci]+Sd1rp1[adci];
		      Sd1rEnergy.push_back ( tempSd1rE/1000. );
		}
		  
            }
        }

        //the sector side of the Sd3 detector not facing the beam - dE layer
        if ( det->TSd1sADC.size() >0 )
        {
           Sd1sMul=0;
            for ( unsigned int adci=0; adci < det->TSd1sADC.size(); adci++ )
            {
                if ( det->TSd1sADC.at ( adci ) >0 ){Sd1sMul++;}
                
                
                double maxSd1s = *max_element(det->TSd1sADC.begin(),det->TSd1sADC.end());
                int channelSd1s = distance(det->TSd1sADC.begin(),max_element(det->TSd1sADC.begin(),det->TSd1sADC.end()));
                
               if(maxSd1s>0){
		 Sd1sEnergyRaw.push_back ( maxSd1s );
                 Sd1sChannel.push_back ( channelSd1s );
                      
                 double tempSd1sE = maxSd1s*Sd1sp0[adci]+Sd1sp1[adci];
                 Sd1sEnergy.push_back ( tempSd1sE/1000. );
	       }
                
            }
        }

        //The Sd2 detector faces the beem with its sector side; This is the treatment of the ring side of the detector
        if ( det->TSd2rADC.size() >0 )
        {
           Sd2rMul=0;
            for ( unsigned int adci=0; adci < det->TSd2rADC.size(); adci++ )
            {
                if ( det->TSd2rADC.at ( adci ) >0 ){Sd2rMul++;}
                
                double maxSd2r = *max_element(det->TSd2rADC.begin(),det->TSd2rADC.end());
                int channelSd2r = distance(det->TSd2rADC.begin(),max_element(det->TSd2rADC.begin(),det->TSd2rADC.end()));
                
		if(maxSd2r>0)
		{
                    Sd2rEnergyRaw.push_back ( maxSd2r );
                    Sd2rChannel.push_back ( channelSd2r );
     
                    double tempSd2rE =  maxSd2r*Sd2rp0[adci]+Sd2rp1[adci];
                    Sd2rEnergy.push_back ( tempSd2rE/1000. );
		}
               }
        }

        //The sector side of the S3 detector facing the beam - E layer
        if ( det->TSd2sADC.size() >0 )
        {
           Sd2sMul=0;
            for ( unsigned int adci=0; adci < det->TSd2sADC.size(); adci++ )
            {
                if ( det->TSd2sADC.at ( adci ) >0 ){Sd2sMul++;}
                
                double maxSd2s = *max_element(det->TSd2sADC.begin(),det->TSd2sADC.end());
                int channelSd2s = distance(det->TSd2sADC.begin(),max_element(det->TSd2sADC.begin(),det->TSd2sADC.end()));
                
               if(maxSd2s>0){
		
		 Sd2sEnergyRaw.push_back ( maxSd2s);
                 Sd2sChannel.push_back ( channelSd2s );
 
                 double tempSd2sE = maxSd2s*Sd2sp0[adci]+Sd2sp1[adci];
                 Sd2sEnergy.push_back ( tempSd2sE/1000. );
	      }
                
            }
        }

        //another S3 detector positioned upstream in front of the Yu detector
        //The ring side of this detector is facing the target and its NOT facing the beam
        if ( det->TSurADC.size() >0 )
        {
            //ring side of the upstream S3 detector facing the target
            SurMul=0;
            for ( unsigned int adci=0; adci < det->TSurADC.size(); adci++ )
            {
                if ( det->TSurADC.at ( adci ) >0 ){SurMul++;}
                
                double maxSur = *max_element(det->TSurADC.begin(),det->TSurADC.end());
                int channelSur = distance(det->TSurADC.begin(),max_element(det->TSurADC.begin(),det->TSurADC.end()));
                
                if(maxSur>0) {
		SurEnergyRaw.push_back ( maxSur);
                SurChannel.push_back ( channelSur );
		}
                
            }
        }

        //sector side of the upstream S3 detector is facing the beam
        if ( det->TSusADC.size() >0 )
        {
          SusMul=0;
            for ( unsigned int adci=0; adci < det->TSusADC.size(); adci++ )
            {
                if ( det->TSusADC.at ( adci ) >0 ){SusMul++;}
                
                double maxSus = *max_element(det->TSusADC.begin(),det->TSusADC.end());
                int channelSus = distance(det->TSusADC.begin(),max_element(det->TSusADC.begin(),det->TSusADC.end()));
                
		if(maxSus>0){
                    SusEnergyRaw.push_back ( maxSus);
                    SusChannel.push_back (channelSus );
		}
                
            }
        }


        tr_out->Fill(); //fill out the tree

	//clear the vectors so you do not multiply the data
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

        ev_num++; //looping over the while loop
    } //end of the main while loop

    f_out->cd(); //access the output file
    tr_out->Write(); //write the output file
    f_out->Close(); //close the output file

    return 0;


} //end of the script
