using namespace std;

#include <TFile.h>
#include "TString.h"
#include "string"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "iostream"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TCutG.h"
#include <vector>
#include <TObject.h>
#include <TClass.h>
#include "sstream"
#include "fstream"
#include "cstdlib"
#include "TRandom2.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TF1.h"
#include "TSpectrum.h"

double analysis()
{
    //open the output file
    TFile *f_out = new TFile("Be00_NEw18_3.root","RECREATE");
  
//Open the input files
    TChain *chain = new TChain ( "AutoTree" );
    
     //test for alphas
    //Check_5225
  /*  chain->Add("/mnt/i/12Be_S1506/Processed_rootFiles/AlphaNo0_5225_1.root");
    chain->Add("/mnt/i/12Be_S1506/Processed_rootFiles/AlphaNo0_5225_2.root");
    chain->Add("/mnt/i/12Be_S1506/Processed_rootFiles/AlphaNo0_5225_3.root");*/
  
/*  chain->Add("/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/AlphaCalLin5225_1.root"); 
  chain->Add("/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/AlphaCalLin5225_2.root"); 
  chain->Add("/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/AlphaCalLin5225_3.root"); */
    
  /*  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Alpha/AlphaAllPar_5225_1.root" );
    chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Alpha/AlphaAllPar_5225_2.root" );
    chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Alpha/AlphaAllPar_5225_3.root" );*/
    
    

    //Carbon data
  /*  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/C_data/C_5001.root" );
    chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/C_data/C_5002.root" );
    chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/C_data/C_5003.root" );
    chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/C_data/C_5004.root" );
    chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/C_data/C_5006.root" );
    chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/C_data/C_5007.root" );
    chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/C_data/C_5009.root" );*/
  
 /* chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5021.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5022.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5023.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5024.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5025.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5027.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5028.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5029.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5030.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5031.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5032.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5033.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5034.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5035.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5036.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5037.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5038.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5039.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5044.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5045.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5048.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5049.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5050.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5051.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5052.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5053.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5054.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5055.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5056.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5057.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5058.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5064.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5065.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5066.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5070.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5071.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5072.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5073.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5074.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5075.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5102.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5103.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5104.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5105.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5106.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5107.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5108.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5109.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5115.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5116.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5117.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5122.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5123.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5124.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5125.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5126.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5128.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5129.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5130.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5131.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5132.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5142.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5143.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5144.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5145.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5146.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5147.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5148.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5162.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5163.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5164.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5165.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5166.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5167.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5169.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5170.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5171.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5172.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5173.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5174.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5188.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5189.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5190.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5191.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5192.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5193.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5194.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5195.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5196.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5197.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5198.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5199.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5201.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5202.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5203.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5204.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5205.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5206.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5207.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5209.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5210.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5211.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5212.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5213.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5214.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5215.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5216.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5217.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5218.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5219.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5220.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/BeAllPed_5221.root");*/
 
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5021.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5022.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5023.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5024.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5025.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5027.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5028.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5029.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5030.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5031.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5032.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5033.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5034.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5035.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5036.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5037.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5038.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5039.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5044.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5045.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5048.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5049.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5050.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5051.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5052.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5053.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5054.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5055.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5056.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5057.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5058.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5064.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5065.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5066.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5070.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5071.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5072.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5073.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5074.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5075.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5102.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5103.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5104.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5105.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5106.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5107.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5108.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5109.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5115.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5116.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5117.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5122.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5123.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5124.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5125.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5126.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5128.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5129.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5130.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5131.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5132.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5142.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5143.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5144.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5145.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5146.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5147.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5148.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5162.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5163.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5164.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5165.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5166.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5167.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5169.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5170.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5171.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5172.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5173.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5174.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5188.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5189.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5190.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5191.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5192.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5193.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5194.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5195.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5196.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5197.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5198.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5199.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5201.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5202.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5203.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5204.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5205.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5206.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5207.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5209.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5210.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5211.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5212.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5213.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5214.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5215.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5216.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5217.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5218.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5219.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5220.root");
  chain->Add ( "/mnt/i/12Be_S1506/Processed_rootFiles/NewFiles/Be00_5221.root");
  
  
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
    vector<int>* ICChannel= new vector<int>();
    vector<double>* ICEnergyRaw= new vector<double>();

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
    geometry.open ( "/mnt/i/12Be_S1506/calib/geometry_s1506.txt" ); //open the calibration file
    if ( !geometry.is_open() )
    {
        cout << " No Geometry file found " << endl;
        return -1;
    }

    string read_geometry;
    istringstream iss;
    string name, dummy;
    float Ydr0,Ydr1,Ydz, Yuz, Sd1z, Sd2z, Sdr0, Sdr1; // Ydz distance from the target, Ydr0 inner radius, Ydz outer radius

    while ( getline ( geometry,read_geometry ) ) //in the calibration file named geometry start reading the lines
    {
        if ( read_geometry.find ( "YD_DISTANCE",0 ) !=string::npos ) //if you find the "YD_DISTANCE" before the end of the file
        {
            iss.clear();
            iss.str ( read_geometry );
            iss >> name >> dummy >> Ydz;

            // cout << " name " << name << " / " << dummy << " number " << Ydz << endl;
        }

        if ( read_geometry.find ( "YD_INNER_RADIUS",0 ) !=string::npos )
        {
            iss.clear();
            iss.str ( read_geometry );
            iss >> name >> dummy >> Ydr0;
        }

        if ( read_geometry.find ( "YD_OUTER_RADIUS",0 ) !=string::npos )
        {
            iss.clear();
            iss.str ( read_geometry );
            iss >> name >> dummy >> Ydr1;
        }

        if ( read_geometry.find ( "YU_DISTANCE",0 ) !=string::npos )
        {
            iss.clear();
            iss.str ( read_geometry );
            iss >> name >> dummy >> Yuz;
        }

        if ( read_geometry.find ( "SD1_DISTANCE",0 ) !=string::npos )
        {
            iss.clear();
            iss.str ( read_geometry );
            iss >> name >> dummy >> Sd1z;
        }

        if ( read_geometry.find ( "SD2_DISTANCE",0 ) !=string::npos )
        {
            iss.clear();
            iss.str ( read_geometry );
            iss >> name >> dummy >> Sd2z;
        }

        if ( read_geometry.find ( "SD_INNER_RADIUS",0 ) !=string::npos )
        {
            iss.clear();
            iss.str ( read_geometry );
            iss >> name >> dummy >> Sdr0;
        }

        if ( read_geometry.find ( "SD_OUTER_RADIUS",0 ) !=string::npos )
        {
            iss.clear();
            iss.str ( read_geometry );
            iss >> name >> dummy >> Sdr1;
        }
    }
    //end of the geometry file
    
    //load the geometrical cuts
    TFile *f_cut = TFile::Open("/mnt/i/12Be_S1506/scripts/PIDcut1r2r.root"); //PID cut for C
    TCutG *pidcut = (TCutG*) f_cut->Get("PID1r2r"); // for C
    
    TFile *f_cutBe1 = TFile::Open("/mnt/i/12Be_S1506/scripts/cutICSd1r2r.root");
    TCutG *cutIc1 = (TCutG*) f_cutBe1->Get("ICSd1r2r");
    
    TFile *f_cutBe2 = TFile::Open("/mnt/i/12Be_S1506/scripts/cutICSd1r2s.root");
    TCutG *cutIc2 = (TCutG*) f_cutBe2->Get("ICSd1r2s");
    
    TFile *f_cutBe3 = TFile::Open("/mnt/i/12Be_S1506/scripts/cutICSd1s2r.root");
    TCutG *cutIc3 = (TCutG*) f_cutBe3->Get("ICSd1s2r");
    
    TFile *f_cutBe4 = TFile::Open("/mnt/i/12Be_S1506/scripts/cutICSd1s2s.root");
    TCutG *cutIc4 = (TCutG*) f_cutBe4->Get("ICSd1s2s");
    
    TFile *f_cutBe5 = TFile::Open("/mnt/i/12Be_S1506/scripts/cutBe_clean.root"); //new cut on Be with PID stuff take into account
    TCutG *cutIc5 = (TCutG*) f_cutBe5->Get("Be");
    
    TFile *f_cutYd1 = TFile::Open("/mnt/i/12Be_S1506/scripts/cut1YdCsI.root");
    TCutG *cutYd1 = (TCutG*) f_cutYd1->Get("band1");
    
    TFile *f_cutYd2 = TFile::Open("/mnt/i/12Be_S1506/scripts/cut2YdCsI.root");
    TCutG *cutYd2 = (TCutG*) f_cutYd2->Get("bend2");
    
    TFile *f_cutYd3 = TFile::Open ( "/mnt/i/12Be_S1506/scripts/cut3YdCsI.root" );
    TCutG *cutYd3 = (TCutG*) f_cutYd3->Get("band3");
    
    TFile *f_cutBe6 = TFile::Open ( "/mnt/i/12Be_S1506/scripts/PIDWideCut.root" );
    TCutG *cutBe6 = (TCutG*) f_cutBe6->Get("PIDcut");
    
    TFile *f_cutYu = TFile::Open ( "/mnt/i/12Be_S1506/scripts/YuBecut.root" );
    TCutG *cutYu = (TCutG*) f_cutYu->Get("Be");
    
    TFile *f_cutBeNew = TFile::Open ( "/mnt/i/12Be_S1506/scripts/Becut1.root");
    TCutG *cutBeNew = (TCutG*) f_cutBeNew->Get("Becut1");
    
    TFile *f_cutBeSmall = TFile::Open ( "/mnt/i/12Be_S1506/scripts/SmallBePIDgate.root");
    TCutG *cutBeSmall = (TCutG*) f_cutBeSmall->Get("SmallBe");
    
   TFile *f_cutBeKine = TFile::Open ( "/mnt/i/12Be_S1506/scripts/BeKinematic.root");
   TCutG *cutBeKine = (TCutG*) f_cutBeKine->Get("CUTG");
    
   
    
    //definition of histograms
    //PID plotts
    TH2D *hSd1rSd2r = new TH2D("hSd1rSd2r","Sd1r vs Sd2r",2048,0,4096,2048,0,4096); //non calibrated Sd1r vs Sd2r
    TH2D *hSd1rSd2rCut = new TH2D("hSd1rSd2rCut","Sd1r vs Sd2r with correlation requirenment",2048,0,4096,2048,0,4096); //non calibrated Sd1r vs Sd2r
    TH2D *hSd1rSd2s = new TH2D("hSd1rSd2s","Sd1r vs Sd2r",2048,0,4096,2048,0,4096); //non calibrated Sd1r vs Sd2r
    TH2D *hSd1sSd2r = new TH2D("hSd1sSd2r","Sd1s vs Sd2r",2048,0,4096,2048,0,4096); //non calibrated Sd1s vs Sd2r
    TH2D *hSd1sSd2s = new TH2D("hSd1sSd2s","Sd1s vs Sd2r",2048,0,4096,2048,0,4096); //non calibrated Sd1s vs Sd2r
    
    TH2D *hCSd1rSd2r = new TH2D("hCSd1rSd2r","Calibrated Sd1r vs Sd2r",1000,0,100,250,0,25); // calibrated Sd1r vs Sd2r
    
    TH2D *hCSd1rSd2s = new TH2D("hCSd1rSd2s","Calibrated Sd1r vs Sd2s",1000,0,100,250,0,25); // calibrated Sd1r vs Sd2r
    
    TH2D *hCSd1sSd2r = new TH2D("hCSd1sSd2r","Calibrated Sd1s vs Sd2r",1000,0,100,250,0,25); // calibrated Sd1s vs Sd2r
    TH2D *hCSd1sSd2s = new TH2D("hCSd1sSd2s","Calibrated Sd1s vs Sd2s",1000,0,100,250,0,25); // calibrated Sd1s vs Sd2r
    
    TH2D *hCSd1rSd2rIC = new TH2D("hCSd1rSd2rIC","Cal Sd1r vs Sd2r, gated by IC",1000,0,100,250,0,25); // calibrated Sd1r vs Sd2r
    TH2D *hCSd1rSd2rICCut = new TH2D("hCSd1rSd2rICCut","Cal Sd1r vs Sd2r, gated by IC gated by the channel corr",1000,0,100,250,0,25); // calibrated Sd1r vs Sd2r
    TH2D *hCSd1rSd2rIC2 = new TH2D("hCSd1rSd2rIC2","Cal Sd1r vs Sd2r, gated by IC2",1000,0,100,250,0,25); // calibrated Sd1r vs Sd2r for second IC peak in Be data
    TH2D *hCSd1rSd2sIC = new TH2D("hCSd1rSd2sIC","Cal Sd1r vs Sd2s, gated by IC",1000,0,100,250,0,25); // calibrated Sd1r vs Sd2r
    TH2D *hCSd1sSd2rIC = new TH2D("hCSd1sSd2rIC","Cal Sd1s vs Sd2r, gated by IC",1000,0,100,250,0,25); // calibrated Sd1s vs Sd2r
    TH2D *hCSd1sSd2sIC = new TH2D("hCSd1sSd2sIC","Cal Sd1s vs Sd2s, gated by IC",1000,0,100,250,0,25); // calibrated Sd1s vs Sd2r
    TH2D *hCSd1rSd2rICYu = new TH2D("hCSd1rSd2rICYU","Cal Sd1r vs Sd2r gated by IC and Yu cut on Be",1000,0,100,250,0,25); // calibrated Sd1s vs Sd2r
    
    TH2D *hCSd1rSd2rYd = new TH2D("hCSd1rSd2rYd","Cal Sd1r vs Sd2r, gated by Yd",1000,0,100,250,0,25); // calibrated Sd1r vs Sd2r
    TH2D *hCSd1rSd2sYd = new TH2D("hCSd1rSd2sYd","Cal Sd1r vs Sd2s, gated by Yd",1000,0,100,250,0,25); // calibrated Sd1r vs Sd2r
    TH2D *hCSd1sSd2rYd = new TH2D("hCSd1sSd2rYd","Cal Sd1s vs Sd2r, gated by Yd",1000,0,100,250,0,25); // calibrated Sd1s vs Sd2r
    TH2D *hCSd1sSd2sYd = new TH2D("hCSd1sSd2sYd","Cal Sd1s vs Sd2s, gated by Yd",1000,0,100,250,0,25); // calibrated Sd1s vs Sd2r
    
    TH2D *hCSd1rSd2rYdIC = new TH2D("hCSd1rSd2rYdIC","Cal Sd1r vs Sd2r, gated by Yd & IC",1000,0,100,250,0,25); // calibrated Sd1r vs Sd2r
    TH2D *hCSd1rSd2rYdICCut = new TH2D("hCSd1rSd2rYdICCut","Cal Sd1r vs Sd2r, gated by Yd & IC & chann corr",1000,0,100,250,0,25); // calibrated Sd1r vs Sd2r
    TH2D *hCSd1rSd2sYdIC = new TH2D("hCSd1rSd2sYdIC","Cal Sd1r vs Sd2s, gated by Yd & IC",1000,0,100,250,0,25); // calibrated Sd1r vs Sd2r
    TH2D *hCSd1sSd2rYdIC = new TH2D("hCSd1sSd2rYdIC","Cal Sd1s vs Sd2r, gated by Yd & IC",1000,0,100,250,0,25); // calibrated Sd1s vs Sd2r
    TH2D *hCSd1sSd2sYdIC = new TH2D("hCSd1sSd2sYdIC","Cal Sd1s vs Sd2s, gated by Yd & IC",1000,0,100,250,0,25); // calibrated Sd1s vs Sd2r
    
    TH2D *hSd1r2rYdBand1 = new TH2D("hSd1r2rYdBand1","hSd1r2rYdBand1",1000,0,100,250,0,25);
    TH2D *hSd1r2rYdBand2 = new TH2D("hSd1r2rYdBand2","hSd1r2rYdBand2",1000,0,100,250,0,25);
    TH2D *hSd1r2rYdBand3 = new TH2D("hSd1r2rYdBand3","hSd1r2rYdBand3",1000,0,100,250,0,25);
    
    TH2D *hSdr1r2Ch = new TH2D("hSdr1r2Ch","Ring 1 vs ring 2",24,0,24,24,0,24);
    TH2D *hSdr1s2Ch = new TH2D("hSdr1s2Ch","Ring 1 vs sector 2",32,0,32,24,0,24);
    TH2D *hSds1r2Ch = new TH2D("hSds1r2Ch","Sector 1 vs ring 2",32,0,32,24,0,24);
    TH2D *hSds1s2Ch = new TH2D("hSds1s2Ch","Sector 1 vs sector 2",32,0,32,32,0,32);
    
    //IC non calibrated
    TH1D *hIC = new TH1D("hIC","IC raw energy",2048,0,4096);
    
    //Yd detector => elastics
    //calculate the bins for this one to correspond to angles we see
    double Ydbins[17]={0}, Ydbintemp[17]={0};
    for (int i=0; i<17;i++){Ydbintemp[i]=TMath::RadToDeg()*TMath::ATan((Ydr0+(16-i)*((Ydr1-Ydr0)/16))/Ydz);
    }
    
    for(int k=16; k>=0; k--)
    {
      int a = TMath::Abs(k-16);
      Ydbins[a]=Ydbintemp[k];
    }
    TH2D *hYdAn = new TH2D("hYdAn","YdE vs Angle", 16,Ydbins,10000,0,5); //no gates at the moment
    TH2D *hYdAnIC = new TH2D("hYdAnIC","YdE vs Angle with an IC gate", 16,Ydbins,10000,0,5); //gate on the IC peak
    TH2D *hYdAnIC2 = new TH2D("hYdAnIC2","YdE vs Angle with an IC2 gate", 16,Ydbins,10000,0,5); //gate on the 2nd IC peak in Be data
    TH2D *hYdAnPID = new TH2D("hYdAnPID","YdE vs Angle with an PID gate", 16,Ydbins,10000,0,5); //gate on the Sds
    TH2D *hYdCsIch = new TH2D("hYdCsIch","Ring Yd vs Sector CsI",16,0,16,16,0,16);
    TH1D *hYuEn = new TH1D("hYuEn","hYuEn",2500,0,10);
    TH2D *hYdCsI[16];
    TH2D *hYdCsI2[16];
    
    for(int i=0; i<16; i++)
    {
     string nameYdCsI = Form("YdCsI_%i",i);
     hYdCsI[i] = new TH2D(nameYdCsI.c_str(),"Yd vs CsI1 sector",2048,0,4096,10000,0,5);
     
     string nameYdCsI2 = Form("YdCsI2_%i",i);
     hYdCsI2[i] = new TH2D(nameYdCsI2.c_str(),"Yd vs CsI2 sector",2048,0,4096,10000,0,5);
    }
        
    //Yu detector => interesting reaction
    //calculate the bins for the Yu detector
    double Yubins[17]={0};
    for (int i=1; i<18; i++){Yubins[i-1]=180-TMath::RadToDeg()*TMath::ATan((50+(16-i)*4.94)/85);}
    
    //for (int j=0; j<17; j++){cout << j << " /" << Yubins[j] << endl;}
    TH2D *hYuAn = new TH2D("hYuAn","YuE vs Angle",16,Yubins,2500,0,10); //no gates at the moment
    TH2D *hYuAnIC = new TH2D("hYuAnIC","YuE vs Angle with a gate on the IC",16,Yubins,60,0,3); //gate on the IC
    TH2D *hYuAnIC2 = new TH2D("hYuAnIC2","YuE vs Angle with a gate on the IC2",16,Yubins,60,0,3); //gate on the second IC peak in the Be data
    TH2D *hYuAnPID1 = new TH2D("hYuAnPID1","YuE vs Angle with a PID1 gate",16,Yubins,60,0,3); //gate on the Sds
    TH2D *hYuAnPID2 = new TH2D("hYuAnPID2","YuE vs Angle with a PID2 gate",16,Yubins,10000,0,10); //gate on the Sds
    TH2D *hYuAnPID3 = new TH2D("hYuAnPID3","YuE vs Angle with a PID3 gate",16,Yubins,10000,0,10); //gate on the Sds
    TH2D *hYuAnPID4 = new TH2D("hYuAnPID4","YuE vs Angle with a PID4 gate",16,Yubins,10000,0,10); //gate on the Sds
    TH2D *hYuAnPIDW = new TH2D("hYuAnPIDW","YuE vs Angle with a PID wide gate",16,Yubins,60,0,3); //gate on the Sds
    
    TH2D *hYuAnPIDWS[8];
    for(int i=0; i<8; i++){
      string nameS = Form("YuEAS_%i",i);
      hYuAnPIDWS[i] = new TH2D(nameS.c_str(),"YuE vs Angle with a PID wide gate for sector",16,Yubins,60,0,3);}
    
    TH1D *hYuES[8];
    TH1D *hYuER[16];
    TH2D *hYuEAS[8];
    TH2D *hYuEAR[16];
    
    TH1D *hYuESIC[8];
    TH1D *hYuERIC[16];
    TH2D *hYuEASIC[8];
    TH2D *hYuEARIC[16];
    
    for(int i=0; i<8; i++)
    {
        string name1 = Form ( "hYuES_%i",i );
        hYuES[i] = new TH1D ( name1.c_str(),"Yu energy for sector",2500,0,10 );
      
        string name2 = Form ( "hYuEAS_%i",i );
	hYuEAS[i]= new TH2D(name2.c_str(),"Yu en vs angle for sector",16,Yubins,2500,0,10);
	
	string name1IC = Form ( "hYuESIC_%i",i );
        hYuESIC[i] = new TH1D ( name1IC.c_str(),"Yu energy for sector gated by IC",60,0,3 );
      
        string name2IC = Form ( "hYuEASIC_%i",i );
	hYuEASIC[i]= new TH2D(name2IC.c_str(),"Yu en vs angle for sector gated by IC",16,Yubins,60,0,3);
    }
    
    for(int i=0; i<16; i++)
    {
        string name3 = Form ( "hYuER_%i",i );
        hYuER[i] = new TH1D ( name3.c_str(),"Yu energy for ring",1250,0,10 );
      
        string name4 = Form ( "hYuEAR_%i",i );
	hYuEAR[i]= new TH2D(name4.c_str(),"Yu en vs angle for ring",16,Yubins,1250,0,10);
	
	string name3IC = Form ( "hYuERIC_%i",i );
        hYuERIC[i] = new TH1D ( name3IC.c_str(),"Yu energy for ring gated by IC",60,0,3 );
      
        string name4IC = Form ( "hYuEARIC_%i",i );
	hYuEARIC[i]= new TH2D(name4IC.c_str(),"Yu en vs angle for ring gated by IC",16,Yubins,60,0,3);
    }
    
    //Q values
    TH1D *hQval = new TH1D("hQval","Q values",100,-10,10);
    
     //start reading the tree
    int ev = chain->GetEntries();
    cout << "Total number of events =" << ev << endl;
    int ev_num=0;
    
    //define the variables for Qvalues calculations and angle calculations
    //stuff for the YY detectors
    float YChWidth = ( Ydr1 - Ydr0 ) /16.;
    float YChMiddle = YChWidth/2;
    float yuM=0, ydM=0;
    double YuthetaM, YdthetaM; //angle for Yu/Yd
    TVector3 beam(0,0,1); //beam vector
    
    //stuff for the Q value calculations
    float amu = 931.5; // MeV
    float massEjec = 938.28; //MeV/c2 proton
    //float kBeam = 114; //kinetic energy of the beam in MeV
    float kBeam = 108; //9 AMeV
    float mbeam = 12 * amu;  // MeV //mass of the beam 12Be or 12C
    float mrecoil = 13 * amu;  // MeV // mass of the recoils => 13Be or 13C
    float mejec = 1 * amu; //mass of the ejectiles => protons
    
    double Qval;
    
    vector<double> *YuAngle = new vector<double>;
    
    TF1* fit_func6 = new TF1 ( "fit_func6","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp(-pow((x-[10]),2)/(2*[2]*[2]))+[11]*TMath::Exp(-pow((x-[12]),2)/(2*[2]*[2]))+[13]*x+[14]",0,10 );
    
    while(ev_num<=ev)
    {
      if ( ev_num%10000==0 ) cout << "Current event = " << ev_num << "\r"<< flush;
      chain->GetEntry ( ev_num );
      
      //fill in the PID histograms
      //non calibrated
      if(Sd1rChannel->size()>0 && Sd2rChannel->size()>0 && Sd1rMul==1 && Sd2rMul==1)
      {//cout << " In channel 1 "<< endl;
	hSd1rSd2r->Fill(Sd2rEnergyRaw->at(0),Sd1rEnergyRaw->at(0));
	if(TMath::Abs(Sd2rChannel->at(0)-Sd1rChannel->at(0))==1){hSd1rSd2rCut->Fill(Sd2rEnergyRaw->at(0),Sd1rEnergyRaw->at(0));}
         hSdr1r2Ch->Fill(Sd2rChannel->at(0),Sd1rChannel->at(0));      } //the non calibrated Sd1r vs Sd2r
      
      if(Sd1rChannel->size()>0 && Sd2sChannel->size()>0 && Sd1rMul==1 && Sd2sMul==1)
      {//cout << " In channel 2 "<< endl;
	hSd1rSd2s->Fill(Sd2sEnergyRaw->at(0),Sd1rEnergyRaw->at(0));
      //  hSdr1s2Ch->Fill(Sd2sChannel->at(0),Sd1rChannel->at(0));  

      } //the non calibrated Sd1r vs Sd2s
      
      if(Sd1sChannel->size()>0 && Sd2rChannel->size()>0 && Sd1sMul==1 && Sd2rMul==1)
      {//cout << " In channel 3 "<< endl;
	hSd1sSd2r->Fill(Sd2rEnergyRaw->at(0),Sd1sEnergyRaw->at(0));
      // hSds1r2Ch->Fill ( Sd2rChannel->at ( 0 ),Sd1sChannel->at ( 0 ) );
      } //the non calibrated Sd1s vs Sd2r
      
      if(Sd1sChannel->size()>0 && Sd2sChannel->size()>0 && Sd1sMul==1 && Sd2sMul==1)
      {//cout << " In channel 4 "<< endl;
	hSd1sSd2s->Fill(Sd2sEnergyRaw->at(0),Sd1sEnergyRaw->at(0));
       hSds1s2Ch->Fill ( Sd2sChannel->at ( 0 ),Sd1sChannel->at ( 0 ) );
      } //the non calibrated Sd1s vs Sd2s
      
      //calibratedBe00_5050.root
      if(Sd1rEnergy->size()>0 && Sd2rEnergy->size()>0 && Sd1rMul==1 && Sd2rMul==1)
      {//cout << " In sd1r sd2r " << endl;
	//if(TMath::Abs(Sd2rChannel->at(0)-Sd1rChannel->at(0))==1){hCSd1rSd2r->Fill(Sd2rEnergy->at(0),Sd1rEnergy->at(0));}
	 if ( ICEnergyRaw->size() >0 ){if ( ICEnergyRaw->at ( 0 ) >600 && ICEnergyRaw->at ( 0 ) <1100){hCSd1rSd2rIC->Fill ( Sd2rEnergy->at ( 0 ),Sd1rEnergy->at ( 0 ) );}}  
            if ( ICEnergyRaw->size() >0 ){if ( ICEnergyRaw->at ( 0 ) >600 && ICEnergyRaw->at ( 0 ) <1100 && TMath::Abs(Sd2rChannel->at(0)-Sd1rChannel->at(0))==1){hCSd1rSd2rICCut->Fill ( Sd2rEnergy->at ( 0 ),Sd1rEnergy->at ( 0 ) );}}   //IC gate for C 1500 < E < 200; for Be 600< E < 1100
            if ( ICEnergyRaw->size() >0 ){if ( ICEnergyRaw->at ( 0 ) <500 ){hCSd1rSd2rIC2->Fill ( Sd2rEnergy->at ( 0 ),Sd1rEnergy->at ( 0 ) );}}   //IC gate for second IC peak in Be data
            if ( YdEnergy->size() >0 && YdMul==1){ if(YdEnergy->at( 0 ) >0){hCSd1rSd2rYd->Fill(Sd2rEnergy->at(0),Sd1rEnergy->at(0));}}
            if ( YdEnergy->size() >0 && YdMul==1 && ICEnergyRaw->size() >0 ){ if(YdEnergy->at( 0 ) >0 && ICEnergyRaw->at ( 0 ) >600 && ICEnergyRaw->at ( 0 ) <1100)
	      {hCSd1rSd2rYdIC->Fill(Sd2rEnergy->at(0),Sd1rEnergy->at(0));}}
	      
	      if ( YdEnergy->size() >0 && YdMul==1 && ICEnergyRaw->size() >0 && TMath::Abs(Sd2rChannel->at(0)-Sd1rChannel->at(0))==1){ if(YdEnergy->at( 0 ) >0 && ICEnergyRaw->at ( 0 ) >600 && ICEnergyRaw->at ( 0 ) <1100)
	      {hCSd1rSd2rYdICCut->Fill(Sd2rEnergy->at(0),Sd1rEnergy->at(0));}}
	      
	      if(YdEnergy->size() >0 && CsI2Channel->size()>0 && CsI2Channel->at(0)==0 && YdMul==1 && ICEnergyRaw->size() >0 && TMath::Abs(Sd2rChannel->at(0)-Sd1rChannel->at(0))==1 && cutYd1->IsInside(CsI2EnergyRaw->at(0),YdEnergy->at(0)))
	      {hSd1r2rYdBand1->Fill(Sd2rEnergy->at(0),Sd1rEnergy->at(0));}
	      
	      if(YdEnergy->size() >0 && CsI2Channel->size()>0 && CsI2Channel->at(0)==0 && YdMul==1 && ICEnergyRaw->size() >0 && TMath::Abs(Sd2rChannel->at(0)-Sd1rChannel->at(0))==1 && cutYd2->IsInside(CsI2EnergyRaw->at(0),YdEnergy->at(0)))
	      {hSd1r2rYdBand2->Fill(Sd2rEnergy->at(0),Sd1rEnergy->at(0));}
	      
	      if(YdEnergy->size() >0 && CsI2Channel->size()>0 && CsI2Channel->at(0)==0 && YdMul==1 && ICEnergyRaw->size() >0 && TMath::Abs(Sd2rChannel->at(0)-Sd1rChannel->at(0))==1 && cutYd3->IsInside(CsI2EnergyRaw->at(0),YdEnergy->at(0)))
	      {hSd1r2rYdBand3->Fill(Sd2rEnergy->at(0),Sd1rEnergy->at(0));}
	
      } 
      
      if(Sd1rEnergy->size()>0 && Sd2sEnergy->size()>0 && Sd1rMul==1 && Sd2sMul==1)
      {//cout << " In sd1r sd2s " << endl;
	hCSd1rSd2s->Fill(Sd2sEnergy->at(0),Sd1rEnergy->at(0));
	if ( ICEnergyRaw->size() >0){ if(ICEnergyRaw->at( 0 ) >600 && ICEnergyRaw->at ( 0 ) <1100 ){hCSd1rSd2sIC->Fill(Sd2sEnergy->at(0),Sd1rEnergy->at(0));}}
	if ( YdEnergy->size() >0 && YdMul==1){ if(YdEnergy->at( 0 ) >0){hCSd1rSd2sYd->Fill(Sd2sEnergy->at(0),Sd1rEnergy->at(0));}}
	if ( YdEnergy->size() >0 && YdMul==1 && ICEnergyRaw->size() >0){ if(YdEnergy->at( 0 ) >0 && ICEnergyRaw->at ( 0 ) >600 && ICEnergyRaw->at ( 0 ) <1100)
	  {hCSd1rSd2sYdIC->Fill(Sd2sEnergy->at(0),Sd1rEnergy->at(0));}}
      } 
      
      if(Sd1sEnergy->size()>0 && Sd2rEnergy->size()>0 && Sd1sMul==1 && Sd2rMul==1)
      {//cout << " In sd1s sd2r " << endl;
	hCSd1sSd2r->Fill(Sd2rEnergy->at(0),Sd1sEnergy->at(0));
	if ( ICEnergyRaw->size() >0){ if(ICEnergyRaw->at( 0 ) >600 && ICEnergyRaw->at ( 0 ) <1100 ){hCSd1sSd2rIC->Fill(Sd2rEnergy->at(0),Sd1sEnergy->at(0));}}
	if ( YdEnergy->size() >0 && YdMul==1){ if(YdEnergy->at( 0 ) >0){hCSd1sSd2rYd->Fill(Sd2rEnergy->at(0),Sd1sEnergy->at(0));}}
	if ( YdEnergy->size() >0 && YdMul==1 && ICEnergyRaw->size() >0){ if(YdEnergy->at( 0 ) >0 && ICEnergyRaw->at ( 0 ) >600 && ICEnergyRaw->at ( 0 ) <1100)
	  {hCSd1sSd2rYdIC->Fill(Sd2rEnergy->at(0),Sd1sEnergy->at(0));}}
	
      } 
      
      if(Sd1sEnergy->size()>0 && Sd2sEnergy->size()>0 && Sd1sMul==1 && Sd2sMul==1)
      {//cout << " In sd1s sd2s " << endl;
	hCSd1sSd2s->Fill(Sd2sEnergy->at(0),Sd1sEnergy->at(0));
	if ( ICEnergyRaw->size() >0){ if(ICEnergyRaw->at( 0 ) >600 && ICEnergyRaw->at ( 0 ) <1100 ){hCSd1sSd2sIC->Fill(Sd2sEnergy->at(0),Sd1sEnergy->at(0));}}
	if ( YdEnergy->size() >0 && YdMul==1){ if(YdEnergy->at( 0 ) >0){hCSd1sSd2sYd->Fill(Sd2sEnergy->at(0),Sd1sEnergy->at(0));}}
	if ( YdEnergy->size() >0 && YdMul==1 && ICEnergyRaw->size() >0){ if(YdEnergy->at( 0 ) >0 && ICEnergyRaw->at ( 0 ) >600 && ICEnergyRaw->at ( 0 ) <1100)
	  {hCSd1sSd2sYdIC->Fill(Sd2sEnergy->at(0),Sd1sEnergy->at(0));}}
      } 
      
      //fill in the IC
      if(ICEnergyRaw->size()>0){hIC->Fill(ICEnergyRaw->at(0));}
      
      //fill in the Yd detector
      if(YdChannel->size()>0 && YdMul==1 && YdEnergy->size()>0)
      {
            if ( YdEnergy->at ( 0 ) >0 && ( YdChannel->at ( 0 ) !=20 && YdChannel->at ( 0 ) !=0 && YdChannel->at ( 0 ) !=55 &&  YdChannel->at ( 0 ) !=127 ) )
            {
	      //cout << " In Yd stuff "<< endl;
                ydM = Ydr0+ ( YdRing->at ( 0 ) *YChWidth )-YChMiddle;
                TVector3 YdposM ( 0,ydM,Ydz );
                YdposM.SetPhi ( TMath::Pi() /2-YdSector->at ( 0 ) *TMath::Pi() /4 ); //Pi/2 because the center of the sector is at 90degrees, Pi/4 is because there is 8 sectors so it is 2Pi/8
                YdthetaM = beam.Angle ( YdposM ) *TMath::RadToDeg();
		
                hYdAn->Fill ( YdthetaM,YdEnergy->at ( 0 ) ); //all, no cuts
                if(ICEnergyRaw->size()>0) { if(ICEnergyRaw->at(0)>600 && ICEnergyRaw->at(0)<1100){hYdAnIC->Fill(YdthetaM,YdEnergy->at ( 0 ));}} //IC gate for C 1500 < E < 200; for Be 600< E < 1100
                if(ICEnergyRaw->size()>0) { if(ICEnergyRaw->at(0)<500){hYdAnIC2->Fill(YdthetaM,YdEnergy->at ( 0 ));}} //IC gate for 2nd Be peak
                if(Sd1rEnergy->size()>0 && Sd2rEnergy->size()>0 && Sd1rMul==1 && Sd2rMul==1){if(cutIc1->IsInside(Sd2rEnergy->at(0),Sd1rEnergy->at(0))){hYdAnPID->Fill ( YdthetaM,YdEnergy->at ( 0 ) );}}
                
                if(CsI1Channel->size()>0)if(CsI1EnergyRaw->at(0)>0 && CsI1Mul==1){{
		  hYdCsIch->Fill(YdChannel->at(0),CsI1Channel->at(0));
		  hYdCsI[CsI1Channel->at(0)]->Fill(CsI1EnergyRaw->at(0),YdEnergy->at(0)); }}
		  
		  if(CsI2Channel->size()>0)if(CsI2EnergyRaw->at(0)>0 && CsI2Mul==1){{
		    hYdCsI2[CsI2Channel->at(0)]->Fill(CsI2EnergyRaw->at(0),YdEnergy->at(0)); }}
            }
      }
      
      //fill in the Yu detector
      if(YuChannel->size()>0  && YuEnergy->size()>0 && YuEnergy->size()==YuChannel->size()) //removed the multiplicity because of alpha test
      {
	//cout << " I am in the loop "<< endl;
	for(unsigned int i=0; i<YuChannel->size(); i++)
	{
	 if (YuEnergy->at ( i ) >0.20 && ( YuChannel->at ( i ) !=80 && YuChannel->at ( i ) !=95 && YuChannel->at ( i ) !=94 && YuChannel->at ( i ) !=82 && YuChannel->at ( i ) !=96 && YuChannel->at ( i ) !=106 && YuChannel->at ( i ) !=111 ))
	 {
	  // if(YuEnergy->at ( i )>1) {cout << "En " << YuEnergy->at ( i ) << " Channel " << YuChannel->at ( i ) << endl;}
                    yuM = Ydr0+ ( YuRing->at ( i ) *YChWidth )-YChMiddle;
                    TVector3 YuposM ( 0,yuM,Yuz );
                    YuposM.SetPhi ( TMath::Pi() /2-YuSector->at ( i ) *TMath::Pi() /4 ); //Pi/2 because the center of the sector is at 90degrees, Pi/4 is because there is 8 sectors so it is 2Pi/8
                    YuthetaM = beam.Angle ( YuposM ) *TMath::RadToDeg();
		    YuAngle->push_back(YuthetaM);
		    
                    hYuES[YuSector->at ( i )]->Fill ( YuEnergy->at ( i ) );
                    hYuER[YuRing->at ( i )]->Fill ( YuEnergy->at ( i ) );
                    hYuEAS[YuSector->at ( i )]->Fill ( YuthetaM,YuEnergy->at ( i ) );
                    hYuEAR[YuRing->at ( i )]->Fill ( YuthetaM,YuEnergy->at ( i ) );
		    
		   hYuAn->Fill ( YuthetaM,YuEnergy->at ( i ) ); //all, no cuts
		   hYuEn->Fill(YuEnergy->at ( i ));
		   
                    if ( ICEnergyRaw->size() >0 )
                    {
		      if ( ICEnergyRaw->at ( 0 ) >600 && ICEnergyRaw->at ( 0 ) <1100 )
                        {
                            hYuAnIC->Fill ( YuthetaM,YuEnergy->at ( i ) );

                            hYuESIC[YuSector->at ( i )]->Fill ( YuEnergy->at ( i ) );
                            hYuERIC[YuRing->at ( i )]->Fill ( YuEnergy->at ( i ) );
                            hYuEASIC[YuSector->at ( i )]->Fill ( YuthetaM,YuEnergy->at ( i ) );
                            hYuEARIC[YuRing->at ( i )]->Fill ( YuthetaM,YuEnergy->at ( i ) );

                        }
                    }
                    
                   if(Sd2rEnergy->size()>0 && Sd1rEnergy->size()>0 && ICEnergyRaw->size()>0)
		    {
                      if ( cutBeNew->IsInside ( Sd2rEnergy->at ( 0 ),Sd1rEnergy->at ( 0 ) ) && ICEnergyRaw->at ( 0 ) >600 && ICEnergyRaw->at ( 0 ) <1100 )
                    {
		      
                            hYuAnPIDWS[YuSector->at ( i )]->Fill ( YuthetaM,YuEnergy->at ( i ) );
			    
                        hYuAnPIDW->Fill ( YuthetaM,YuEnergy->at ( i ) );
                        Qval = ( 1+mejec/mrecoil ) * ( YuEnergy->at ( i ) ) - ( 1 - mbeam/mrecoil ) * ( kBeam ) - 2 * TMath::Sqrt ( mbeam*mejec* ( YuEnergy->at ( i ) ) * ( kBeam ) ) / mrecoil * TMath::Cos ( YuthetaM * TMath::Pi() / 180. );
                        hQval->Fill ( Qval );
    
                    }}
		   
	 }
	  
	}
	//old stuff
           /* if ( YuEnergy->at ( 0 ) >0 && ( YuChannel->at ( 0 ) !=82 && YuChannel->at ( 0 ) !=96 && YuChannel->at ( 0 ) !=106  ) ) //<112 for C and <111 for Be data
            {
	      //cout << "In Yu stuff " << endl;
                yuM = Ydr0+ ( YuRing->at ( 0 ) *YChWidth )-YChMiddle;
                TVector3 YuposM ( 0,yuM,Yuz );
                YuposM.SetPhi ( TMath::Pi() /2-YuSector->at ( 0 ) *TMath::Pi() /4 ); //Pi/2 because the center of the sector is at 90degrees, Pi/4 is because there is 8 sectors so it is 2Pi/8
                YuthetaM = beam.Angle ( YuposM ) *TMath::RadToDeg();
               // cout << "Ring " << YuRing->at ( 0 ) << " angle " <<  YuthetaM << endl;

                hYuAn->Fill ( YuthetaM,YuEnergy->at ( 0 ) ); //all, no cuts
		hYuEn->Fill(YuEnergy->at ( 0 ));
		
                hYuES[YuSector->at ( 0 )]->Fill ( YuEnergy->at ( 0 ) );
                hYuER[YuRing->at ( 0 )]->Fill ( YuEnergy->at ( 0 ) );
                hYuEAS[YuSector->at ( 0 )]->Fill ( YuthetaM,YuEnergy->at ( 0 ) );
                hYuEAR[YuRing->at ( 0 )]->Fill ( YuthetaM,YuEnergy->at ( 0 ) );
		
                if ( ICEnergyRaw->size() >0){ if(ICEnergyRaw->at( 0 ) >600 && ICEnergyRaw->at ( 0 ) <1100 ){hYuAnIC->Fill ( YuthetaM,YuEnergy->at ( 0 ) );
		  
		hYuESIC[YuSector->at ( 0 )]->Fill ( YuEnergy->at ( 0 ) );
                hYuERIC[YuRing->at ( 0 )]->Fill ( YuEnergy->at ( 0 ) );
                hYuEASIC[YuSector->at ( 0 )]->Fill ( YuthetaM,YuEnergy->at ( 0 ) );
                hYuEARIC[YuRing->at ( 0 )]->Fill ( YuthetaM,YuEnergy->at ( 0 ) );
		  
		}}
  
                if ( ICEnergyRaw->size() >0){ if(ICEnergyRaw->at ( 0 ) <500 ){hYuAnIC2->Fill ( YuthetaM,YuEnergy->at ( 0 ) );}}
                
                if(ICEnergyRaw->size() >0 && Sd1rEnergy->size()>0 && Sd2rEnergy->size()>0 && Sd1rMul==1 && Sd2rMul==1 ){if(cutIc5->IsInside(Sd2rEnergy->at(0),Sd1rEnergy->at(0))&& ICEnergyRaw->at( 0 ) >600 && ICEnergyRaw->at ( 0 ) <1100){hYuAnPID1->Fill ( YuthetaM,YuEnergy->at ( 0 ) );
               
                }
                    if ( cutBe6->IsInside ( Sd2rEnergy->at ( 0 ),Sd1rEnergy->at ( 0 ) ) && ICEnergyRaw->at ( 0 ) >600 && ICEnergyRaw->at ( 0 ) <1100 )
                    {
                        hYuAnPIDW->Fill ( YuthetaM,YuEnergy->at ( 0 ) );
                        Qval = ( 1+mejec/mrecoil ) * ( YuEnergy->at ( 0 ) ) - ( 1 - mbeam/mrecoil ) * ( kBeam ) - 2 * TMath::Sqrt ( mbeam*mejec* ( YuEnergy->at ( 0 ) ) * ( kBeam ) ) / mrecoil * TMath::Cos ( YuthetaM * TMath::Pi() / 180. );
                        hQval->Fill ( Qval );
    
                    }

                }
                
               // if(Sd1rEnergy->size()>0 && Sd2sEnergy->size()>0 && Sd1rMul==1 && Sd2sMul==1){if(cutIc2->IsInside(Sd2sEnergy->at(0),Sd1rEnergy->at(0))){hYuAnPID2->Fill ( YuthetaM,YuEnergy->at ( 0 ) );}}
               // if(Sd1sEnergy->size()>0 && Sd2rEnergy->size()>0 && Sd1sMul==1 && Sd2rMul==1){if(cutIc3->IsInside(Sd2rEnergy->at(0),Sd1sEnergy->at(0))){hYuAnPID3->Fill ( YuthetaM,YuEnergy->at ( 0 ) );}}
               // if(Sd1sEnergy->size()>0 && Sd2sEnergy->size()>0 && Sd1sMul==1 && Sd2sMul==1){if(cutIc4->IsInside(Sd2sEnergy->at(0),Sd1sEnergy->at(0))){hYuAnPID4->Fill ( YuthetaM,YuEnergy->at ( 0 ) );}}
            }*/
      }
      
      //look into the PID detector with an Yu cut on Be
        if(Sd1rEnergy->size()>0 && Sd2rEnergy->size()>0 && Sd1rMul==1 && Sd2rMul==1 && YuEnergy->size()>0)
      { for(int i=0; i< Sd1rEnergy->size(); i++)
	{
	if ( ICEnergyRaw->size() >0 ){if ( ICEnergyRaw->at ( 0 ) >600 && ICEnergyRaw->at ( 0 ) <1100 && TMath::Abs(Sd2rChannel->at(0)-Sd1rChannel->at(0))==1 && cutBeKine->IsInside(YuAngle->at(i),YuEnergy->at(i))){
      hCSd1rSd2rICYu->Fill(Sd2rEnergy->at(0),Sd1rEnergy->at(0));
      }}}}
      
      ev_num++;
    } //end of the main while loop
    
    YuAngle->clear();
    
    
    //Search for peaks in Yd strips
    int npeaks=3;
    TSpectrum* s = new TSpectrum ( 2*npeaks );
    int nfound=0;
    double* txp = s->GetPositionX();
    
    double xp[3];
    double aS[8],bS[8],cS[8],aR[16], bR[16], cR[16];

    //serch in the sectors
   /* for ( int j=0; j<8; j++ )
    {
        nfound = s->Search ( hYuES[j],2,"",0.1 );
      
            for ( int p=0; p<nfound; p++ ){xp[p] = txp[p];}
          
     aS[j] = min ( min (xp[0],xp[1]), xp[2]);
     cS[j] = max ( max (xp[0],xp[1]), xp[2]);
    
    if(xp[0]!=aS[j] && xp[0]!=cS[j]) {bS[j] = xp[0];}
    else if(xp[1]!=aS[j] && xp[1]!=cS[j]) {bS[j] = xp[1];}
    else if(xp[2]!=aS[j] && xp[2]!=cS[j]) {bS[j] = xp[2];}
       
          // cout << j << "pedestal " << ped[j] << " peak 1 " << a[j] << " peak 2 " << b[j] << " peak 3 " << c[j] << endl;
    }
    
      for ( int i=0; i<8; i++ )
    {          
        fit_func6->SetParLimits ( 1,aS[i]-0.08,aS[i]-0.025 );
        fit_func6->SetParLimits ( 4,aS[i]-0.025,aS[i]+0.03 );
        fit_func6->SetParLimits ( 6,bS[i]-0.08,bS[i]-0.025 );
        fit_func6->SetParLimits ( 8,bS[i]-0.025,bS[i]+0.03 );
        fit_func6->SetParLimits ( 10,cS[i]-0.08,cS[i]-0.025 );
        fit_func6->SetParLimits ( 12,cS[i]-0.025,cS[i]+0.03 );
        fit_func6->SetParLimits ( 2,0,0.03 );
	
	 hYuES[i]->Fit ( "fit_func6","","",aS[i]-0.5,cS[i]+0.5 );
    }*/
    
    //search in the rings
 /*    for ( int j=0; j<16; j++ )
    {
        nfound = s->Search ( hYuER[j],2,"",0.1 );
      
        for ( int p=0; p<nfound; p++ ){xp[p] = txp[p];}
          
     aR[j] = min ( min (xp[0],xp[1]), xp[2]);
     cR[j] = max ( max (xp[0],xp[1]), xp[2]);
    
    if(xp[0]!=aR[j] && xp[0]!=cR[j]) {bR[j] = xp[0];}
    else if(xp[1]!=aR[j] && xp[1]!=cR[j]) {bR[j] = xp[1];}
    else if(xp[2]!=aR[j] && xp[2]!=cR[j]) {bR[j] = xp[2];}
       
          // cout << j << "pedestal " << ped[j] << " peak 1 " << a[j] << " peak 2 " << b[j] << " peak 3 " << c[j] << endl;
    }
    
      for ( int i=0; i<16; i++ )
    {          
        fit_func6->SetParLimits ( 1,aR[i]-0.08,aR[i]-0.025 );
        fit_func6->SetParLimits ( 4,aR[i]-0.025,aR[i]+0.03 );
        fit_func6->SetParLimits ( 6,bR[i]-0.08,bR[i]-0.025 );
        fit_func6->SetParLimits ( 8,bR[i]-0.025,bR[i]+0.03 );
        fit_func6->SetParLimits ( 10,cR[i]-0.08,cR[i]-0.025 );
        fit_func6->SetParLimits ( 12,cR[i]-0.025,cR[i]+0.03 );
        fit_func6->SetParLimits ( 2,0,0.03 );
	
	 hYuER[i]->Fit ( "fit_func6","","",aR[i]-0.5,cR[i]+0.5 );
    }*/
    
    
    f_out->cd();
    
    //write the histograms
    //PID 
   // hSdr1r2Ch->Write();
  //  hSdr1s2Ch->Write();
    //hSds1r2Ch->Write();
  //  hSds1s2Ch->Write();
    
 //   hSd1rSd2r->Write(); 
 //   hSd1rSd2rCut->Write(); 
  /*  hSd1rSd2s->Write();
    hSd1sSd2r->Write();
    hSd1sSd2s->Write();*/

    
    hCSd1rSd2s->Write();
    hCSd1rSd2sIC->Write();
   /* hCSd1sSd2r->Write();
    hCSd1sSd2s->Write();*/
   
    hCSd1rSd2r->Write();
    hCSd1rSd2rIC->Write();
    hCSd1rSd2rICCut->Write();
  /*  hSd1r2rYdBand1->Write();
    hSd1r2rYdBand2->Write();
    hSd1r2rYdBand3->Write();
    
    hCSd1rSd2rIC2->Write();*/
 //   hCSd1rSd2sIC->Write();
 //   hCSd1sSd2rIC->Write();
 //   hCSd1sSd2sIC->Write();
    
  /*  hCSd1rSd2rYd->Write();
    hCSd1rSd2sYd->Write();
    hCSd1sSd2rYd->Write();
    hCSd1sSd2sYd->Write();
    
    hCSd1rSd2rYdIC->Write();
    hCSd1rSd2rYdICCut->Write();
    hCSd1rSd2sYdIC->Write();
    hCSd1sSd2rYdIC->Write();
    hCSd1sSd2sYdIC->Write();*/
    
    
    //IC
    hIC->Write(); 
    
    //Yd
    hYdAn->Write(); 
    hYdAnIC->Write(); 
  /*  hYdAnIC2->Write(); 
    hYdAnPID->Write();
    hYdCsIch->Write();
    for(int i=0; i<16; i++){hYdCsI[i]->Write();}
    for(int i=0; i<16; i++){hYdCsI2[i]->Write();}*/
    
    //Yu
    hYuAn->Write(); 
    hYuAnIC->Write(); 
  //  hYuAnIC2->Write(); 
  //  hYuAnPID1->Write(); 
    hYuAnPIDW->Write(); 
 //   hYuAnPID2->Write(); 
 //   hYuAnPID3->Write(); 
 //   hYuAnPID4->Write(); 
     hYuEn->Write();
    hQval->Write();
    for(int i=0; i<8; i++){hYuAnPIDWS[i]->Write();}
   
    hCSd1rSd2rICYu->Write();
    
  /*  for(int i=0; i<8; i++){hYuES[i]->Write();}
   for(int i=0; i<7; i++){hYuESIC[i]->Write();}
    for(int i=0; i<8; i++){ hYuEAS[i]->Write();}
    for(int i=0; i<7; i++){ hYuEASIC[i]->Write();}
    for ( int j=0; j<16; j++ ){hYuER[j]->Write();}
    for ( int j=0; j<16; j++ ){hYuERIC[j]->Write();}
    for ( int j=0; j<16; j++ ){hYuEAR[j]->Write();}
    for ( int j=0; j<16; j++ ){hYuEARIC[j]->Write();}*/
    
    f_out->Close();
    
    return 0;
    
} //end of the program