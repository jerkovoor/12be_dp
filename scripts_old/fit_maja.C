#include "TMath.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TChain.h"
#include "TString.h"
#include <iostream>
#include <vector>
#include <string>

vector<TH1F*> hcoinc;

//Count the number of events in a gaussian fit
Double_t sqrt2pi = TMath::Sqrt(2*TMath::Pi());

//Log half life determination

//log fit with background taken into account
TF1 *log_halfLife = new TF1("log_halfLife","[0]/[1]*TMath::Exp(x*TMath::Log(2))*TMath::Exp(-1/[1]*TMath::Exp(x*TMath::Log(2)))*TMath::Exp(-1/[3]*TMath::Exp(x*TMath::Log(2)))+[2]/[3]*TMath::Exp(x*TMath::Log(2))*TMath::Exp(-1/[3]*TMath::Exp(x*TMath::Log(2)))",3,20);

//[0] amplitude of the signal
//[1] searched half life
//[2] amplitude of the background
//[3] background half life

//log fit without a background
TF1 *log_hl_1comp = new TF1("log_hl_1comp","[0]/[1]*TMath::Exp(x*TMath::Log(2))*TMath::Exp(-1/[1]*TMath::Exp(x*TMath::Log(2)))",3,20);

//[0] amplitude of the signal
//[1] searched half life

//gaussian fits

// 1 gaussian
TF1 *gaus_pol1 = new TF1("gaus_pol1","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]+[4]*x",0,6000);


//2 gaussians
TF1 *gaus2_pol1 = new TF1("gaus2_pol1","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]+[6]*x",0,6000);

//3 gaussians
TF1 *gaus3_pol1 = new TF1("gaus3_pol1","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]+[8]*x",0,6000);

//4 gaussians
TF1 *gaus4_pol1 = new TF1("gaus4_pol1","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]+[10]*x",0,6000);

//5 gaussians
TF1 *gaus5_pol1 = new TF1("gaus5_pol1","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp(-pow((x-[10]),2)/(2*[2]*[2]))+[11]+[12]*x",0,6000);

//6 gaussians
TF1 *gaus6_pol1 = new TF1("gaus6_pol1","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp(-pow((x-[10]),2)/(2*[2]*[2]))+[11]*TMath::Exp(-pow((x-[12]),2)/(2*[2]*[2]))+[13]+[14]*x",0,6000);

//7 gaussians
TF1 *gaus7_pol1 = new TF1("gaus7_pol1","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp(-pow((x-[10]),2)/(2*[2]*[2]))+[11]*TMath::Exp(-pow((x-[12]),2)/(2*[2]*[2]))+[13]*TMath::Exp(-pow((x-[14]),2)/(2*[2]*[2]))+[15]+[16]*x",0,6000);

//9 gaussians
TF1 *gaus9_pol1 = new TF1("gaus9_pol1","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp(-pow((x-[10]),2)/(2*[2]*[2]))+[11]*TMath::Exp(-pow((x-[12]),2)/(2*[2]*[2]))+[13]*TMath::Exp(-pow((x-[14]),2)/(2*[2]*[2]))+[15]*TMath::Exp(-pow((x-[16]),2)/(2*[2]*[2]))+[17]*TMath::Exp(-pow((x-[18]),2)/(2*[2]*[2]))+[19]+[20]*x",0,6000);

//exponential fit

//simple exponential
TF1 *expo_simple = new TF1("expo_simple","[0] * TMath::Log(2)/[1] * TMath::Exp(-TMath::Log(2)/[1] * x)",0,3000);
//[0] amplitude
//[1] half life (or equivalent)

//simple exponential with a linear fit
TF1 *expo_lin = new TF1("expo_lin","[0] * TMath::Log(2)/[1] * TMath::Exp(-TMath::Log(2)/[1] * x) + [2] * x + [3]",0,100000);
//[0] amplitude
//[1] half life (or equivalent)
//[2] 
//[3] background offset on the y-axis

//2 combined exponentials
TF1 *expo_wbdf = new TF1("expo_wbdf","[0] * TMath::Log(2)/[1] * TMath::Exp(-TMath::Log(2)/[1] * x) + [2] * TMath::Log(2)/[3] * TMath::Exp(-TMath::Log(2)/[3] * x)",0,3000);
//[0] amplitude first 
//[1] half life (or equivalent) first
//[2] amplitude background
//[3] half life background

//2 half lives - beta with background
TF1 *hl_feed = new TF1("hl_feed","([0]/([0]-[1]))*[2]*(TMath::Exp(-x*TMath::Log(2)/[0]) - TMath::Exp(-x*TMath::Log(2)/[1]))+[3]*TMath::Exp(-x*TMath::Log(2)/[1])+ [4]*x + [5]",0,30000);

//2 half lives no background
TF1 *hl_feed_no = new TF1("hl_feed_no","([0]/([0]-[1]))*[2]*(TMath::Exp(-x*TMath::Log(2)/[0]) - TMath::Exp(-x*TMath::Log(2)/[1]))+[3]*TMath::Exp(-x*TMath::Log(2)/[1])",0,30000);



