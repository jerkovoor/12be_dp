//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 30 12:05:27 2019 by ROOT version 6.19/01
// from TTree AutoTree/AutoTree
// found on file: ../../Analysis/Be_newdecode/decodeBe_Yupedestal_5021.root
//////////////////////////////////////////////////////////

#ifndef Reduce_h
#define Reduce_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

using namespace std;

class Reduce {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           YdMul;
   vector<int>     *YdChannel;
   vector<double>  *YdEnergyRaw;
   vector<double>  *YdEnergy;
   vector<int>     *YdRing;
   vector<int>     *YdSector;
   Int_t           YuMul;
   vector<int>     *YuChannel;
   vector<double>  *YuEnergyRaw;
   vector<double>  *YuEnergy;
   vector<int>     *YuRing;
   vector<int>     *YuSector;
   Int_t           CsI1Mul;
   vector<int>     *CsI1Channel;
   vector<double>  *CsI1EnergyRaw;
   Int_t           CsI2Mul;
   vector<int>     *CsI2Channel;
   vector<double>  *CsI2EnergyRaw;
   Int_t           ICChannel;
   Double_t        ICEnergyRaw;
   Int_t           Sd1rMul;
   vector<int>     *Sd1rChannel;
   vector<double>  *Sd1rEnergyRaw;
   vector<double>  *Sd1rEnergy;
   Int_t           Sd1sMul;
   vector<int>     *Sd1sChannel;
   vector<double>  *Sd1sEnergyRaw;
   vector<double>  *Sd1sEnergy;
   Int_t           Sd2rMul;
   vector<int>     *Sd2rChannel;
   vector<double>  *Sd2rEnergyRaw;
   vector<double>  *Sd2rEnergy;
   Int_t           Sd2sMul;
   vector<int>     *Sd2sChannel;
   vector<double>  *Sd2sEnergyRaw;
   vector<double>  *Sd2sEnergy;
   Int_t           SurMul;
   vector<int>     *SurChannel;
   vector<double>  *SurEnergyRaw;
   Int_t           SusMul;
   vector<int>     *SusChannel;
   vector<double>  *SusEnergyRaw;

   // List of branches
   TBranch        *b_YdMul;   //!
   TBranch        *b_YdChannel;   //!
   TBranch        *b_YdEnergyRaw;   //!
   TBranch        *b_YdEnergy;   //!
   TBranch        *b_YdRing;   //!
   TBranch        *b_YdSector;   //!
   TBranch        *b_YuMul;   //!
   TBranch        *b_YuChannel;   //!
   TBranch        *b_YuEnergyRaw;   //!
   TBranch        *b_YuEnergy;   //!
   TBranch        *b_YuRing;   //!
   TBranch        *b_YuSector;   //!
   TBranch        *b_CsI1Mul;   //!
   TBranch        *b_CsI1Channel;   //!
   TBranch        *b_CsI1EnergyRaw;   //!
   TBranch        *b_CsI2Mul;   //!
   TBranch        *b_CsI2Channel;   //!
   TBranch        *b_CsI2EnergyRaw;   //!
   TBranch        *b_ICChannel;   //!
   TBranch        *b_ICEnergyRaw;   //!
   TBranch        *b_Sd1rMul;   //!
   TBranch        *b_Sd1rChannel;   //!
   TBranch        *b_Sd1rEnergyRaw;   //!
   TBranch        *b_Sd1rEnergy;   //!
   TBranch        *b_Sd1sMul;   //!
   TBranch        *b_Sd1sChannel;   //!
   TBranch        *b_Sd1sEnergyRaw;   //!
   TBranch        *b_Sd1sEnergy;   //!
   TBranch        *b_Sd2rMul;   //!
   TBranch        *b_Sd2rChannel;   //!
   TBranch        *b_Sd2rEnergyRaw;   //!
   TBranch        *b_Sd2rEnergy;   //!
   TBranch        *b_Sd2sMul;   //!
   TBranch        *b_Sd2sChannel;   //!
   TBranch        *b_Sd2sEnergyRaw;   //!
   TBranch        *b_Sd2sEnergy;   //!
   TBranch        *b_SurMul;   //!
   TBranch        *b_SurChannel;   //!
   TBranch        *b_SurEnergyRaw;   //!
   TBranch        *b_SusMul;   //!
   TBranch        *b_SusChannel;   //!
   TBranch        *b_SusEnergyRaw;   //!

   Reduce(TTree *tree=0);
   virtual ~Reduce();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Reduce_cxx
Reduce::Reduce(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../Analysis/Be_newdecode/decodeBe_Yupedestal_5021.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../Analysis/Be_newdecode/decodeBe_Yupedestal_5021.root");
      }
      f->GetObject("AutoTree",tree);

   }
   Init(tree);
}

Reduce::~Reduce()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Reduce::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Reduce::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Reduce::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   YdChannel = 0;
   YdEnergyRaw = 0;
   YdEnergy = 0;
   YdRing = 0;
   YdSector = 0;
   YuChannel = 0;
   YuEnergyRaw = 0;
   YuEnergy = 0;
   YuRing = 0;
   YuSector = 0;
   CsI1Channel = 0;
   CsI1EnergyRaw = 0;
   CsI2Channel = 0;
   CsI2EnergyRaw = 0;
   Sd1rChannel = 0;
   Sd1rEnergyRaw = 0;
   Sd1rEnergy = 0;
   Sd1sChannel = 0;
   Sd1sEnergyRaw = 0;
   Sd1sEnergy = 0;
   Sd2rChannel = 0;
   Sd2rEnergyRaw = 0;
   Sd2rEnergy = 0;
   Sd2sChannel = 0;
   Sd2sEnergyRaw = 0;
   Sd2sEnergy = 0;
   SurChannel = 0;
   SurEnergyRaw = 0;
   SusChannel = 0;
   SusEnergyRaw = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("YdMul", &YdMul, &b_YdMul);
   fChain->SetBranchAddress("YdChannel", &YdChannel, &b_YdChannel);
   fChain->SetBranchAddress("YdEnergyRaw", &YdEnergyRaw, &b_YdEnergyRaw);
   fChain->SetBranchAddress("YdEnergy", &YdEnergy, &b_YdEnergy);
   fChain->SetBranchAddress("YdRing", &YdRing, &b_YdRing);
   fChain->SetBranchAddress("YdSector", &YdSector, &b_YdSector);
   fChain->SetBranchAddress("YuMul", &YuMul, &b_YuMul);
   fChain->SetBranchAddress("YuChannel", &YuChannel, &b_YuChannel);
   fChain->SetBranchAddress("YuEnergyRaw", &YuEnergyRaw, &b_YuEnergyRaw);
   fChain->SetBranchAddress("YuEnergy", &YuEnergy, &b_YuEnergy);
   fChain->SetBranchAddress("YuRing", &YuRing, &b_YuRing);
   fChain->SetBranchAddress("YuSector", &YuSector, &b_YuSector);
   fChain->SetBranchAddress("CsI1Mul", &CsI1Mul, &b_CsI1Mul);
   fChain->SetBranchAddress("CsI1Channel", &CsI1Channel, &b_CsI1Channel);
   fChain->SetBranchAddress("CsI1EnergyRaw", &CsI1EnergyRaw, &b_CsI1EnergyRaw);
   fChain->SetBranchAddress("CsI2Mul", &CsI2Mul, &b_CsI2Mul);
   fChain->SetBranchAddress("CsI2Channel", &CsI2Channel, &b_CsI2Channel);
   fChain->SetBranchAddress("CsI2EnergyRaw", &CsI2EnergyRaw, &b_CsI2EnergyRaw);
   fChain->SetBranchAddress("ICChannel", &ICChannel, &b_ICChannel);
   fChain->SetBranchAddress("ICEnergyRaw", &ICEnergyRaw, &b_ICEnergyRaw);
   fChain->SetBranchAddress("Sd1rMul", &Sd1rMul, &b_Sd1rMul);
   fChain->SetBranchAddress("Sd1rChannel", &Sd1rChannel, &b_Sd1rChannel);
   fChain->SetBranchAddress("Sd1rEnergyRaw", &Sd1rEnergyRaw, &b_Sd1rEnergyRaw);
   fChain->SetBranchAddress("Sd1rEnergy", &Sd1rEnergy, &b_Sd1rEnergy);
   fChain->SetBranchAddress("Sd1sMul", &Sd1sMul, &b_Sd1sMul);
   fChain->SetBranchAddress("Sd1sChannel", &Sd1sChannel, &b_Sd1sChannel);
   fChain->SetBranchAddress("Sd1sEnergyRaw", &Sd1sEnergyRaw, &b_Sd1sEnergyRaw);
   fChain->SetBranchAddress("Sd1sEnergy", &Sd1sEnergy, &b_Sd1sEnergy);
   fChain->SetBranchAddress("Sd2rMul", &Sd2rMul, &b_Sd2rMul);
   fChain->SetBranchAddress("Sd2rChannel", &Sd2rChannel, &b_Sd2rChannel);
   fChain->SetBranchAddress("Sd2rEnergyRaw", &Sd2rEnergyRaw, &b_Sd2rEnergyRaw);
   fChain->SetBranchAddress("Sd2rEnergy", &Sd2rEnergy, &b_Sd2rEnergy);
   fChain->SetBranchAddress("Sd2sMul", &Sd2sMul, &b_Sd2sMul);
   fChain->SetBranchAddress("Sd2sChannel", &Sd2sChannel, &b_Sd2sChannel);
   fChain->SetBranchAddress("Sd2sEnergyRaw", &Sd2sEnergyRaw, &b_Sd2sEnergyRaw);
   fChain->SetBranchAddress("Sd2sEnergy", &Sd2sEnergy, &b_Sd2sEnergy);
   fChain->SetBranchAddress("SurMul", &SurMul, &b_SurMul);
   fChain->SetBranchAddress("SurChannel", &SurChannel, &b_SurChannel);
   fChain->SetBranchAddress("SurEnergyRaw", &SurEnergyRaw, &b_SurEnergyRaw);
   fChain->SetBranchAddress("SusMul", &SusMul, &b_SusMul);
   fChain->SetBranchAddress("SusChannel", &SusChannel, &b_SusChannel);
   fChain->SetBranchAddress("SusEnergyRaw", &SusEnergyRaw, &b_SusEnergyRaw);
   Notify();
}

Bool_t Reduce::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Reduce::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Reduce::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Reduce_cxx
