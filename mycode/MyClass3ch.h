//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 20 13:21:06 2019 by ROOT version 6.12/06
// from TTree Pico/Analysis Output
// found on file: 201901-A3.root
//////////////////////////////////////////////////////////

#ifndef MyClass3ch_h
#define MyClass3ch_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class MyClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Baseline_Window;
   Double_t        MCP1_baseline_level;
   Double_t        MCP1_baseline_rms;
   Double_t        MCP1_global_maximum_y;
   Double_t        MCP1_global_maximum_x;
   Double_t        MCP1_start_x;
   Double_t        MCP1_end_x;
   Double_t        MCP1_invert_maximum_x;
   Double_t        MCP1_invert_maximum_y;
   Double_t        MCP1_all_charge[4];
   Double_t        MCP1_rise_time;
   Double_t        MCP1_CFDtime[8];
   Double_t        MCP1_CFDfrac[8];
   Bool_t          MCP1_CFDfailed[8];
   Double_t        MCP1_LEDtime[14];
   Double_t        MCP1_LEDthrd[14];
   Bool_t          MCP1_LEDfailed[14];
   Double_t        MCP2_baseline_level;
   Double_t        MCP2_baseline_rms;
   Double_t        MCP2_global_maximum_y;
   Double_t        MCP2_global_maximum_x;
   Double_t        MCP2_start_x;
   Double_t        MCP2_end_x;
   Double_t        MCP2_invert_maximum_x;
   Double_t        MCP2_invert_maximum_y;
   Double_t        MCP2_all_charge[4];
   Double_t        MCP2_rise_time;
   Double_t        MCP2_CFDtime[8];
   Double_t        MCP2_CFDfrac[8];
   Bool_t          MCP2_CFDfailed[8];
   Double_t        MCP2_LEDtime[14];
   Double_t        MCP2_LEDthrd[14];
   Bool_t          MCP2_LEDfailed[14];
   Double_t        MCP3_baseline_level;
   Double_t        MCP3_baseline_rms;
   Double_t        MCP3_global_maximum_y;
   Double_t        MCP3_global_maximum_x;
   Double_t        MCP3_start_x;
   Double_t        MCP3_end_x;
   Double_t        MCP3_invert_maximum_x;
   Double_t        MCP3_invert_maximum_y;
   Double_t        MCP3_all_charge[4];
   Double_t        MCP3_rise_time;
   Double_t        MCP3_CFDtime[8];
   Double_t        MCP3_CFDfrac[8];
   Bool_t          MCP3_CFDfailed[8];
   Double_t        MCP3_LEDtime[14];
   Double_t        MCP3_LEDthrd[14];
   Bool_t          MCP3_LEDfailed[14];

   // List of branches
   TBranch        *b_Baseline_Window;   //!
   TBranch        *b_MCP1_baseline_level;   //!
   TBranch        *b_MCP1_baseline_rms;   //!
   TBranch        *b_MCP1_global_maximum_y;   //!
   TBranch        *b_MCP1_global_maximum_x;   //!
   TBranch        *b_MCP1_start_x;   //!
   TBranch        *b_MCP1_end_x;   //!
   TBranch        *b_MCP1_invert_maximum_x;   //!
   TBranch        *b_MCP1_invert_maximum_y;   //!
   TBranch        *b_MCP1_all_charge;   //!
   TBranch        *b_MCP1_rise_time;   //!
   TBranch        *b_MCP1_CFDtime;   //!
   TBranch        *b_MCP1_CFDfrac;   //!
   TBranch        *b_MCP1_CFDfailed;   //!
   TBranch        *b_MCP1_LEDtime;   //!
   TBranch        *b_MCP1_LEDthrd;   //!
   TBranch        *b_MCP1_LEDfailed;   //!
   TBranch        *b_MCP2_baseline_level;   //!
   TBranch        *b_MCP2_baseline_rms;   //!
   TBranch        *b_MCP2_global_maximum_y;   //!
   TBranch        *b_MCP2_global_maximum_x;   //!
   TBranch        *b_MCP2_start_x;   //!
   TBranch        *b_MCP2_end_x;   //!
   TBranch        *b_MCP2_invert_maximum_x;   //!
   TBranch        *b_MCP2_invert_maximum_y;   //!
   TBranch        *b_MCP2_all_charge;   //!
   TBranch        *b_MCP2_rise_time;   //!
   TBranch        *b_MCP2_CFDtime;   //!
   TBranch        *b_MCP2_CFDfrac;   //!
   TBranch        *b_MCP2_CFDfailed;   //!
   TBranch        *b_MCP2_LEDtime;   //!
   TBranch        *b_MCP2_LEDthrd;   //!
   TBranch        *b_MCP2_LEDfailed;   //!
   TBranch        *b_MCP3_baseline_level;   //!
   TBranch        *b_MCP3_baseline_rms;   //!
   TBranch        *b_MCP3_global_maximum_y;   //!
   TBranch        *b_MCP3_global_maximum_x;   //!
   TBranch        *b_MCP3_start_x;   //!
   TBranch        *b_MCP3_end_x;   //!
   TBranch        *b_MCP3_invert_maximum_x;   //!
   TBranch        *b_MCP3_invert_maximum_y;   //!
   TBranch        *b_MCP3_all_charge;   //!
   TBranch        *b_MCP3_rise_time;   //!
   TBranch        *b_MCP3_CFDtime;   //!
   TBranch        *b_MCP3_CFDfrac;   //!
   TBranch        *b_MCP3_CFDfailed;   //!
   TBranch        *b_MCP3_LEDtime;   //!
   TBranch        *b_MCP3_LEDthrd;   //!
   TBranch        *b_MCP3_LEDfailed;   //!

   MyClass(const char *rootname = "", TTree *tree = 0);
   virtual ~MyClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyClass_cxx
MyClass::MyClass(const char *rootname, TTree *tree) : fChain(0)
{
   // if parameter tree is not specified (or zero), connect the file
   // used to generate this class and read the Tree.
   if (tree == 0)
   {
      TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject(rootname);
      if (!f || !f->IsOpen())
      {
         f = new TFile(rootname);
      }
      f->GetObject("Pico", tree);
   }
   Init(tree);
}

MyClass::~MyClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyClass::LoadTree(Long64_t entry)
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

void MyClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Baseline_Window", &Baseline_Window, &b_Baseline_Window);
   fChain->SetBranchAddress("MCP1_baseline_level", &MCP1_baseline_level, &b_MCP1_baseline_level);
   fChain->SetBranchAddress("MCP1_baseline_rms", &MCP1_baseline_rms, &b_MCP1_baseline_rms);
   fChain->SetBranchAddress("MCP1_global_maximum_y", &MCP1_global_maximum_y, &b_MCP1_global_maximum_y);
   fChain->SetBranchAddress("MCP1_global_maximum_x", &MCP1_global_maximum_x, &b_MCP1_global_maximum_x);
   fChain->SetBranchAddress("MCP1_start_x", &MCP1_start_x, &b_MCP1_start_x);
   fChain->SetBranchAddress("MCP1_end_x", &MCP1_end_x, &b_MCP1_end_x);
   fChain->SetBranchAddress("MCP1_invert_maximum_x", &MCP1_invert_maximum_x, &b_MCP1_invert_maximum_x);
   fChain->SetBranchAddress("MCP1_invert_maximum_y", &MCP1_invert_maximum_y, &b_MCP1_invert_maximum_y);
   fChain->SetBranchAddress("MCP1_all_charge", MCP1_all_charge, &b_MCP1_all_charge);
   fChain->SetBranchAddress("MCP1_rise_time", &MCP1_rise_time, &b_MCP1_rise_time);
   fChain->SetBranchAddress("MCP1_CFDtime", MCP1_CFDtime, &b_MCP1_CFDtime);
   fChain->SetBranchAddress("MCP1_CFDfrac", MCP1_CFDfrac, &b_MCP1_CFDfrac);
   fChain->SetBranchAddress("MCP1_CFDfailed", MCP1_CFDfailed, &b_MCP1_CFDfailed);
   fChain->SetBranchAddress("MCP1_LEDtime", MCP1_LEDtime, &b_MCP1_LEDtime);
   fChain->SetBranchAddress("MCP1_LEDthrd", MCP1_LEDthrd, &b_MCP1_LEDthrd);
   fChain->SetBranchAddress("MCP1_LEDfailed", MCP1_LEDfailed, &b_MCP1_LEDfailed);
   fChain->SetBranchAddress("MCP2_baseline_level", &MCP2_baseline_level, &b_MCP2_baseline_level);
   fChain->SetBranchAddress("MCP2_baseline_rms", &MCP2_baseline_rms, &b_MCP2_baseline_rms);
   fChain->SetBranchAddress("MCP2_global_maximum_y", &MCP2_global_maximum_y, &b_MCP2_global_maximum_y);
   fChain->SetBranchAddress("MCP2_global_maximum_x", &MCP2_global_maximum_x, &b_MCP2_global_maximum_x);
   fChain->SetBranchAddress("MCP2_start_x", &MCP2_start_x, &b_MCP2_start_x);
   fChain->SetBranchAddress("MCP2_end_x", &MCP2_end_x, &b_MCP2_end_x);
   fChain->SetBranchAddress("MCP2_invert_maximum_x", &MCP2_invert_maximum_x, &b_MCP2_invert_maximum_x);
   fChain->SetBranchAddress("MCP2_invert_maximum_y", &MCP2_invert_maximum_y, &b_MCP2_invert_maximum_y);
   fChain->SetBranchAddress("MCP2_all_charge", MCP2_all_charge, &b_MCP2_all_charge);
   fChain->SetBranchAddress("MCP2_rise_time", &MCP2_rise_time, &b_MCP2_rise_time);
   fChain->SetBranchAddress("MCP2_CFDtime", MCP2_CFDtime, &b_MCP2_CFDtime);
   fChain->SetBranchAddress("MCP2_CFDfrac", MCP2_CFDfrac, &b_MCP2_CFDfrac);
   fChain->SetBranchAddress("MCP2_CFDfailed", MCP2_CFDfailed, &b_MCP2_CFDfailed);
   fChain->SetBranchAddress("MCP2_LEDtime", MCP2_LEDtime, &b_MCP2_LEDtime);
   fChain->SetBranchAddress("MCP2_LEDthrd", MCP2_LEDthrd, &b_MCP2_LEDthrd);
   fChain->SetBranchAddress("MCP2_LEDfailed", MCP2_LEDfailed, &b_MCP2_LEDfailed);
   fChain->SetBranchAddress("MCP3_baseline_level", &MCP3_baseline_level, &b_MCP3_baseline_level);
   fChain->SetBranchAddress("MCP3_baseline_rms", &MCP3_baseline_rms, &b_MCP3_baseline_rms);
   fChain->SetBranchAddress("MCP3_global_maximum_y", &MCP3_global_maximum_y, &b_MCP3_global_maximum_y);
   fChain->SetBranchAddress("MCP3_global_maximum_x", &MCP3_global_maximum_x, &b_MCP3_global_maximum_x);
   fChain->SetBranchAddress("MCP3_start_x", &MCP3_start_x, &b_MCP3_start_x);
   fChain->SetBranchAddress("MCP3_end_x", &MCP3_end_x, &b_MCP3_end_x);
   fChain->SetBranchAddress("MCP3_invert_maximum_x", &MCP3_invert_maximum_x, &b_MCP3_invert_maximum_x);
   fChain->SetBranchAddress("MCP3_invert_maximum_y", &MCP3_invert_maximum_y, &b_MCP3_invert_maximum_y);
   fChain->SetBranchAddress("MCP3_all_charge", MCP3_all_charge, &b_MCP3_all_charge);
   fChain->SetBranchAddress("MCP3_rise_time", &MCP3_rise_time, &b_MCP3_rise_time);
   fChain->SetBranchAddress("MCP3_CFDtime", MCP3_CFDtime, &b_MCP3_CFDtime);
   fChain->SetBranchAddress("MCP3_CFDfrac", MCP3_CFDfrac, &b_MCP3_CFDfrac);
   fChain->SetBranchAddress("MCP3_CFDfailed", MCP3_CFDfailed, &b_MCP3_CFDfailed);
   fChain->SetBranchAddress("MCP3_LEDtime", MCP3_LEDtime, &b_MCP3_LEDtime);
   fChain->SetBranchAddress("MCP3_LEDthrd", MCP3_LEDthrd, &b_MCP3_LEDthrd);
   fChain->SetBranchAddress("MCP3_LEDfailed", MCP3_LEDfailed, &b_MCP3_LEDfailed);
   Notify();
}

Bool_t MyClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyClass_cxx
