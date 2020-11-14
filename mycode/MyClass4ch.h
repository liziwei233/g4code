//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 19 19:22:53 2020 by ROOT version 6.12/06
// from TTree Pico/Analysis Output
// found on file: KT0864-2150V-16Anode-A1A2A4.root
//////////////////////////////////////////////////////////

#ifndef MyClass_h
#define MyClass_h

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
   Double_t        TR1_baseline_level;
   Double_t        TR1_baseline_rms;
   Double_t        TR1_global_maximum_y;
   Double_t        TR1_global_maximum_x;
   Double_t        TR1_start_x;
   Double_t        TR1_end_x;
   Double_t        TR1_invert_maximum_x;
   Double_t        TR1_invert_maximum_y;
   Double_t        TR1_secondinvertpeak_x;
   Double_t        TR1_secondinvertpeak_y;
   Double_t        TR1_all_charge[4];
   Double_t        TR1_rise_time[4];
   Double_t        TR1_width;
   Double_t        TR1_CFDtime[8];
   Double_t        TR1_CFDfrac[8];
   Double_t        TR1_CFDr[8];
   Bool_t          TR1_CFDfailed[8];
   Double_t        MCP2_baseline_level;
   Double_t        MCP2_baseline_rms;
   Double_t        MCP2_global_maximum_y;
   Double_t        MCP2_global_maximum_x;
   Double_t        MCP2_start_x;
   Double_t        MCP2_end_x;
   Double_t        MCP2_invert_maximum_x;
   Double_t        MCP2_invert_maximum_y;
   Double_t        MCP2_secondinvertpeak_x;
   Double_t        MCP2_secondinvertpeak_y;
   Double_t        MCP2_all_charge[4];
   Double_t        MCP2_rise_time[4];
   Double_t        MCP2_width;
   Double_t        MCP2_CFDtime[8];
   Double_t        MCP2_CFDfrac[8];
   Double_t        MCP2_CFDr[8];
   Bool_t          MCP2_CFDfailed[8];
   Double_t        MCP3_baseline_level;
   Double_t        MCP3_baseline_rms;
   Double_t        MCP3_global_maximum_y;
   Double_t        MCP3_global_maximum_x;
   Double_t        MCP3_start_x;
   Double_t        MCP3_end_x;
   Double_t        MCP3_invert_maximum_x;
   Double_t        MCP3_invert_maximum_y;
   Double_t        MCP3_secondinvertpeak_x;
   Double_t        MCP3_secondinvertpeak_y;
   Double_t        MCP3_all_charge[4];
   Double_t        MCP3_rise_time[4];
   Double_t        MCP3_width;
   Double_t        MCP3_CFDtime[8];
   Double_t        MCP3_CFDfrac[8];
   Double_t        MCP3_CFDr[8];
   Bool_t          MCP3_CFDfailed[8];
   Double_t        MCP4_baseline_level;
   Double_t        MCP4_baseline_rms;
   Double_t        MCP4_global_maximum_y;
   Double_t        MCP4_global_maximum_x;
   Double_t        MCP4_start_x;
   Double_t        MCP4_end_x;
   Double_t        MCP4_invert_maximum_x;
   Double_t        MCP4_invert_maximum_y;
   Double_t        MCP4_secondinvertpeak_x;
   Double_t        MCP4_secondinvertpeak_y;
   Double_t        MCP4_all_charge[4];
   Double_t        MCP4_rise_time[4];
   Double_t        MCP4_width;
   Double_t        MCP4_CFDtime[8];
   Double_t        MCP4_CFDfrac[8];
   Double_t        MCP4_CFDr[8];
   Bool_t          MCP4_CFDfailed[8];

   // List of branches
   TBranch        *b_Baseline_Window;   //!
   TBranch        *b_TR1_baseline_level;   //!
   TBranch        *b_TR1_baseline_rms;   //!
   TBranch        *b_TR1_global_maximum_y;   //!
   TBranch        *b_TR1_global_maximum_x;   //!
   TBranch        *b_TR1_start_x;   //!
   TBranch        *b_TR1_end_x;   //!
   TBranch        *b_TR1_invert_maximum_x;   //!
   TBranch        *b_TR1_invert_maximum_y;   //!
   TBranch        *b_TR1_secondinvertpeak_x;   //!
   TBranch        *b_TR1_secondinvertpeak_y;   //!
   TBranch        *b_TR1_all_charge;   //!
   TBranch        *b_TR1_rise_time;   //!
   TBranch        *b_TR1_width;   //!
   TBranch        *b_TR1_CFDtime;   //!
   TBranch        *b_TR1_CFDfrac;   //!
   TBranch        *b_TR1_CFDr;   //!
   TBranch        *b_TR1_CFDfailed;   //!
   TBranch        *b_MCP2_baseline_level;   //!
   TBranch        *b_MCP2_baseline_rms;   //!
   TBranch        *b_MCP2_global_maximum_y;   //!
   TBranch        *b_MCP2_global_maximum_x;   //!
   TBranch        *b_MCP2_start_x;   //!
   TBranch        *b_MCP2_end_x;   //!
   TBranch        *b_MCP2_invert_maximum_x;   //!
   TBranch        *b_MCP2_invert_maximum_y;   //!
   TBranch        *b_MCP2_secondinvertpeak_x;   //!
   TBranch        *b_MCP2_secondinvertpeak_y;   //!
   TBranch        *b_MCP2_all_charge;   //!
   TBranch        *b_MCP2_rise_time;   //!
   TBranch        *b_MCP2_width;   //!
   TBranch        *b_MCP2_CFDtime;   //!
   TBranch        *b_MCP2_CFDfrac;   //!
   TBranch        *b_MCP2_CFDr;   //!
   TBranch        *b_MCP2_CFDfailed;   //!
   TBranch        *b_MCP3_baseline_level;   //!
   TBranch        *b_MCP3_baseline_rms;   //!
   TBranch        *b_MCP3_global_maximum_y;   //!
   TBranch        *b_MCP3_global_maximum_x;   //!
   TBranch        *b_MCP3_start_x;   //!
   TBranch        *b_MCP3_end_x;   //!
   TBranch        *b_MCP3_invert_maximum_x;   //!
   TBranch        *b_MCP3_invert_maximum_y;   //!
   TBranch        *b_MCP3_secondinvertpeak_x;   //!
   TBranch        *b_MCP3_secondinvertpeak_y;   //!
   TBranch        *b_MCP3_all_charge;   //!
   TBranch        *b_MCP3_rise_time;   //!
   TBranch        *b_MCP3_width;   //!
   TBranch        *b_MCP3_CFDtime;   //!
   TBranch        *b_MCP3_CFDfrac;   //!
   TBranch        *b_MCP3_CFDr;   //!
   TBranch        *b_MCP3_CFDfailed;   //!
   TBranch        *b_MCP4_baseline_level;   //!
   TBranch        *b_MCP4_baseline_rms;   //!
   TBranch        *b_MCP4_global_maximum_y;   //!
   TBranch        *b_MCP4_global_maximum_x;   //!
   TBranch        *b_MCP4_start_x;   //!
   TBranch        *b_MCP4_end_x;   //!
   TBranch        *b_MCP4_invert_maximum_x;   //!
   TBranch        *b_MCP4_invert_maximum_y;   //!
   TBranch        *b_MCP4_secondinvertpeak_x;   //!
   TBranch        *b_MCP4_secondinvertpeak_y;   //!
   TBranch        *b_MCP4_all_charge;   //!
   TBranch        *b_MCP4_rise_time;   //!
   TBranch        *b_MCP4_width;   //!
   TBranch        *b_MCP4_CFDtime;   //!
   TBranch        *b_MCP4_CFDfrac;   //!
   TBranch        *b_MCP4_CFDr;   //!
   TBranch        *b_MCP4_CFDfailed;   //!

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
   fChain->SetBranchAddress("TR1_baseline_level", &TR1_baseline_level, &b_TR1_baseline_level);
   fChain->SetBranchAddress("TR1_baseline_rms", &TR1_baseline_rms, &b_TR1_baseline_rms);
   fChain->SetBranchAddress("TR1_global_maximum_y", &TR1_global_maximum_y, &b_TR1_global_maximum_y);
   fChain->SetBranchAddress("TR1_global_maximum_x", &TR1_global_maximum_x, &b_TR1_global_maximum_x);
   fChain->SetBranchAddress("TR1_start_x", &TR1_start_x, &b_TR1_start_x);
   fChain->SetBranchAddress("TR1_end_x", &TR1_end_x, &b_TR1_end_x);
   fChain->SetBranchAddress("TR1_invert_maximum_x", &TR1_invert_maximum_x, &b_TR1_invert_maximum_x);
   fChain->SetBranchAddress("TR1_invert_maximum_y", &TR1_invert_maximum_y, &b_TR1_invert_maximum_y);
   fChain->SetBranchAddress("TR1_secondinvertpeak_x", &TR1_secondinvertpeak_x, &b_TR1_secondinvertpeak_x);
   fChain->SetBranchAddress("TR1_secondinvertpeak_y", &TR1_secondinvertpeak_y, &b_TR1_secondinvertpeak_y);
   fChain->SetBranchAddress("TR1_all_charge", TR1_all_charge, &b_TR1_all_charge);
   fChain->SetBranchAddress("TR1_rise_time", TR1_rise_time, &b_TR1_rise_time);
   fChain->SetBranchAddress("TR1_width", &TR1_width, &b_TR1_width);
   fChain->SetBranchAddress("TR1_CFDtime", TR1_CFDtime, &b_TR1_CFDtime);
   fChain->SetBranchAddress("TR1_CFDfrac", TR1_CFDfrac, &b_TR1_CFDfrac);
   fChain->SetBranchAddress("TR1_CFDr", TR1_CFDr, &b_TR1_CFDr);
   fChain->SetBranchAddress("TR1_CFDfailed", TR1_CFDfailed, &b_TR1_CFDfailed);
   fChain->SetBranchAddress("MCP2_baseline_level", &MCP2_baseline_level, &b_MCP2_baseline_level);
   fChain->SetBranchAddress("MCP2_baseline_rms", &MCP2_baseline_rms, &b_MCP2_baseline_rms);
   fChain->SetBranchAddress("MCP2_global_maximum_y", &MCP2_global_maximum_y, &b_MCP2_global_maximum_y);
   fChain->SetBranchAddress("MCP2_global_maximum_x", &MCP2_global_maximum_x, &b_MCP2_global_maximum_x);
   fChain->SetBranchAddress("MCP2_start_x", &MCP2_start_x, &b_MCP2_start_x);
   fChain->SetBranchAddress("MCP2_end_x", &MCP2_end_x, &b_MCP2_end_x);
   fChain->SetBranchAddress("MCP2_invert_maximum_x", &MCP2_invert_maximum_x, &b_MCP2_invert_maximum_x);
   fChain->SetBranchAddress("MCP2_invert_maximum_y", &MCP2_invert_maximum_y, &b_MCP2_invert_maximum_y);
   fChain->SetBranchAddress("MCP2_secondinvertpeak_x", &MCP2_secondinvertpeak_x, &b_MCP2_secondinvertpeak_x);
   fChain->SetBranchAddress("MCP2_secondinvertpeak_y", &MCP2_secondinvertpeak_y, &b_MCP2_secondinvertpeak_y);
   fChain->SetBranchAddress("MCP2_all_charge", MCP2_all_charge, &b_MCP2_all_charge);
   fChain->SetBranchAddress("MCP2_rise_time", MCP2_rise_time, &b_MCP2_rise_time);
   fChain->SetBranchAddress("MCP2_width", &MCP2_width, &b_MCP2_width);
   fChain->SetBranchAddress("MCP2_CFDtime", MCP2_CFDtime, &b_MCP2_CFDtime);
   fChain->SetBranchAddress("MCP2_CFDfrac", MCP2_CFDfrac, &b_MCP2_CFDfrac);
   fChain->SetBranchAddress("MCP2_CFDr", MCP2_CFDr, &b_MCP2_CFDr);
   fChain->SetBranchAddress("MCP2_CFDfailed", MCP2_CFDfailed, &b_MCP2_CFDfailed);
   fChain->SetBranchAddress("MCP3_baseline_level", &MCP3_baseline_level, &b_MCP3_baseline_level);
   fChain->SetBranchAddress("MCP3_baseline_rms", &MCP3_baseline_rms, &b_MCP3_baseline_rms);
   fChain->SetBranchAddress("MCP3_global_maximum_y", &MCP3_global_maximum_y, &b_MCP3_global_maximum_y);
   fChain->SetBranchAddress("MCP3_global_maximum_x", &MCP3_global_maximum_x, &b_MCP3_global_maximum_x);
   fChain->SetBranchAddress("MCP3_start_x", &MCP3_start_x, &b_MCP3_start_x);
   fChain->SetBranchAddress("MCP3_end_x", &MCP3_end_x, &b_MCP3_end_x);
   fChain->SetBranchAddress("MCP3_invert_maximum_x", &MCP3_invert_maximum_x, &b_MCP3_invert_maximum_x);
   fChain->SetBranchAddress("MCP3_invert_maximum_y", &MCP3_invert_maximum_y, &b_MCP3_invert_maximum_y);
   fChain->SetBranchAddress("MCP3_secondinvertpeak_x", &MCP3_secondinvertpeak_x, &b_MCP3_secondinvertpeak_x);
   fChain->SetBranchAddress("MCP3_secondinvertpeak_y", &MCP3_secondinvertpeak_y, &b_MCP3_secondinvertpeak_y);
   fChain->SetBranchAddress("MCP3_all_charge", MCP3_all_charge, &b_MCP3_all_charge);
   fChain->SetBranchAddress("MCP3_rise_time", MCP3_rise_time, &b_MCP3_rise_time);
   fChain->SetBranchAddress("MCP3_width", &MCP3_width, &b_MCP3_width);
   fChain->SetBranchAddress("MCP3_CFDtime", MCP3_CFDtime, &b_MCP3_CFDtime);
   fChain->SetBranchAddress("MCP3_CFDfrac", MCP3_CFDfrac, &b_MCP3_CFDfrac);
   fChain->SetBranchAddress("MCP3_CFDr", MCP3_CFDr, &b_MCP3_CFDr);
   fChain->SetBranchAddress("MCP3_CFDfailed", MCP3_CFDfailed, &b_MCP3_CFDfailed);
   fChain->SetBranchAddress("MCP4_baseline_level", &MCP4_baseline_level, &b_MCP4_baseline_level);
   fChain->SetBranchAddress("MCP4_baseline_rms", &MCP4_baseline_rms, &b_MCP4_baseline_rms);
   fChain->SetBranchAddress("MCP4_global_maximum_y", &MCP4_global_maximum_y, &b_MCP4_global_maximum_y);
   fChain->SetBranchAddress("MCP4_global_maximum_x", &MCP4_global_maximum_x, &b_MCP4_global_maximum_x);
   fChain->SetBranchAddress("MCP4_start_x", &MCP4_start_x, &b_MCP4_start_x);
   fChain->SetBranchAddress("MCP4_end_x", &MCP4_end_x, &b_MCP4_end_x);
   fChain->SetBranchAddress("MCP4_invert_maximum_x", &MCP4_invert_maximum_x, &b_MCP4_invert_maximum_x);
   fChain->SetBranchAddress("MCP4_invert_maximum_y", &MCP4_invert_maximum_y, &b_MCP4_invert_maximum_y);
   fChain->SetBranchAddress("MCP4_secondinvertpeak_x", &MCP4_secondinvertpeak_x, &b_MCP4_secondinvertpeak_x);
   fChain->SetBranchAddress("MCP4_secondinvertpeak_y", &MCP4_secondinvertpeak_y, &b_MCP4_secondinvertpeak_y);
   fChain->SetBranchAddress("MCP4_all_charge", MCP4_all_charge, &b_MCP4_all_charge);
   fChain->SetBranchAddress("MCP4_rise_time", MCP4_rise_time, &b_MCP4_rise_time);
   fChain->SetBranchAddress("MCP4_width", &MCP4_width, &b_MCP4_width);
   fChain->SetBranchAddress("MCP4_CFDtime", MCP4_CFDtime, &b_MCP4_CFDtime);
   fChain->SetBranchAddress("MCP4_CFDfrac", MCP4_CFDfrac, &b_MCP4_CFDfrac);
   fChain->SetBranchAddress("MCP4_CFDr", MCP4_CFDr, &b_MCP4_CFDr);
   fChain->SetBranchAddress("MCP4_CFDfailed", MCP4_CFDfailed, &b_MCP4_CFDfailed);
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
