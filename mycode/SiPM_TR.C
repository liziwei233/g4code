//
//* Environment: 
//**     class LZWfunc
//**     class DrawMyfunc
//* Data: 
//*      .root file obtaied by .trc files of Oscilloscope.
//* Function:
//**     Caculate the time resolution of SiPM test.
//* Date: 2019.3.19
//*
//***********************************
//
#include <string>
#include <time.h>
#include <TString.h>

#include "Include/DrawMyfunc.h"
#include "Include/LZWfunc.h"
#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>
#include <TH2.h>
#include <TF1.h>

using namespace std;

int SiPM_TR(const char* rootname="/mnt/f/SiPM/labtest/Scintillator/8-13/8-13source_bar2_v2.1_re"){
    gStyle->SetOptFit(1111);
    char name[100];
    char buff[1024];
    //const char* rootname="/mnt/f/SiPM/labtest/Scintillator/8-13/8-13source_bar2_v2.1_re";
    sprintf(name,"%s",rootname);
        sprintf(buff,"%s.root",name);
        TFile *f = new TFile(buff,"READ");
        TTree *t3 = (TTree*)f->Get("Pico");

    EVENT A,B;
    vector<EVENT*> ch; 
    ch.push_back(&A);
    ch.push_back(&B);
    //ch.push_back(&MCP);

    double frac[20];
    double thrd[20];
    

	//t3->SetBranchAddress("MCP2_all_charge",&ch.at(0).Q);
	//t3->SetBranchAddress("MCP2_all_charge",&A.Q);
 
	t3->SetBranchAddress("MCP1_all_charge",A.charge);
	t3->SetBranchAddress("MCP2_all_charge",B.charge);
	//t3->SetBranchAddress("MCP3_all_charge",&B.Q);
 
    t3->SetBranchAddress("MCP1_CFDtime",A.CFD);
	t3->SetBranchAddress("MCP2_CFDtime",B.CFD);
	//t3->SetBranchAddress("MCP3_CFDtime",B.CFD);
    
    t3->SetBranchAddress("MCP1_CFDfrac",frac);


	t3->SetBranchAddress("MCP1_global_maximum_y",&A.Amp);
	t3->SetBranchAddress("MCP2_global_maximum_y",&B.Amp);
	//t3->SetBranchAddress("MCP3_global_maximum_y",&B.Amp);

	t3->SetBranchAddress("MCP1_rise_time",&A.rise);
	t3->SetBranchAddress("MCP2_rise_time",&B.rise);
	//t3->SetBranchAddress("MCP3_rise_time",&B.rise);
	
	t3->SetBranchAddress("MCP1_baseline_rms",&A.blrms);
	t3->SetBranchAddress("MCP2_baseline_rms",&B.blrms);
	//t3->SetBranchAddress("MCP3_baseline_rms",&B.blrms);

//*
//*** Set your parameter ***
	//** Hist range !! **
    RANGE def={-9999,9999};
	RANGE x={0,200};
	//RANGE y={-5e-3,100e-3};
	RANGE y={-5e-3,500e-3};
    //RANGE q={-0.1,1};
    RANGE q={0,400};
	RANGE r={0,10};
    RANGE t={-5,5};
	int rbt=1,rbU=1;
    double fac=1.8;
    int iter=1;

    charRANGE ARange={
		x,
        y,
        q,
        r,
        def,
        def,
        t
	};//x,y,q,r,bl,blrms,t
	charRANGE BRange=ARange;
    //BRange.q.R=350;
    vector<charRANGE> range;
    range.push_back(ARange);
    range.push_back(BRange);
	//**
	//** CUT RANGE !!**
    charRANGE Acut={
		def,
        def,
        def,
        def,
        def,
        def,
        def
	};//blrms,y 
    charRANGE Bcut=Acut;
    vector<charRANGE> cut;
    cut.push_back(Acut);
    cut.push_back(Bcut); 

	OPTION opt={rbt,rbU,iter,fac};
//*
//**** end ****
//*
    

    LZWfunc TRfunc;
    DrawMyfunc draw;


	//TRfunc.Set_name(buff);
    //double L=sizeof(frac) / sizeof(frac[0]);
    //cout<<"the size of frac "<<L<<endl;
    int peArray[]={60,90,120,150,180,210,240,270,300};
    int N = sizeof(peArray)/sizeof(peArray[0]);
    for(int i=0; i<N;i++){
        cut.at(0).q.L=peArray[i]*0.8;
        cut.at(1).q.L=peArray[i]*0.9;
    sprintf(buff,"%s_TR_thred%dpe",name,peArray[i]);
    TRfunc.CH2Correction(t3,ch,frac,range,cut,opt,buff);
    }
   
   
   //delete ch;

/*
    //return 0;
    sprintf(buff,"%s_MCP_AB",name);
	TRfunc.Set_name(buff);
    TRfunc.CH3Correction(t3,ch,frac,rbU,rbt,fac,4);
*/
    //U.clear();
    //cut.clear();
    //ch.clear();

 return 0;   


 
 
}

int main(){
    return SiPM_TR();
}