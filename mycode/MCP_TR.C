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

int MCP_TR(){
    gStyle->SetOptFit(1111);
    char name[100];
    char buff[1024];
    const char* rootname="rawdata/th100";
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
 
	//t3->SetBranchAddress("MCP1_all_charge",&MCP.Q);
	t3->SetBranchAddress("MCP1_all_charge",&A.Q);
	t3->SetBranchAddress("MCP2_all_charge",&B.Q);
 
    //t3->SetBranchAddress("MCP1_CFDtime",MCP.CFD);
	t3->SetBranchAddress("MCP1_CFDtime",A.CFD);
	t3->SetBranchAddress("MCP2_CFDtime",B.CFD);
    
    t3->SetBranchAddress("MCP1_CFDfrac",frac);


	//t3->SetBranchAddress("MCP1_global_maximum_y",&MCP.Amp);
	t3->SetBranchAddress("MCP1_global_maximum_y",&A.Amp);
	t3->SetBranchAddress("MCP2_global_maximum_y",&B.Amp);

	//t3->SetBranchAddress("MCP1_rise_time",&MCP.rise);
	t3->SetBranchAddress("MCP1_rise_time",&A.rise);
	t3->SetBranchAddress("MCP2_rise_time",&B.rise);
	
	//t3->SetBranchAddress("MCP1_baseline_rms",&MCP.blrms);
	t3->SetBranchAddress("MCP1_baseline_rms",&A.blrms);
	t3->SetBranchAddress("MCP2_baseline_rms",&B.blrms);


    RANGE def={-9999,9999};
    //* set cut vector 
    //*
    vector<CUT> cut;
    //cout<<cut<<endl;
    //RANGE tcut={-0.8,-0.1};
    //CUT cutA={-1e4,0.99,0,1.2e3,0,5};    
    CUT cutA={
        def,
        def,
        def,
        def,
        def,
        def,
        def
    };    
    //CUT cutA={-1e4,1e4,0,1.2e3,0,5};    
    CUT cutB=cutA;    
    //CUT cutMCP={-1e4,1e4,0,20,0,5};
    cut.push_back(cutA);    
    cut.push_back(cutB);    

    //* set charactor range vector
    //*
    RANGE y={0,2.5};
    RANGE q={0,30};
    RANGE t={-1,1};
    vector<charRANGE> range; 
    charRANGE charA={
        def,
        y,
        q,
        def,
        def,
        def,
        t
    };
    charRANGE charB=charA;
    range.push_back(charA);
    range.push_back(charB);

    //Correct option set
    int rbt=16,rbU=8;
    int iter=4;
    double fac=2;
    OPTION opt={rbt,rbU,iter,fac};

    LZWfunc TRfunc;
    //TRfunc.Set_cut(cut);
    //TRfunc.Set_t(t);
    //TRfunc.Set_U(U);


    sprintf(buff,"%s_A-B",name);
	//TRfunc.Set_name(buff);
    //double L=sizeof(frac) / sizeof(frac[0]);
    //cout<<"the size of frac "<<L<<endl;
    TRfunc.CH2Correction(t3,ch,frac,range,cut,opt,name);
  
   
   
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
    return MCP_TR();
}