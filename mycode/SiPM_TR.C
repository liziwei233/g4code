#include "/mnt/c/Subsys/work/g4code/mycode/Include/DrawMyfunc.h"
#include "/mnt/c/Subsys/work/g4code/mycode/Include/LZWfunc.h"
#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>
#include <TH2.h>
#include <TF1.h>

using namespace std;

int SiPM_TR(){
    gStyle->SetOptFit(1111);
    char name[100];
    char buff[1024];
    const char* rootname="/mnt/f/experiment/FTOF/SIPM/labtest/2019-1-14/cos";
    sprintf(name,"%s",rootname);
        sprintf(buff,"%s.root",name);
        TFile *f = new TFile(buff,"READ");
        TTree *t3 = (TTree*)f->Get("Pico");

    EVENT A,B,MCP;
    vector<EVENT*> ch; 
    ch.push_back(&A);
    ch.push_back(&B);
    ch.push_back(&MCP);

    double frac[20];
    double thrd[20];
    

	//t3->SetBranchAddress("MCP2_all_charge",&ch.at(0).Q);
	//t3->SetBranchAddress("MCP2_all_charge",&A.Q);
 
	t3->SetBranchAddress("MCP1_all_charge",&MCP.Q);
	t3->SetBranchAddress("MCP2_all_charge",&A.Q);
	t3->SetBranchAddress("MCP3_all_charge",&B.Q);
 
    t3->SetBranchAddress("MCP1_CFDtime",MCP.CFD);
	t3->SetBranchAddress("MCP2_CFDtime",A.CFD);
	t3->SetBranchAddress("MCP3_CFDtime",B.CFD);
    
    t3->SetBranchAddress("MCP3_CFDfrac",frac);


	t3->SetBranchAddress("MCP1_global_maximum_y",&MCP.Amp);
	t3->SetBranchAddress("MCP2_global_maximum_y",&A.Amp);
	t3->SetBranchAddress("MCP3_global_maximum_y",&B.Amp);

	t3->SetBranchAddress("MCP1_rise_time",&MCP.rise);
	t3->SetBranchAddress("MCP2_rise_time",&A.rise);
	t3->SetBranchAddress("MCP3_rise_time",&B.rise);
	
	t3->SetBranchAddress("MCP1_baseline_rms",&MCP.blrms);
	t3->SetBranchAddress("MCP2_baseline_rms",&A.blrms);
	t3->SetBranchAddress("MCP3_baseline_rms",&B.blrms);


    
    vector<CUT> cut;
    //CUT cutA={-1e4,0.99,0,1.2e3,0,5};    
    CUT cutA={0.9,1e4,0,1.2e3,0,5};    
    //CUT cutA={-1e4,1e4,0,1.2e3,0,5};    
    CUT cutB=cutA;    
    CUT cutMCP={-1e4,0.23,0,20,0,5};
    //CUT cutMCP={-1e4,1e4,0,20,0,5};
    cut.push_back(cutA);    
    cut.push_back(cutB);    
    cut.push_back(cutMCP);    

    RANGE t={-5,2};

    vector<RANGE> U; 
    RANGE UA={0,1200};
    RANGE UB=UA;
    RANGE UMCP={0,12};
    U.push_back(UA);
    U.push_back(UB);
    U.push_back(UMCP);

    //RANGE t(-2,2);
    int rbt=4,rbU=4;
    double fac=2;

    LZWfunc TRfunc;
    TRfunc.Set_cut(cut);
    TRfunc.Set_t(t);
    TRfunc.Set_U(U);


    sprintf(buff,"%s_A-B",name);
	TRfunc.Set_name(buff);
    //double L=sizeof(frac) / sizeof(frac[0]);
    //cout<<"the size of frac "<<L<<endl;
    TRfunc.CH2Correction(t3,ch,frac,rbU,rbt,fac,4);
  
   
   
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