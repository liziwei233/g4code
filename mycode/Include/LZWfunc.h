
#ifndef __LZWfunc_h__
#define __LZWfunc_h__

#include "TCanvas.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "DrawMyfunc.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"


using namespace std;

ofstream output;
struct EVENT{
    double time,CFD[20],LED[20],Amp,Q,rise,bl,blrms,x;
    int npe;
};
struct CUT{
    double Amplow,Ampup;
    double Qlow,Qup;
    double riselow,riseup;
};
struct POSITION{
    double x;
    double y;
};
struct RANGE{
    double L;
    double R;
};
struct charRANGE{
    RANGE x,y,q,r,bl,blrms,t;
};


class LZWfunc
{
    public:
    LZWfunc();
    ~LZWfunc();
    

    Double_t fpeaks(Double_t *x, Double_t *par);
    TF1* gausfit(TH1 *h, int rbU, double fac, RANGE U);
    TF1* gausfit(TH1 *h, int rbU, double fac, RANGE* U);
    TF1* twoguasfit(TH1 *ht, double fac, int rbt, RANGE* t);
    TF1* twoguasfit(TH1 *ht, double fac, int rbt, RANGE t);
    TF1* profilefit(TH2 *Rt, double rbU, double rbt, RANGE t, RANGE U, char *name);
    TF1* CH3Correction(TTree *t1, EVENT *A, EVENT *B, EVENT *MCP, CUT basecut, CUT selcut, int rbU, int rbt, char *name, double fac, int iter, RANGE t, RANGE U);
    void CH3Correction(TTree *t1, vector<EVENT *> ch, double* p, int rbU, int rbt, double fac, int iter);
    TF1* CH2Correction(TTree *t1, EVENT *A, EVENT *B, EVENT *MCP, CUT basecut, CUT selcut, int rbU, int rbt, char *name, double fac, int iter, RANGE t, RANGE U);
    TF1* CH2Correction(TTree *t1, EVENT *A, EVENT *B, double Uth, int rbU, int rbt, char *name, double fac, int iter, RANGE t, RANGE U);
    TF1* CH2Correction(TTree * t1, EVENT * A, EVENT * B, POSITION * pos, double range, int rbU, int rbt, char *name, double fac, int iter, RANGE t,RANGE U);
    TF1* CH2Correction(TTree* t1,EVENT* A,EVENT* B,double* p,int rbU,int rbt,char* name,double fac,int iter, RANGE t, RANGE U);
    void CH2Correction(TTree* t1,vector<EVENT*> ch, double* p,int rbU,int rbt,double fac,int iter);
    
    TF1* CH1Correction(TTree * t1, EVENT * A, EVENT * B, double *p, int rbU, int rbt, char *name, double fac, int iter, RANGE t, RANGE U);
    void CH1Correction(TTree * t1, EVENT * A, double * t0, RANGE t,RANGE U, int rbU, int rbt,  double fac, int iter);

    void drawcharacter(TTree *t3, int chN, string* chname, vector<charRANGE> chR);

    string getTime();


    void Set_t(RANGE formal_t){
        t=formal_t;
    };
    void Set_U(vector<RANGE> formal_U){
        U=formal_U;
    };
    void Set_name(string formal_name){
        name=formal_name;
    };
    void Set_cut(vector<CUT> formal_cut){
        cut=formal_cut;
    };
    void Set_charcut(charRANGE formal_charcut){
        charcut=formal_charcut;
    };
    void Set_charRange(vector<charRANGE> formal_chR){
        chR=formal_chR;
    }
    

    RANGE Get_t(){
        return t;
    };
    RANGE Get_U(int Uid){
        return U.at(Uid);
    };
    string Get_name(){
        return name;
    };
    CUT Get_cut(int cutid){
        return cut.at(cutid);
    };

    charRANGE Get_charcut(){
        return charcut;
    }

    charRANGE Get_charRange(int chN){
        return chR.at(chN);
    };
   
    

    private:
    vector<CUT> cut;
    vector<RANGE> U;
    RANGE t;
    string name;
    charRANGE charcut;
    vector<charRANGE> chR;
    char buff[1024];
    /*TCanvas* c1;
    float MSize;
    float MStyle;
    Color_t MColor;
    Color_t LColor;
    float LWidth;
    int LStyle;
    Color_t FColor; */

};
#endif