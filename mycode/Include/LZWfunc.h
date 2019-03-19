
#ifndef __LZWfunc_h__
#define __LZWfunc_h__

#include "TTree.h"
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
#include "TPaveStats.h"
#include "TList.h"

using namespace std;
//*
//* if run for geant4 data, disable the notation below
#define G4_FLAG




ofstream output;







typedef struct POSITION{
    double x;
    double y;
    double z;
}DIRECTION;

struct RANGE{
    double L;
    double R;
};

struct OPTION{
    int rbt,rbu,iter;
    double fac;
};

struct gausPAR{
    double h; //height 
    double m; //mean
    double s;  //sigma
    bool operator < (const gausPAR & gP) const {//symbol overloading
            return m<gP.m;
    }
}; 

#ifndef G4_FLAG
struct EVENT{
    double time,CFD[20],LED[20],Amp,Q,rise,bl,blrms,x;
};
typedef struct charRANGE{
    RANGE x,y,q,r,bl,blrms,t;
}CUT;
#else
struct EVENT{
    double time,CFD[20],LED[20],Amp;
    int npe;
};
typedef struct charRANGE{
    RANGE y,t,npe;
}CUT;
#endif

class LZWfunc
{
    public:
    LZWfunc();
    ~LZWfunc();
    

    Double_t fpeaks(Double_t *x, Double_t *par);
    TF1* fpeaksfit(TH1 *ha,int npeaks,double res,double sigma, double thrd);
    TF1* gausfit(TH1 *h, int rbU, double fac, RANGE U);
    TF1* gausfit(TH1 *h, int rbU, double fac, RANGE* U);
    TF1* twoguasfit(TH1 *ht, double fac, int rbt, RANGE* t);
    TF1* twoguasfit(TH1 *ht, double fac, int rbt, RANGE t);
    TF1* profilefit(TH2 *Rt, double rbU, double rbt, RANGE t, RANGE U, char *name);
    
    void CH3Correction(TTree *t1, vector<EVENT *> ch, double *p, vector<charRANGE> range, vector<CUT> cut, OPTION opt, string name);
    TF1 *CH3Correction(TTree *t1, vector<EVENT *> ch, vector<charRANGE> range, vector<CUT> cut, OPTION opt, string name);
    TF1 *CH2Correction(TTree *t1, vector<EVENT *> ch, vector<charRANGE> range, vector<CUT> cut, OPTION opt, string name);
    void CH2Correction(TTree *t1, vector<EVENT *> ch, double *p, vector<charRANGE> range, vector<CUT> cut, OPTION opt, string name);
    TF1 *CH1Correction(TTree *t1, EVENT *A, charRANGE range, CUT cut, OPTION opt, string name);
    TF1 *CH1Correction(TTree *t1, EVENT *A, double *t0, charRANGE range, CUT cut, OPTION opt, string name);
    void drawcharacter(TTree *t3, int chN, string* chname, vector<charRANGE> chR);

    bool ifstat(EVENT ch, CUT cut);
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
   
    int Get_npk(){
        return npk;
    }

    gausPAR* Get_PeakPar(){
        return gPar;
    }
    void SetADD(bool val){
        ADD=val;
    }

    private:
    
    bool ADD;
    
    vector<CUT> cut;
    vector<RANGE> U;
    RANGE t;
    string name;
    charRANGE charcut;
    vector<charRANGE> chR;
    char buff[1024];
    int npk;
    gausPAR gPar[100];
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