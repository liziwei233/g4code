#ifndef __drawmyfunc_h__
#define __drawmyfunc_h__

#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TVirtualPad.h"

class DrawMyfunc
{
    public:
    DrawMyfunc();
    ~DrawMyfunc();
    
    void Graph(TGraph *datagraph, const char *xtitle, const char *ytitle, Float_t MSize=1, Int_t MStyle =28, Color_t MColor=1, Float_t LWidth=1, Int_t LStyle=1, Color_t LColor=1, Color_t FColor=16);
    void Hist(TH1 *datahist, const char *xtitle,const char *ytitle, Float_t LWidth=1, Int_t LStyle=1, Color_t LColor=1, Color_t TitleColor=1);
    void Hist(TH2 *datahist, const char *xtitle,const char *ytitle, Float_t LWidth=1, Int_t LStyle=1, Color_t LColor=1, Color_t TitleColor=1);
    void Pad(TVirtualPad *pad,const char* xname,const char* yname,float x1,float y1,float x2,float y2);
    TLegend* Leg(Double_t xlow=0.2, Double_t ylow=0.2, Double_t xup=0.5, Double_t yup=0.5, Int_t textFont=62, Size_t textSize=0.05);
    TLatex* Latex( Double_t x=0.65, Double_t y=0.5,const char* text="", Int_t textFont=62, Size_t textSize=0.05, Color_t colorIndex=2);
    void SetPad(TVirtualPad *pad,float left, float bottom, float right, float top);
    void SetXTitleOffset(double x){ XTitleOffset=x;};
    void SetYTitleOffset(double y){ XTitleOffset=y;};
    void Setstat(TH1* datahist, double x1,double y1,double x2,double y2);
    TPad* BuildPad(const int Npad, double* ratio);
    private:
    double XTitleOffset;
    double YTitleOffset;
    
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