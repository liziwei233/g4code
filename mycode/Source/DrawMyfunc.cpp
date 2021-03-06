//#include "/mnt/c/Subsys/work/g4code/mycode/Include/DrawMyfunc.h"
//#include "/mnt/c/Subsys/work/g4code/mycode/Include/LZWfunc.h"
//#include "/Users/liziwei/learning/g4code/mycode/Include/LZWfunc.h"
//#include "../Include/LZWfunc.h"
#include "../Include/DrawMyfunc.h"
#include <iostream>
using namespace std;

DrawMyfunc::DrawMyfunc(){
    XTitleOffset=1.0;
    YTitleOffset=1.0;
}

DrawMyfunc::~DrawMyfunc()
{
    cout<<"The destructor is called"<<endl;
}
/*
void DrawMyfunc::Style(){
    // global style 
	   gStyle->SetOptStat(0);
        gStyle->SetOptFit(1111);
	   gStyle->SetOptTitle(0);
        gStyle->SetStatX(0.9);
        gStyle->SetStatY(0.9);
        gStyle->SetStatH(0.15);
        gStyle->SetStatW(0.2);
        gStyle->SetStatStyle(0);
        gStyle->SetPadGridX(1);
        gStyle->SetPadGridY(1);
        //gStyle->SetTitleFont(132,"XYZ");
        //gStyle->SetTitleFont(132,"T");
        gStyle->SetLabelFont(132,"XYZ");
        gStyle->SetLabelSize(0.05,"XYZ");
        gStyle->SetTitleSize(0.05,"XYZ");
	   gStyle->SetTitleSize(0.08,"T");
        //gStyle->SetTitleX(0.55);
        //gStyle->SetFillColor(kGreen);
        //gStyle->SetLabelColor(kRed,"XYZ");
        //gStyle->SetTitleColor(kRed,"XYZ");
        // gStyle->SetTitleColor(kRed,"T");
        //gStyle->SetAxisColor(kRed,"XYZ");
        //gStyle->SetFrameLineColor(kRed);
        gStyle->SetFuncWidth(2);
        gStyle->SetHistLineWidth(2);
        gStyle->SetFrameFillStyle(3003);
        //gStyle->SetFrameFillColor(kGreen);
        //gStyle->SetFuncColor(kBlue);
	   //gStyle->SetTitleSize(0.07,"xyz");
        gStyle->SetTitleOffset(0.95,"xyz");
	   gStyle->SetLabelOffset(0.01,"xyz");
}*/
void DrawMyfunc::Graph(TGraph *datagraph, const char *xtitle, const char *ytitle, Float_t MSize, Int_t MStyle, Color_t MColor, Float_t LWidth, Int_t LStyle, Color_t LColor, Color_t FColor){
     datagraph->SetLineColor( LColor );
     datagraph->SetLineWidth( LWidth );
     datagraph->SetLineStyle( LStyle );
     datagraph->SetMarkerSize( MSize );
     datagraph->SetMarkerStyle( MStyle );
     datagraph->SetMarkerColor( MColor );
     datagraph->SetFillColor( FColor );
     //datagraph->SetFillStyle( FStyle );
     datagraph->GetXaxis()->SetTitle( xtitle);
     datagraph->GetYaxis()->SetTitle( ytitle);
     datagraph->GetXaxis()->SetAxisColor(1);
     datagraph->GetYaxis()->SetAxisColor(1);
     datagraph->GetXaxis()->SetLabelColor(1);
     datagraph->GetYaxis()->SetLabelColor(LColor);
     datagraph->GetXaxis()->SetLabelFont( 42 );
     datagraph->GetYaxis()->SetLabelFont( 42 );
     datagraph->GetXaxis()->SetLabelSize( 0.05 );
     datagraph->GetYaxis()->SetLabelSize( 0.05 );
     datagraph->GetXaxis()->SetLabelOffset( 0.01 );
     datagraph->GetYaxis()->SetLabelOffset( 0.01 );
     datagraph->GetXaxis()->SetTitleFont( 42 );
     datagraph->GetYaxis()->SetTitleFont( 42 );
     //datagraph->GetXaxis()->SetTitleColor( TitleColor);
     //datagraph->GetYaxis()->SetTitleColor( TitleColor );
     datagraph->GetXaxis()->SetTitleSize(0.06);
     datagraph->GetYaxis()->SetTitleSize(0.06);
     datagraph->GetXaxis()->SetTitleOffset(XTitleOffset);
     datagraph->GetYaxis()->SetTitleOffset(YTitleOffset);
     datagraph->GetXaxis()->SetNdivisions(505);     
     datagraph->GetYaxis()->SetNdivisions(510);
}

void DrawMyfunc::Hist(TH1 *datahist, const char *xtitle,const char *ytitle, Float_t LWidth, Int_t LStyle, Color_t LColor, Color_t TitleColor){
     
     datahist->SetLineColor( LColor );
     datahist->SetLineWidth( LWidth );
     datahist->SetLineStyle( LStyle );
     datahist->GetXaxis()->SetTitle( xtitle);
     datahist->GetYaxis()->SetTitle( ytitle);
     datahist->GetXaxis()->SetAxisColor(1);
     datahist->GetYaxis()->SetAxisColor(1);
     datahist->GetXaxis()->SetLabelColor(1);
     datahist->GetYaxis()->SetLabelColor(1);
     datahist->GetXaxis()->SetLabelFont( 42 );
     datahist->GetYaxis()->SetLabelFont( 42 );
     datahist->GetXaxis()->SetLabelSize( 0.06 );
     datahist->GetYaxis()->SetLabelSize( 0.06 );
     datahist->GetXaxis()->SetLabelOffset( 0.01 );
     datahist->GetYaxis()->SetLabelOffset( 0.01 );
     datahist->GetXaxis()->SetTitleFont( 42 );
     datahist->GetYaxis()->SetTitleFont( 42 );
     datahist->GetXaxis()->SetTitleColor( TitleColor);
     datahist->GetYaxis()->SetTitleColor( TitleColor );
     datahist->GetXaxis()->SetTitleSize(0.06);
     datahist->GetYaxis()->SetTitleSize(0.06);
     datahist->GetXaxis()->SetTitleOffset(XTitleOffset);
     datahist->GetYaxis()->SetTitleOffset(YTitleOffset);
     //datahist->GetXaxis()->SetBorderSize(5);
     datahist->GetXaxis()->SetNdivisions(505);     
     datahist->GetYaxis()->SetNdivisions(510);
     datahist->GetXaxis()->CenterTitle();
     datahist->GetYaxis()->CenterTitle();
     
     
}
void DrawMyfunc::Hist(TH2 *datahist, const char *xtitle,const char *ytitle, Float_t LWidth, Int_t LStyle, Color_t LColor, Color_t TitleColor){
     datahist->SetLineColor( LColor );
     datahist->SetLineWidth( LWidth );
     datahist->GetXaxis()->SetTitle( xtitle);
     datahist->GetYaxis()->SetTitle( ytitle);
     datahist->GetXaxis()->SetAxisColor(1);
     datahist->GetYaxis()->SetAxisColor(1);
     datahist->GetXaxis()->SetLabelColor(1);
     datahist->GetYaxis()->SetLabelColor(1);
     datahist->GetXaxis()->SetLabelFont( 42 );
     datahist->GetYaxis()->SetLabelFont( 42 );
     datahist->GetXaxis()->SetLabelSize( 0.06 );
     datahist->GetYaxis()->SetLabelSize( 0.06 );
     datahist->GetXaxis()->SetLabelOffset( 0.01 );
     datahist->GetYaxis()->SetLabelOffset( 0.01 );
     datahist->GetXaxis()->SetTitleFont( 42 );
     datahist->GetYaxis()->SetTitleFont( 42 );
     datahist->GetXaxis()->SetTitleColor( TitleColor);
     datahist->GetYaxis()->SetTitleColor( TitleColor );
     datahist->GetXaxis()->SetTitleSize(0.06);
     datahist->GetYaxis()->SetTitleSize(0.06);
     datahist->GetXaxis()->SetTitleOffset(XTitleOffset);
     datahist->GetYaxis()->SetTitleOffset(YTitleOffset);
     //datahist->GetXaxis()->SetBorderSize(5);
     datahist->GetXaxis()->SetNdivisions(510);     
     datahist->GetYaxis()->SetNdivisions(510);
     datahist->GetXaxis()->CenterTitle();
     datahist->GetYaxis()->CenterTitle();
}

void DrawMyfunc::Pad(TVirtualPad *pad,const char* xname,const char* yname,float x1,float y1,float x2,float y2){


        TH1F* hpad = pad->DrawFrame(x1,y1,x2,y2);
        hpad->GetXaxis()->SetTitle( xname);
        hpad->GetYaxis()->SetTitle( yname);
        hpad->GetXaxis()->SetAxisColor(1);
        hpad->GetYaxis()->SetAxisColor(1);
        hpad->GetXaxis()->SetLabelColor(1);
        hpad->GetYaxis()->SetLabelColor(1);
        hpad->GetXaxis()->SetLabelFont( 42 );
        hpad->GetYaxis()->SetLabelFont( 42 );
        hpad->GetXaxis()->SetLabelSize( 0.05 );
        hpad->GetYaxis()->SetLabelSize( 0.05 );
        hpad->GetXaxis()->SetLabelOffset( 0.01 );
        hpad->GetYaxis()->SetLabelOffset( 0.01 );
        hpad->GetXaxis()->SetTitleFont( 42 );
        hpad->GetYaxis()->SetTitleFont( 42 );
        //hpad->GetXaxis()->SetTitleColor( TitleColor);
        //hpad->GetYaxis()->SetTitleColor( TitleColor );
        hpad->GetXaxis()->SetTitleSize(0.06);
        hpad->GetYaxis()->SetTitleSize(0.06);
        hpad->GetXaxis()->SetTitleOffset(XTitleOffset);
        hpad->GetYaxis()->SetTitleOffset(YTitleOffset);
        pad->Modified();
        pad->Update();

}

void DrawMyfunc::SetPad(TVirtualPad *pad,float left, float bottom, float right, float top){
  pad->SetFillColor(10);
  pad->SetBorderMode(0);
  pad->SetBorderSize(0);
  pad->SetFrameFillColor(10);
  //pad->SetFrameFillStyle(3003);
  pad->SetFrameBorderMode(0);
  pad->SetFrameBorderSize(0);
  pad->SetLeftMargin(left);
  pad->SetRightMargin(right);
  pad->SetTopMargin(top);
  pad->SetBottomMargin(bottom);
}

TLegend* DrawMyfunc::Leg(Double_t xlow, Double_t ylow, Double_t xup, Double_t yup, Int_t textFont, Size_t textSize){
  TLegend *leg = new TLegend(xlow,ylow,xup,yup);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetFillColor(10);
  leg->SetTextFont(textFont);
  leg->SetTextSize(textSize);
  //leg->Draw("same");
  return leg;
}
TLatex* DrawMyfunc::Latex( Double_t x, Double_t y,const char* text, Int_t textFont, Size_t textSize, Color_t colorIndex){
  TLatex *latex = new TLatex(x,y,text);
  latex->SetNDC();
  latex->SetTextFont(textFont);
  latex->SetTextSize(textSize);
  latex->SetTextColor(colorIndex);
  //latex->Draw("same");
  return latex;
}
void DrawMyfunc::Setstat(TH1* datahist, double x1,double y1,double x2,double y2){
     
     TPaveStats *st = (TPaveStats*)datahist->FindObject("stats");
     st->SetX1NDC(x1);
     st->SetX2NDC(x2);
     st->SetY1NDC(y1);
     st->SetY2NDC(y2);
}

TPad* DrawMyfunc::BuildPad(const int Npad, double* ratio){
  //* example : build 3 pad ,and the ratio is 1:3:5
  //* so set Npad=3, ratio[3]={1,3,5}
  char buff[1024];
  TPad* pad[Npad];
  double y[Npad];
  double hp=0.98;
  for(int i=0;i<Npad;i++)
  {
    y[i]=hp-(0.98-0.01)*ratio[i]/TMath::Mean(Npad,ratio)*Npad;
    sprintf(buff,"pad%d",i);
    pad[i]=new TPad(buff,"",0.01,y[i],0.98,hp);
    hp=y[i];
  }
  return pad[Npad];
}

void DrawMyfunc::drawline(float x1, float x2){

	TLine* line1 = new TLine(x1,gPad->VtoPixel(gPad->GetUymin()),x1,gPad->VtoPixel(gPad->GetUymax()));
    TLine* line2 = new TLine(x2,gPad->VtoPixel(gPad->GetUymin()),x2,gPad->VtoPixel(gPad->GetUymax()));
    
    line1->SetLineStyle(2);
    line1->SetLineColor(6);
    line1->Draw("same");
    line2->SetLineStyle(2);
    line2->SetLineColor(6);
    line2->Draw("same");
}