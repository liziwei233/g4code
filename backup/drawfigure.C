#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TAttAxis.h>
#include <TStyle.h>
#include <TRandom.h>

void drawfigure(){

 
void DrawMyGraph(TGraph *datagraph, const char *xtitle, const char *ytitle, Float_t MSize=1, Int_t MStyle =28, Color_t MColor=1, Color_t LColor=1, Float_t LWidth=1, Int_t LStyle=1, Color_t FColor=16);
void DrawMyHist1(TH1 *datahist, const char *xtitle,const char *ytitle, Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);
void DrawMyHist2(TH2 *datahist, const char *xtitle,const char *ytitle, Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);
void SetMyPad(TVirtualPad *pad,float left, float right, float top, float bottom);
TLegend* DrawMyLeg(Double_t xlow=0.2, Double_t ylow=0.2, Double_t xup=0.5, Double_t yup=0.5, Int_t textFont=62, Size_t textSize=0.05);

gStyle->SetOptTitle(0);
gStyle->SetOptStat(1111);

float x1[]={-75,-30,0,5,10,15,20,25,30,75}; //deep=width changed
float y1[]={1.7,2,4,5,18,152,140,139,145,151};
float y2[]={2400,2322,2011,2019,2043,2083,2139,2218,2322,2399};
float y3[]={452,238,148,119,21,6,7,7,7,3};

//float y1[]={62.79,33.52,22.88,22.27,21.18};
//float x2[]={5,10,20,30,40}; //d=10mm,width changed
//float y2[]={33.42,33.52,32.76,33.29,32.46};
//float x3[]={30,45,60,90};
//float y3[]={,24.26,,33.52};
//float ys1[]={}
int n=sizeof(x1)/sizeof((x1)[0]);

TGraph *gTR = new TGraph(n,x1,y1); 
TGraph *gLY = new TGraph(n,x1,y2); 
TGraph *gHit = new TGraph(n,x1,y3);

DrawMyGraph(gTR,"theta (deg)","Time resolution (ps)",1.5,20,4,4,2,7);
DrawMyGraph(gLY,"theta (deg)","Time resolution (ps)",2.,22,2,2,2,7);
DrawMyGraph(gHit,"theta (deg)","Time resolution (ps)",2.,22,2,2,2,7);
/*
DrawMyGraph(gR1,"Rate","Time resolution (ps)",1.5,20,1,1,2,1);
DrawMyGraph(gR2,"Rate","Time resolution (ps)",1,21,2,2,2,1);
DrawMyGraph(gth,"threshold (mV)","Time resolution (ps)",1.5,20,1,1,2,1);
DrawMyGraph(gth_TA,"threshold (mV)","Time resolution (ps)",1,21,2,2,2,1);
DrawMyGraph(gAn,"Angle (deg) ","Time resolution (ps)",1.5,20,1,1,2,1);
//DrawMyGraph(gH,"RF2","",1,4,8,8,3,1);
*/
TCanvas *c1 = new TCanvas("c1","c1",800,600);
c1->cd(); 
c1->Divide(1,3);
//gR1->Draw("ALP");
c1->cd(3);
SetMyPad(gpad,0.01,0.01,0.01,0.01);
gTR->Draw("ALP");
//gh->GetYaxis()->SetRangeUser(10,75);
c1->cd(1);
SetMyPad(gpad,0.01,0.01,0.01,0.01);
gLY->Draw("ALP");
c1->cd(2);
SetMyPad(gpad,0.01,0.01,0.01,0.1);
gHit->Draw("ALP");

/*
gth->Draw("ALP");
gth->GetYaxis()->SetRangeUser(5,25);
gth_TA->Draw("LPsame");
*/



/*	
TGaxis *axis = new TGaxis(7.6,1,7.6,6,1400,1550,510,"L");
axis->SetLabelColor(8);
axis->SetLabelFont(42);
axis->SetLabelOffset(-0.02);
axis->SetLabelSize(0.05);
axis->Draw();


TAxis* a = gR->GetXaxis();
a->SetNdivisions(110);
a->ChangeLabel(1,-1,-1,-1,-1,-1,"330K");
a->ChangeLabel(2,-1,-1,-1,-1,-1,"680K");
a->ChangeLabel(3,-1,-1,-1,-1,-1,"750K");
a->ChangeLabel(4,-1,-1,-1,-1,-1,"820K");
a->ChangeLabel(5,-1,-1,-1,-1,-1,"1000K");
a->ChangeLabel(6,-1,-1,-1,-1,-1,"1200K");
a->ChangeLabel(7,-1,-1,-1,-1,-1,"1500K");
gR->Draw("same");
*/

		
TLegend* leg;
leg = DrawMyLeg(0.55,0.68,0.65,0.85);
//leg->AddEntry(gR1,"with threshold = 0.2mV","lp");
//leg->AddEntry(gR2,"without ","lp");
leg->AddEntry(gLY,"Lightyieled","lp");
leg->AddEntry(gTR,"Time Res","lp");
leg->AddEntry(gHit,"Hits","lp");
leg->Draw();


/*leg->AddEntry(gF,"FWHM (ns)","lp");
leg->AddEntry(ge,"#mu/#sigma (*10)","lp");
leg->AddEntry(gH,"workHV (V)","lp");
leg->Draw();
*/

c1->SaveAs("TR.png");


}

void DrawMyGraph(TGraph *datagraph, const char *xtitle, const char *ytitle, Float_t MSize=1, Int_t MStyle =28, Color_t MColor=1, Color_t LColor=1, Float_t LWidth=1, Int_t LStyle=1, Color_t FColor=16){
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
     datagraph->GetYaxis()->SetLabelColor(1);
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
     datagraph->GetXaxis()->SetTitleOffset(0.8);
     datagraph->GetYaxis()->SetTitleOffset(0.8);
}

void DrawMyHist1(TH1 *datahist, const char *xtitle,const char *ytitle, Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1){
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
     datahist->GetXaxis()->SetTitleOffset(1.0);
     datahist->GetYaxis()->SetTitleOffset(0.8);
     //datahist->GetXaxis()->SetBorderSize(5);
     datahist->GetXaxis()->SetNdivisions(510);     
     datahist->GetYaxis()->SetNdivisions(510);
     datahist->GetXaxis()->CenterTitle();
     datahist->GetYaxis()->CenterTitle();
}
void DrawMyHist2(TH2 *datahist, const char *xtitle,const char *ytitle, Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1){
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
     datahist->GetXaxis()->SetTitleOffset(1.0);
     datahist->GetYaxis()->SetTitleOffset(0.8);
     //datahist->GetXaxis()->SetBorderSize(5);
     datahist->GetXaxis()->SetNdivisions(510);     
     datahist->GetYaxis()->SetNdivisions(510);
     datahist->GetXaxis()->CenterTitle();
     datahist->GetYaxis()->CenterTitle();
}
void SetMyPad(TVirtualPad *pad,float left, float right, float top, float bottom){
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
TLegend* DrawMyLeg(Double_t xlow=0.2, Double_t ylow=0.2, Double_t xup=0.5, Double_t yup=0.5, Int_t textFont=62, Size_t textSize=0.05){
  TLegend *leg = new TLegend(xlow,ylow,xup,yup);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetFillColor(10);
  leg->SetTextFont(textFont);
  leg->SetTextSize(textSize);
  //leg->Draw("same");
  return leg;
}
