#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TAttAxis.h>
#include <TStyle.h>
#include <TRandom.h>

void output(){
    gStyle->SetOptStat(1111);
 
    void DrawMyGraph(TGraph *datagraph, const char *xtitle, const char *ytitle, Float_t MSize=1, Int_t MStyle =28, Int_t MColor=1, Int_t LColor=1, Float_t LWidth=1, Int_t LStyle=1, Int_t FColor=16);
    void DrawMyHist1(TH1 *datahist, const char *xtitle,const char *ytitle, Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);
    void DrawMyHist2(TH2 *datahist, const char *xtitle,const char *ytitle, Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);
    void SetMyPad(TPad *pad,float left, float right, float top, float bottom);
    TLegend* DrawMyLeg(Double_t xlow=0.2, Double_t ylow=0.2, Double_t xup=0.5, Double_t yup=0.5, Int_t textFont=62, Size_t textSize=0.05);
    Int_t ID = 0,B = 0,binN = 0;
    Double_t L = 0;// path length (cm)
    Double_t t = 0; // hit time (ns)
    Double_t wvL=0; // waveLength (nm)
    bool flag=0;
    Int_t temp=0;
    Int_t detected = 0; //val=0,undetected; val=1,detected by L; val=2,detected by R;

    Int_t counter1 = 0,counter2 = 0;
    Double_t x1[300]={0},x2[300]={0};
    Double_t yB1[300]={0},yB2[300]={0};
    Double_t yL1[300]={0},yL2[300]={0};
    Double_t yW1[300]={0},yW2[300]={0};
    Double_t yT1[300]={0},yT2[300]={0};

    TFile *f1 = new TFile("hist.root","RECREATE");
    TTree *t1 = new TTree("t1","absorbed photons");
    TTree *t2 = new TTree("t2","hit photons");
    t1->Branch("ID",&ID,"ID/I");
    t1->Branch("B",&B,"B/I");
    t1->Branch("L",&L,"L/D");
    t1->Branch("t",&t,"t/D");
    t1->Branch("wvL",&wvL,"wvL/D");


    t2->Branch("ID",&ID,"ID/I");
    t2->Branch("B",&B,"B/I");
    t2->Branch("L",&L,"L/D");
    t2->Branch("t",&t,"t/D");
    t2->Branch("wvL",&wvL,"wvL/D");

    TH1D *HB1 = new TH1D("B_abs","",1e3,0,2e4);
    TH1D *HB2 = new TH1D("B_hit","",1e3,0,2e4);
    TH1D *HL1 = new TH1D("L_abs","",1e3,0,2e4);
    TH1D *HL2 = new TH1D("L_hit","",1e3,0,2e4);
    TH1D *HW1 = new TH1D("wvl_abs","",200,100,9e3);
    TH1D *HW2 = new TH1D("wvl_hit","",200,100,9e3);
    TH1D *HT1 = new TH1D("T_abs","",4e3,0,2e3);
    TH1D *HT2 = new TH1D("T_hit","",4e3,0,2e3);

    TH2D *H1 = new TH2D("Absorbed","",2e4,0,2e4,800,100,900);
    TH2D *H2 = new TH2D("hitpmt","",2e4,0,2e4,800,100,900);

    DrawMyHist1(HT1,"Time (ns)","Counts",1,2);
    DrawMyHist1(HT2,"Time (ns)","Counts",2,2);

    DrawMyHist1(HL1,"path length (cm)","Counts",1,2);
    DrawMyHist1(HL2,"path length (cm)","Counts",2,2);

    DrawMyHist1(HW1,"wavelength (nm)","Counts",1,2);
    DrawMyHist1(HW2,"wavelength (nm)","Counts",2,2);

    DrawMyHist1(HB1,"Bounce","Counts",1,2);
    DrawMyHist1(HB2,"Bounce","Counts",2,2);

    DrawMyHist2(H1,"time (ns)","wavelength (nm)",1,2);
    DrawMyHist2(H2,"time (ns)","wavelength (nm)",2,2);
    H1->SetMarkerSize( 1 );
    H1->SetMarkerStyle( 28 );
    H1->SetMarkerColor( 1 );

    H2->SetMarkerSize( 1 );
    H2->SetMarkerStyle( 26 );
    H2->SetMarkerColor( 2 );

    ifstream indata;
    indata.open("PhotonsInform.dat");
    while(indata&&(!indata.eof())){
       
		indata>>ID>>wvL>>B>>L>>t>>detected;

		
         // undetected 
        if (detected==0) {
        x1[counter1]=ID;
        yB1[counter1]=B;
        yL1[counter1]=L;
        yW1[counter1]=wvL;
        yT1[counter1]=t;
        HB1->Fill(B);
        HL1->Fill(L);
        HW1->Fill(wvL);
        HT1->Fill(t);
        H1->Fill(t,wvL);
        t1->Fill();
        //cout<<x1[counter1]<<endl;
        //cout<<yB1[counter1]<<endl;
        counter1++;
        }
        else {
        x2[counter2]=ID;
        yB2[counter2]=B;
        yL2[counter2]=L;
        yW2[counter2]=wvL;
        yT2[counter2]=t;

        HB2->Fill(B);
        HL2->Fill(L);
        HW2->Fill(wvL);
        HT2->Fill(t);
        H2->Fill(t,wvL);
        t2->Fill();
        counter2++;

         }
    }

    t1->Write();
    t2->Write();

    H1->Write();
    H2->Write();

    cout<<" detected = "<<counter2<<endl;
    cout<<" undetected = "<<counter1<<endl;

    indata.close();
    
    /*
    TGraph *gB1 = new TGraph(100,x1,yB1);
    
    TGraph *gL1 = new TGraph(counter1,x1,yL1);
    TGraph *gW1 = new TGraph(counter1,x1,yW1);
    TGraph *gT1 = new TGraph(counter1,x1,yT1);

    TGraph *gB2 = new TGraph(counter2,x2,yB2);
    TGraph *gL2 = new TGraph(counter2,x2,yL2);
    TGraph *gW2 = new TGraph(counter2,x2,yW2);
    TGraph *gT2 = new TGraph(counter2,x2,yT2);

    DrawMyGraph(gT1,"ID", "Time", 1,28, 1,1, 1,1);
    DrawMyGraph(gT2,"ID", "Time", 1,26, 2,2, 1,1);

    TCanvas *c = new TCanvas("c","",1000,600);
    c->cd();
    gPad->SetLogy();
    gT1->Draw("AP");
    gT1->GetYaxis()->SetRangeUser(0,100);
    
    gB1->
    gT2->Draw("Psame");
    */

    TCanvas *c2 = new TCanvas("c2","",1000,600);
    SetMyPad(c2,0.1, 0.1,0.1, 0.15);
    c2->cd();
    c2->SetLogx();
    H2->Draw();
    H1->Draw("SAMES");
    TGaxis::SetMaxDigits(3);
    TLegend* leg; 
	leg = DrawMyLeg(0.6,0.2,0.75,0.32);
	leg->AddEntry(H2,"Hit PMT","lp");
	leg->AddEntry(H1,"Absorbed","lp");
	leg->Draw();

    gPad->Update(); // inevitable!!
    TPaveStats *stL = (TPaveStats*)H2->FindObject("stats");
	//stR->SetName("PMT_Right");
	stL->SetY1NDC(0.4);
	stL->SetY2NDC(0.68);
	//gPad->Modified();
	//gPad->Update();
    
}
void DrawMyGraph(TGraph *datagraph, const char *xtitle, const char *ytitle, Float_t MSize=1, Int_t MStyle =28, Int_t MColor=1, Int_t LColor=1, Float_t LWidth=1, Int_t LStyle=1, Int_t FColor=16){
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
void SetMyPad(TPad *pad,float left, float right, float top, float bottom){
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