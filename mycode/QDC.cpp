#include <string>
#include <time.h>
#include <TString.h>
#include <TStyle.h>
using namespace std;
#define verbose
char buff[1024];
char line[180];
int counter = 0;
char path[1024] = "/mnt/c/Users/liziwei/OneDrive/Documents/paper/picture/Picture_program/saveroot";
TFile *sfile = NULL;
TFile *prefile = NULL;

void setgStyle()
{

    gStyle->SetFrameLineWidth(3);
    //gStyle->SetFrameBorderSize(2);
    gStyle->SetTickLength(0.04);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(1111);
    gStyle->SetErrorX(1);
    gStyle->SetEndErrorSize(4);
}

void DrawMyGraph(TGraph *datagraph, char *xtitle, char *ytitle, Size_t MSize = 1, Style_t MStyle = 28, Color_t MColor = 1, Color_t LColor = 1, Width_t LWidth = 1, Style_t LStyle = 1, Color_t FColor = 0)
{
    datagraph->SetLineColor(LColor);
    datagraph->SetLineWidth(LWidth);
    datagraph->SetLineStyle(LStyle);
    datagraph->SetMarkerSize(MSize);
    datagraph->SetMarkerStyle(MStyle);
    datagraph->SetMarkerColor(MColor);
    datagraph->SetFillColor(FColor);
    //datagraph->SetFillStyle( FStyle );
    datagraph->GetXaxis()->SetTitle(xtitle);
    datagraph->GetYaxis()->SetTitle(ytitle);
    datagraph->GetXaxis()->SetAxisColor(1);
    datagraph->GetYaxis()->SetAxisColor(1);
    datagraph->GetXaxis()->SetLabelColor(1);
    datagraph->GetYaxis()->SetLabelColor(1);
    datagraph->GetXaxis()->SetLabelFont(42);
    datagraph->GetYaxis()->SetLabelFont(42);
    datagraph->GetXaxis()->SetLabelSize(0.07);
    datagraph->GetYaxis()->SetLabelSize(0.07);
    datagraph->GetXaxis()->SetLabelOffset(0.01);
    datagraph->GetYaxis()->SetLabelOffset(0.01);
    datagraph->GetXaxis()->SetTitleFont(42);
    datagraph->GetYaxis()->SetTitleFont(42);
    //datagraph->GetXaxis()->SetTitleColor( TitleColor);
    //datagraph->GetYaxis()->SetTitleColor( TitleColor );
    datagraph->GetXaxis()->SetTitleSize(0.07);
    datagraph->GetYaxis()->SetTitleSize(0.07);
    datagraph->GetXaxis()->SetTitleOffset(1.0);
    datagraph->GetYaxis()->SetTitleOffset(1.0);
    datagraph->GetXaxis()->CenterTitle(1);
    datagraph->GetYaxis()->CenterTitle(1);
}
void DrawMyHist(TH1 *datahist, char *xtitle, char *ytitle, Color_t LColor = 1, Width_t LWidth = 3, Color_t TitleColor = 1)
{
    datahist->SetLineColor(LColor);
    datahist->SetLineWidth(LWidth);

    datahist->GetXaxis()->SetTitle(xtitle);
    datahist->GetYaxis()->SetTitle(ytitle);
    datahist->GetXaxis()->SetAxisColor(1);
    datahist->GetYaxis()->SetAxisColor(1);
    datahist->GetXaxis()->SetLabelColor(1);
    datahist->GetYaxis()->SetLabelColor(1);
    datahist->GetXaxis()->SetLabelFont(42);
    datahist->GetYaxis()->SetLabelFont(42);
    datahist->GetXaxis()->SetLabelSize(0.07);
    datahist->GetYaxis()->SetLabelSize(0.07);
    datahist->GetXaxis()->SetLabelOffset(0.01);
    datahist->GetYaxis()->SetLabelOffset(0.01);
    datahist->GetXaxis()->SetTitleFont(42);
    datahist->GetYaxis()->SetTitleFont(42);
    datahist->GetXaxis()->SetTitleColor(TitleColor);
    datahist->GetYaxis()->SetTitleColor(TitleColor);
    datahist->GetXaxis()->SetTitleSize(0.07);
    datahist->GetYaxis()->SetTitleSize(0.07);
    datahist->GetXaxis()->SetTitleOffset(1.0);
    datahist->GetYaxis()->SetTitleOffset(1.0);
    //datahist->GetXaxis()->SetBorderSize(5);

    datahist->GetXaxis()->SetNdivisions(510);
    datahist->GetYaxis()->SetNdivisions(510);
    datahist->GetXaxis()->CenterTitle(1);
    datahist->GetYaxis()->CenterTitle(1);

     
}

void DrawMy2dHist(TH2 *datahist, char *xtitle, char *ytitle, Color_t LColor = 1, Width_t LWidth = 3, Color_t TitleColor = 1)
{
    datahist->SetLineColor(LColor);
    datahist->SetLineWidth(LWidth);
      datahist->GetXaxis()->SetTitle(xtitle);
    datahist->GetYaxis()->SetTitle(ytitle);
    datahist->GetXaxis()->SetAxisColor(1);
    datahist->GetYaxis()->SetAxisColor(1);
    datahist->GetXaxis()->SetLabelColor(1);
    datahist->GetYaxis()->SetLabelColor(1);
    datahist->GetXaxis()->SetLabelFont(42);
    datahist->GetYaxis()->SetLabelFont(42);
    datahist->GetXaxis()->SetLabelSize(0.07);
    datahist->GetYaxis()->SetLabelSize(0.07);
    datahist->GetXaxis()->SetLabelOffset(0.01);
    datahist->GetYaxis()->SetLabelOffset(0.01);
    datahist->GetXaxis()->SetTitleFont(42);
    datahist->GetYaxis()->SetTitleFont(42);
    datahist->GetXaxis()->SetTitleColor(TitleColor);
    datahist->GetYaxis()->SetTitleColor(TitleColor);
    datahist->GetXaxis()->SetTitleSize(0.07);
    datahist->GetYaxis()->SetTitleSize(0.07);
    datahist->GetXaxis()->SetTitleOffset(1.0);
    datahist->GetYaxis()->SetTitleOffset(1.0);
    //datahist->GetXaxis()->SetBorderSize(5);
    datahist->GetXaxis()->SetNdivisions(510);
    datahist->GetYaxis()->SetNdivisions(510);
    datahist->GetXaxis()->CenterTitle(1);
    datahist->GetYaxis()->CenterTitle(1);
     
}
void DrawMyfun(TF1 *datafunc, char *xtitle, char *ytitle, Color_t LColor = 1, Width_t LWidth = 3, Style_t LStyle = 1)
{
    datafunc->SetLineColor(LColor);
    datafunc->SetLineWidth(LWidth);
    datafunc->SetLineStyle(LStyle);
    datafunc->GetXaxis()->SetTitle(xtitle);
    datafunc->GetYaxis()->SetTitle(ytitle);
    datafunc->GetXaxis()->SetAxisColor(1);
    datafunc->GetYaxis()->SetAxisColor(1);
    datafunc->GetXaxis()->SetLabelColor(1);
    datafunc->GetYaxis()->SetLabelColor(1);
    datafunc->GetXaxis()->SetLabelFont(42);
    datafunc->GetYaxis()->SetLabelFont(42);
    datafunc->GetXaxis()->SetLabelSize(0.07);
    datafunc->GetYaxis()->SetLabelSize(0.07);
    datafunc->GetXaxis()->SetLabelOffset(0.01);
    datafunc->GetYaxis()->SetLabelOffset(0.01);
    datafunc->GetXaxis()->SetTitleFont(42);
    datafunc->GetYaxis()->SetTitleFont(42);
    datafunc->GetXaxis()->SetTitleColor(1);
    datafunc->GetYaxis()->SetTitleColor(1);
    datafunc->GetXaxis()->SetTitleSize(0.07);
    datafunc->GetYaxis()->SetTitleSize(0.07);
    datafunc->GetXaxis()->SetTitleOffset(1.0);
    datafunc->GetYaxis()->SetTitleOffset(1.0);
    //datahist->GetXaxis()->SetBorderSize(5);
    datafunc->GetXaxis()->SetNdivisions(510);
    datafunc->GetYaxis()->SetNdivisions(510);
    datafunc->GetXaxis()->CenterTitle(1);
    datafunc->GetYaxis()->CenterTitle(1);
     
}
void Drawline(float x1, float x2, Width_t LWidth = 1, Style_t LStyle = 1, Color_t LColor = 6)
{

    //TLine* line1 = new TLine(x1,gPad->VtoPixel(gPad->GetUymin()),x1,gPad->VtoPixel(gPad->GetUymax()));
    //TLine* line2 = new TLine(x2,gPad->VtoPixel(gPad->GetUymin()),x2,gPad->VtoPixel(gPad->GetUymax()));
    gPad->Update();
    gPad->Modified();
    TLine *line1 = new TLine(x1, gPad->GetUymin(), x1, gPad->GetUymax());
    TLine *line2 = new TLine(x2, gPad->GetUymin(), x2, gPad->GetUymax());

    line1->SetLineWidth(LWidth);
    line1->SetLineStyle(LStyle);
    line1->SetLineColor(LColor);
    line1->Draw("same");
    line2->SetLineWidth(LWidth);
    line2->SetLineStyle(LStyle);
    line2->SetLineColor(LColor);
    line2->Draw("same");
}

TLegend *DrawMyLeg(Double_t xlow = 0.3, Double_t ylow = 0.6, Double_t xup = 0.6, Double_t yup = 0.9, Int_t textFont = 62, Double_t textSize = 0.03)
{
    TLegend *leg = new TLegend(xlow, ylow, xup, yup);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetFillColor(10);
    leg->SetTextFont(textFont);
    leg->SetTextSize(textSize);
    //leg->Draw("same");
    return leg;
}
void SetMyPad(TVirtualPad *pad, float left, float right, float top, float bottom)
{
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
TLatex *DrawMyLatex(char *text, Double_t x = 0.65, Double_t y = 0.5, Int_t textFont = 42, Size_t textSize = 0.06, Color_t colorIndex = 1)
{
    TLatex *latex = new TLatex(x, y, text);
    latex->SetNDC();
    latex->SetTextFont(textFont);
    latex->SetTextSize(textSize);
    latex->SetTextColor(colorIndex);
    latex->Draw("same");
    return latex;
}

TCanvas *cdC(int n)
{
    sprintf(buff, "c%d", n);
    TCanvas *c = new TCanvas(buff, buff, 800, 600);
    c->cd();
    gPad->SetGrid();
    SetMyPad(gPad, 0.15, 0.05, 0.1, 0.14);
    return c;
}

TF1 *gausfit(TH1 *h, double facleft, double facright, int rbU, double UL, double UR)
{
    double mean = 0;
    double sigma = 0;
    TH1 *hU = (TH1 *)h->Clone();
    hU->Draw();
    hU->Rebin(rbU);
    hU->GetXaxis()->SetRangeUser(UL, UR);
    TF1 *fitU = new TF1("fitU", "gaus", UL, UR);
    mean = hU->GetBinCenter(hU->GetMaximumBin());
    fitU->SetParameter(1, mean);
    cout << mean << "\t" << sigma << endl;
    hU->Fit(fitU, "R");
    mean = fitU->GetParameter(1);
    sigma = fitU->GetParameter(2);

    cout << mean << "\t" << sigma << endl;

    TFitResultPtr failed = hU->Fit(fitU, "", "", mean - facleft * sigma, mean + facright * sigma);
    //failed =1 means fit failed
    if (failed)
        return fitU = 0;

    else
    {

        if (UL < mean - 20 * sigma)
            UL = mean - 20 * sigma;
        if (UR > mean + 20 * sigma)
            UR = mean + 20 * sigma;

        hU->GetXaxis()->SetRangeUser(UL, UR);
        sfile->WriteTObject(hU);
        return fitU;
    }
}

TF1 *twogausfit(TH1 *ht, double fac, double rangefac, int rbt, double tL, double tR)
{
    //First fit for ensuring the rangement of histgram;
    TH1 *h = (TH1 *)ht->Clone();
    h->Rebin(rbt);
    double mean = h->GetBinCenter(h->GetMaximumBin());
    double sigma = h->GetRMS();
    TF1 *fit = new TF1("fit", "gaus", tL, tR);
    h->GetXaxis()->SetRangeUser(tL, tR);
    fit->SetParameter(1, mean);
    //fit->SetParameter(2,sigma);
    h->Fit(fit);
    mean = fit->GetParameter(1);
    sigma = TMath::Abs(fit->GetParameter(2));
    if (tL < mean - 5 * sigma || sigma > 1)
    {
        tL = mean - 5 * sigma;
        tR = mean + 10 * sigma;
    }
    cout << h->GetName() << "\t" << tL << "\t" << tR << endl;

    h->GetXaxis()->SetRangeUser(tL, tR);

    TF1 *fit2 = new TF1("fit2", "gaus(0)+gaus(3)", tL, tR);
    fit2->SetParNames("C_{TR}", "#mu_{TR}", "#sigma_{TR}", "C_{bkgnd}", "#mu_{bkgnd}", "#sigma_{bkgnd}");
    fit2->SetParameter(1, mean);
    fit2->SetParameter(2, sigma);
    fit2->SetParLimits(3, 0, fit->GetParameter(0) * fac);
    fit2->SetParameter(4, mean);
    fit2->SetParameter(5, 2 * sigma);

    //h->Fit(fit2);
    h->Fit(fit2, "", "", mean - rangefac * sigma, mean + rangefac * sigma);
    TF1 *fit_tr = new TF1("fit_tr", "gaus", tL, tR);
    fit_tr->SetParameter(0, fit2->GetParameter(0));
    fit_tr->SetParameter(1, fit2->GetParameter(1));
    fit_tr->SetParameter(2, fit2->GetParameter(2));
    TF1 *fit_bg = new TF1("fit_bg", "gaus", tL, tR);
    fit_bg->SetParameter(0, fit2->GetParameter(3));
    fit_bg->SetParameter(1, fit2->GetParameter(4));
    fit_bg->SetParameter(2, fit2->GetParameter(5));
    fit_tr->SetLineColor(3);
    fit_tr->SetLineStyle(7);
    fit_bg->SetLineColor(3);
    fit_bg->SetLineStyle(7);
    fit_tr->Draw("same");
    fit_bg->Draw("same");

    if (fit2->GetParameter(3) >= fit2->GetParameter(0) && TMath::Abs(fit2->GetParameter(5)) >= 0.05)
    {
        mean = fit2->GetParameter(4);
        sigma = TMath::Abs(fit2->GetParameter(5));
    }
    else if (TMath::Abs(fit2->GetParameter(2)) >= 0.05)
    {
        mean = fit2->GetParameter(1);
        sigma = TMath::Abs(fit2->GetParameter(2));
    }

    if (tL < mean - 10 * sigma)
    {
        tL = mean - 10 * sigma;
    }
    if (tR > mean + 10 * sigma)
    {
        tR = mean + 10 * sigma;
    }

    cout << h->GetName() << "\t" << tL << "\t" << tR << endl;
    h->GetXaxis()->SetRangeUser(tL, tR);
    return fit2;
}

TF1 *profilefit(TH2 *Rt, double rbU, double rbt, double tL, double tR, double UL, double UR, char *name)
{

    TCanvas *cAT = new TCanvas("cAT", "cAT", 800, 600);
    TCanvas *cpfx = new TCanvas("cpfx", "cpfx", 800, 600);
    cAT->Clear();
    cpfx->Clear();

    cAT->cd();
    TH2 *Qt = (TH2 *)Rt->Clone();
    //TH2* Qt = (TH2*) Rt->Clone("tmp");
    Qt->Draw("colz");

    Qt->RebinX(rbU * 2);
    Qt->RebinY(rbt * 2);
    Qt->GetYaxis()->SetRangeUser(tL, tR);
    Qt->GetXaxis()->SetRangeUser(UL, UR);
    Qt->ProfileX();

    cpfx->cd();
    TH1 *Qpfx = Qt->ProfileX();
    //Qpfx->Reset();

    //Qpfx=(TH1*)gDirectory->Get("Qt_pfx");
    Qpfx->Draw();
    Qpfx->GetYaxis()->SetRangeUser(tL, tR);
    Qpfx->GetYaxis()->SetTitle("Timediff (ps)");
    Qpfx->GetXaxis()->SetRangeUser(UL, UR);
    Qpfx->GetXaxis()->SetTitle("TOT (ps)");

    //TF1 *fitQt = new TF1("fitQt", "[0]+[1]/TMath::Sqrt(abs(x))+[2]/abs(x)+[3]/abs(x)/TMath::Sqrt(abs(x))+[4]/abs(x)/abs(x)", UL, UR);
    TF1 *fitQt = new TF1("fitQt", "pol5", UL, UR);

    fitQt->SetNpx(1000000);
    Qpfx->Fit(fitQt, "R");

    sprintf(buff, "%s_pfx.png", name);
    cpfx->SaveAs(buff);
    sprintf(buff, "%s_ATrelationship.png", name);
    cAT->SaveAs(buff);
    sfile->WriteTObject(Qpfx);
    sfile->WriteTObject(Qt);

    Qpfx->Reset();
    //delete Qt;
    //delete Qpfx;
    //delete c5;
    return fitQt;
}
/*
void read792N(const char *name)
{
    double qdc;

    sprintf(buff, "%s.root", name);
    TFile *file = new TFile(buff, "recreate");
    TTree *tree = new TTree("tree", "");
    tree->Branch("qdc", &qdc, "qdc/D");
    sprintf(buff, "%s.txt", name);
    ifstream input;
    input.open(buff);
    if (!input)
    {
        cout << "file can't be found!" << endl;
        return;
    }
    while (input && (!input.eof()))
    {
        input >> qdc;
        tree->Fill();
    }
    input.close();
    if (file->IsOpen())
    {
        tree->Write();
    }
}
void readV965(const char *name)
{
    double qdc[32];

    sprintf(buff, "%s.root", name);
    TFile *file = new TFile(buff, "recreate");
    TTree *tree = new TTree("tree", "");
    tree->Branch("qdc", qdc, "qdc[32]/D");
    sprintf(buff, "%s.txt", name);
    ifstream input;
    input.open(buff);
    if (!input)
    {
        cout << "file can't be found!" << endl;
        return;
    }
   

    while (input && (!input.eof()))
    {
        for (int i = 0; i < 32; i++)
            input >> qdc[i];
        tree->Fill();
    }
    input.close();
    if (file->IsOpen())
    {
        tree->Write();
    }
}

void get792NHist(const char *name)
{


    double qdc;
    sprintf(buff, "%s.root", name);
    if (gSystem->AccessPathName(buff))
    {
        read792N(name);
    }
    TFile *file = new TFile(buff, "read");
    TTree *tree = (TTree *)file->Get("tree");
    tree->SetBranchAddress("qdc", &qdc);

    TH1D *hq = new TH1D("hq", "", 400, 0, 400);

    int N = tree->GetEntries();
    for (int i = 0; i < N; i++)
    {
        tree->GetEntry(i);
        hq->Fill(qdc);
    }
    DrawMyHist(hq, "qdc (100fC)", "Counts", 1, 2);
    TCanvas *c1 = cdC(1);
    hq->Draw();
}

void getV965Hist(const char *name, int chNum=0)
{
    

    double qdc[32];
    sprintf(buff, "%s.root", name);
    if (gSystem->AccessPathName(buff))
    {
        readV965(name);
    }
    TFile *file = new TFile(buff, "read");
    TTree *tree = (TTree *)file->Get("tree");
    tree->SetBranchAddress("qdc", qdc);

    TH1D *hq = new TH1D("hq", "200fC", 4000, 0, 4000);
    TH1D *hq2 = new TH1D("hq2", "25fC", 4000, 0, 4000);

    int N = tree->GetEntries();
    for (int i = 0; i < N; i++)
    {
        tree->GetEntry(i);
        hq->Fill(qdc[chNum]);
        hq2->Fill(qdc[chNum+16]);
    }
    DrawMyHist(hq, "qdc (200fC)", "Counts", 1, 2);
    DrawMyHist(hq2, "qdc (25fC)", "Counts", 1, 2);
    TCanvas *c1 = cdC(1);
    hq->Draw();
    TCanvas *c2 = cdC(2);
    hq2->Draw();
}
*/
void readQDC(const char *name,const int totalNum=16)
{
    
    double qdc[totalNum];

    sprintf(buff, "%s.root", name);
    TFile *file = new TFile(buff, "recreate");
    TTree *tree = new TTree("tree", "");
    sprintf(buff,"qdc[%d]/D",totalNum);
    tree->Branch("qdc", qdc, buff);

    sprintf(buff, "%s.txt", name);
    ifstream input;
    input.open(buff);
    if (!input)
    {
        cout << "file can't be found!" << endl;
        return;
    }
   

    while (input && (!input.eof()))
    {
        for (int i = 0; i < totalNum; i++)
            input >> qdc[i];
        tree->Fill();
    }
    input.close();
    if (file->IsOpen())
    {
        tree->Write();
    }
}

void getHist(const char *name,const int totalNum, int chNum=0)
{
    
    string LSB[3] = {"200fC","25fC","100fC"};
    double qdc[totalNum];
    sprintf(buff, "%s.root", name);
    if (gSystem->AccessPathName(buff))
    {
        readQDC(name,totalNum);
        sprintf(buff, "%s.root", name);
    }
    TFile *file = new TFile(buff, "read");
    TTree *tree = (TTree *)file->Get("tree");
    tree->SetBranchAddress("qdc", qdc);
    if(totalNum>16) sprintf(buff,"qdc (%s)",LSB[0].data());
    else sprintf(buff,"qdc (%s)",LSB[2].data());
    TH1D *hq = new TH1D("hq", buff, 4000, 0, 4000);
    TH1D *hq2 = new TH1D("hq2", "qdc (25fC)", 4000, 0, 4000);

    int N = tree->GetEntries();
    for (int i = 0; i < N; i++)
    {
        tree->GetEntry(i);
        hq->Fill(qdc[chNum]);
        if (totalNum>16)
        hq2->Fill(qdc[chNum+16]);
    }
    DrawMyHist(hq, buff, "Counts", 1, 2);
    DrawMyHist(hq2, "qdc (25fC)", "Counts", 1, 2);
    TCanvas *c1 = cdC(1);
    hq->Draw();
    TCanvas *c2 = cdC(2);
    if (totalNum>16) hq2->Draw();
}