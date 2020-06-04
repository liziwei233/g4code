#include "Include/DrawMyClass.h"
#define MyClass_cxx
#include "MyClass3ch.h"
//#include "MyClass4ch.h"

//char path[1024] = "/mnt/d/Experiment/labtest/XGS_MCP-PMT/12-2";
//char path[1024] = "/mnt/f/MCP/12-22";
char path[1024] = "/mnt/f/LPZ/R3809-HVSCAN";
//char name[1024] = "201901-A3";
char name[1024] = "1244-3100v";
TH1F *hq = new TH1F("hq", ";charge (pC);Counts",11e3, -1, 10);
TH1F *hq2 = new TH1F("hq2", ";charge (pC);Counts", 11e3, -1, 10);

TH1F *ha = new TH1F("ha", ";Amp(V);Counts", 1e3, -1, 200e-3);
TH1F *hr = new TH1F("hr", ";risetime (ns);Counts", 1e3, 0, 1);
TH1F *ht = new TH1F("ht", ";time (ns);Counts", 50e3, 0, 50); // 1ps/bin
TH1F *hbl = new TH1F("hbl", ";baseline (V);Counts", 400, -10e-3, 10e-3);
TH1F *hblrms = new TH1F("hblrms", ";baselineRMS (V);Counts", 2e3, 0, 20e-3);
TH1F *hctratio = new TH1F("hctratio", ";Crosstalk ratio;Counts", 11e3, -1, 10);
/*
TH1F *hq = new TH1F("hq", ";charge (pC);Counts", 2e3, -1, 20);
TH1F *hq2 = new TH1F("hq2", ";charge (pC);Counts", 2e3, -1, 20);

TH1F *ha = new TH1F("ha", ";Amp(V);Counts", 1e3, -1, 1);
TH1F *hr = new TH1F("hr", ";risetime (ns);Counts", 1e3, 0, 1);
TH1F *ht = new TH1F("ht", ";time (ns);Counts", 50e3, 0, 50); // 1ps/bin
TH1F *hbl = new TH1F("hbl", ";baseline (V);Counts", 400, -10e-3, 10e-3);
TH1F *hblrms = new TH1F("hblrms", ";baselineRMS (V);Counts", 2e3, 0, 20e-3);
TH1F *hctratio = new TH1F("hctratio", ";Crosstalk ratio;Counts", 11e3, -1, 10);
 */
TTree *t1 = new TTree();

void getrootname(const char *rootname = "201901-A3")
{
    sprintf(name,"%s",rootname);
    cout<< "your file name is: "<<name<<endl; 
}


void gethist()
{
    TGaxis::SetMaxDigits(3);

    double fcharge;
    double fcharge2;
    double frise;
    double fbaseline;
    double fbaselinerms;
    double famplitude;
    double fampnerbor;
    double ftime;
    double freftime;
    hq->Reset();
    hq2->Reset();  
    ha->Reset();
    hr->Reset();
    ht->Reset();
    hbl->Reset();
    hblrms->Reset();
    hctratio->Reset();
    double riseth = -999;   //Unit: ns.
    double chargeth = 0.02; //Unit: pC.
    sprintf(buff, "%s/%s.root", path, name);
    MyClass t(buff);
    t1 = (TTree *)t.fChain;
    int N = t1->GetEntries();
    for (int i = 0; i < N; i++)
    {
        t1->GetEntry(i);
        fcharge =       t.MCP2_all_charge[0];
        fcharge2 =      t.MCP2_all_charge[1];
        frise =         t.MCP2_rise_time;
        fbaseline =     t.MCP2_baseline_level;
        fbaselinerms =  t.MCP2_baseline_rms;
        famplitude =    t.MCP2_global_maximum_y;
        ftime =         t.MCP2_CFDtime[3];

        freftime =      t.MCP3_CFDtime[3];
        fampnerbor =    t.MCP1_global_maximum_y;

        if (frise > riseth)
        {
            hq->Fill(fcharge);
            hq2->Fill(fcharge2);
            ha->Fill(famplitude);
            //cout<<Q[1]<<endl;
            hbl->Fill(fbaseline);
            hblrms->Fill(fbaselinerms);
            if (fcharge > chargeth)
            {
                hr->Fill(frise);
                hctratio->Fill(fampnerbor/famplitude);
                if(frise<0.4) ht->Fill(ftime - freftime);
            }
        }
    }
}

double pmtfun(double *x, double *par)
{
    Double_t val = 0.;
    Double_t amp = par[0];
    Double_t lambda = par[1];
    Double_t ped = par[2];
    Double_t pedSigma = par[3];
    Double_t peak = par[4];
    Double_t sigma = par[5];
    Double_t bw = par[6];

    //TF1 *myGaus = new TF1("myGaus","[0]*exp(-pow((x-[1])/[2],2)/2.)/sqrt(2.*TMath::Pi()*pow([2],2))",0,4000);
    //myGaus->SetParameters(TMath::PoissonI(0,lambda), ped, pedSigma);
    //val += myGaus->Eval(x[0]);
    val += TMath::PoissonI(0, lambda) * TMath::Gaus(x[0], ped, pedSigma, kTRUE);

    for (int i = 1; i < 4; i++)
    {
        //myGaus->SetParameters(TMath::PoissonI(i,lambda), ped+i*peak, sigma*sqrt(i));// CHANGED!!
        //val += myGaus->Eval(x[0]);
        val += TMath::PoissonI(i, lambda) * TMath::Gaus(x[0], ped + i * peak, sigma * sqrt(i), kTRUE);
    }

    //delete myGaus;
    return amp * val * bw;
}
/*
double HVfun(double *x, double *par)
{
    double val = 0.;
    double A = par[1];
    double alpha = par[2];
    double C = par[0];
    //val = A*TMath::Power(x[0],beta);
    val = TMath::Exp(C + x[0] * A * alpha);
    return val;
}
 */
TH1 *SPSfit(TH1 *h, int rbq, RANGE u, double fac)
{
    TH1 *hqdc = (TH1 *)h->Clone();
    hqdc->Draw();
    hqdc->Rebin(rbq);
    hqdc->SetLineColor(1);
    TF1 *myGaus = new TF1("myGaus", "gaus", u.L, u.R);
    int ibin = hqdc->GetMaximumBin();
    double pedMean = hqdc->GetBinCenter(ibin);
    double pedfitrangeleft = hqdc->GetBinCenter(ibin - 1);
    double pedfitrangeright = hqdc->GetBinCenter(ibin + 1);
    cout << "pedstal first fit range:" << pedfitrangeleft << "\t" << pedfitrangeright << endl;
    hqdc->GetXaxis()->SetRangeUser(pedfitrangeleft, pedfitrangeright);
    //
    //* try to fit pedstal
    hqdc->Fit(myGaus, "", "", pedfitrangeleft, pedfitrangeright); //???????????
    pedMean = myGaus->GetParameter(1);
    double pedSigma = myGaus->GetParameter(2);
    //return hqdc;

    //
    //* find the position of SPE
    hqdc->GetXaxis()->SetRangeUser(pedMean + fac * pedSigma, u.R);
    ibin = hqdc->GetMaximumBin();
    double mean = hqdc->GetBinCenter(ibin) - pedMean;
    //double sigma = hqdc->GetStdDev()/10;
    double sigma = hqdc->GetStdDev();
    double mean2 = hqdc->GetMean() - pedMean;
    if (mean2 / mean > 5)
        mean = mean2;

    myGaus->SetParameter(1, mean);
    myGaus->SetParameter(2, sigma);
    hqdc->Fit(myGaus, "R", "", pedMean + fac * pedSigma, pedMean + fac * pedSigma + 2 * mean);
    mean = myGaus->GetParameter(1);
    sigma = myGaus->GetParameter(2);
    cout << " init. par.: SPEmean = " << mean << "; SPEsigma = " << sigma << endl;
    //return hqdc;

    hqdc->GetXaxis()->SetRangeUser(u.L, u.R);

    TF1 *myFun = new TF1("myFun", "pmtfun", u.L, u.R, 7);
    //TF1 *myFun = new TF1("myFun", myfunc, &LZWfunc::mcpfun, u.L, u.R, 8);
    myFun->SetParNames("N", "#lambda", "#mu_{ped}", "#sigma_{ped}", "#mu", "#sigma", "BW");

    //Int_t ibin =  h->GetMaximumBin();
    myFun->SetParameters(hqdc->GetEntries(), 0.01, pedMean, pedSigma, mean, 0.8*sigma, hqdc->GetBinWidth(5));
    myFun->FixParameter(0, hqdc->GetEntries());                       //fix total yield
    myFun->FixParameter(6, hqdc->GetBinWidth(5));                     //fix bin width
    myFun->SetParLimits(2, pedMean - 0.0001, pedMean + 0.0001);       //fix pedestal mean
    myFun->SetParLimits(3, pedSigma - 0.000001, pedSigma + 0.000001); //fix pedestal sigma
    //myFun->SetParLimits(4, pedMean+fac*pedSigma, 5*mean);//>10 for 1400V

    //myFun->FixParameter(4, mean); //>10 for 1400V
    myFun->SetParLimits(4, mean - pedMean, 1.5 * (mean - pedMean)); //>10 for 1400V
    myFun->SetParLimits(5, 0.5 * sigma, 1. * sigma);
    //if(pedSigma>5) myFun->SetRange(pedMean+12*pedSigma,  u.R);
    //else
    //myFun->SetRange(pedMean+fac*pedSigma,  mean+1*sigma);
    myFun->SetRange(pedMean + fac * pedSigma, pedMean + mean + 1 * sigma);
    cout << " first fitting..." << endl;
    hqdc->Fit(myFun, "R");
    //return hqdc;
    mean = myFun->GetParameter(4); //return;
    sigma = myFun->GetParameter(5);
    //myFun->SetRange(pedMean+mean-1.3*TMath::Abs(sigma),  u.R);
    myFun->SetRange(pedMean + fac * pedSigma, pedMean + 3 * mean);
    //myFun->SetParLimits(4,mean/2., u.R);
    myFun->SetParameter(4, mean);
    myFun->SetParameter(5, sigma);
    myFun->SetParLimits(4, 0.9*(mean-pedMean), 1.2*(mean-pedMean));
    myFun->SetParLimits(5, 0.5*sigma, 1. * sigma);

    //myFun->SetRange(u.L,  u.R);
    //cout<<" second fitting..."<<endl;
    //myFun->SetParLimits(4,mean/2., u.R);
    hqdc->Fit(myFun, "R");
    //return hqdc;
    cout << " fitting done" << endl;
    mean = myFun->GetParameter(4);
    sigma = myFun->GetParameter(5);
    Double_t par[10], parErr[10];
    for (int j = 0; j < 7; j++)
    {
        par[j] = myFun->GetParameter(j);
        parErr[j] = myFun->GetParError(j);
    }

    Double_t xmin = u.L, xmax = 0;
    if ((pedMean - 15 * pedSigma) > 0)
        xmin = pedMean - 12 * pedSigma;
    if (par[1] > 2)
        hqdc->GetXaxis()->SetRangeUser(xmin, pedMean + 10 * mean + 3*sigma);
    else if (par[1] > 1)
        hqdc->GetXaxis()->SetRangeUser(xmin, pedMean + 7 * mean + 3*sigma);
    else if (par[1] > 0.1)
        hqdc->GetXaxis()->SetRangeUser(xmin, pedMean + 5 * mean + 3*sigma);
    else if (par[1] > 0.05)
        hqdc->GetXaxis()->SetRangeUser(xmin, pedMean + 4 * mean + 3*sigma);
    else
        hqdc->GetXaxis()->SetRangeUser(xmin, pedMean + 4 * mean + 3*sigma);

    hqdc->SetLineWidth(1.5);
    myGaus->SetLineWidth(2);
    myGaus->SetLineColor(4);
    //myGaus->Draw("same");

    myFun->SetNpx(1000);
    //myFun->SetLineWidth(1);
    myFun->SetLineColor(1);
    //myFun->SetLineStyle(7);
    myFun->GetRange(xmin, xmax);
    myFun->SetRange(0, xmax);
    myFun->DrawCopy("same");
    myFun->SetLineStyle(1);
    myFun->SetLineWidth(3);
    myFun->SetRange(pedMean - 4 * pedSigma, u.R);
    //myFun->Draw("same");

    TF1 *peakFun[10];
    for (int j = 1; j < 4; j++)
    {
        sprintf(buff, "peakFun%d", j);
        peakFun[j] = new TF1(buff, "[0]*exp(-pow((x-[1])/[2],2)/2.)/sqrt(2.*TMath::Pi()*pow([2],2))", u.L, u.R);
        //peakFun[j] = new TF1(buff,"gaus",u.L,u.R);
        //if(j==0) peakFun[j]->SetParameters(par[1]*par[0]*par[7],par[3]+par[5], par[6]);
        //else
        peakFun[j]->SetParameters(TMath::PoissonI(j, par[1]) * par[0] * par[6], par[2] + j * par[4], par[5] * sqrt(1.0 * j)); //CHANGED!!
        //if(j==0) peakFun[j]->SetLineColor(6);
        //else
        peakFun[j]->SetLineColor(kRed + 1);
        peakFun[j]->SetLineWidth(2);
        peakFun[j]->Draw("same");
    }

    myFun->Draw("same");
    return hqdc;
}
TF1 *gausfit(TH1 *h, double sigma, double facleft, double facright, int rbU, double UL, double UR)
{
    double mean = 0;
    //double sigma = 0;
    TH1 *hU = (TH1 *)h->Clone();
    hU->Draw();
    hU->Rebin(rbU);
    hU->GetXaxis()->SetRangeUser(UL + 1, UR - 1);
    mean = hU->GetBinCenter(hU->GetMaximumBin());
    //sigma = hU->GetRMS();
    TF1 *fitU = new TF1("fitU", "gaus", mean - facleft * sigma, mean + facright * sigma);
    hU->GetXaxis()->SetRangeUser(mean - facleft * sigma, mean + facright * sigma);
    cout << mean << "\t" << sigma << endl;

    // fitU->SetParLimits(0, 0,hU->GetMaximum()*1.1);
    fitU->SetParameter(1, mean);
    hU->Fit(fitU, "Q");
    mean = fitU->GetParameter(1);
    sigma = fitU->GetParameter(2);

    cout << mean << "\t" << sigma << endl;
    //return NULL;
    TFitResultPtr failed = hU->Fit(fitU, "Q", "", mean - facleft * sigma, mean + facright * sigma);
    //failed =1 means fit failed
    if (failed)
        //return hU = NULL;
        return fitU = NULL;

    else
    {

        if (UL < mean - 20 * sigma)
            UL = mean - 20 * sigma;
        if (UR > mean + 20 * sigma)
            UR = mean + 20 * sigma;

        hU->GetXaxis()->SetRangeUser(UL, UR);
        //return hU;
        return fitU;
    }
}
void drawMPE(int CanvasNum=1,double Gain=1.939e6)
{
    //if (!hq->GetEntries())
        gethist();
    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    c1->SetLogy();
    DrawMyHist(hq, "", "",1,3);
    TF1 *fq;
    fq = gausfit(hq,1,1.8,1.5,40,-1,8);
    
    double MGain = fq->GetParameter(1)* 1e-12 / 1.6e-19;
    double NPE = MGain/Gain;
    //double Gain=fq->GetParameter(4)*1e-12/1.6e-19;
    cout << "NPE=" << NPE << endl;
    sprintf(buff, "MGain=%0.2e", MGain);
    TLatex *l = DrawMyLatex(buff, 0.3, 0.6);
    l->Draw();
    sprintf(buff, "NPE=%0.2f", NPE);
    l = DrawMyLatex(buff, 0.3, 0.2);
    l->Draw();
    sprintf(buff, "%s/%sMPEcharge.png", path,name);
    c1->SaveAs(buff);
}
void drawMPE2(int CanvasNum=1,double NPE=18.43,double leftrange=2,double rightrange=1.2,double sigma = 1)
{
    //if (!hq->GetEntries())
        gethist();
    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    c1->SetLogy();
    DrawMyHist(hq, "", "",1,3);
    TF1 *fq;
    fq = gausfit(hq,1,leftrange,rightrange,40,-1,8);
    
    double MGain = fq->GetParameter(1)* 1e-12 / 1.6e-19;
    double Gain = MGain/NPE;
    //double Gain=fq->GetParameter(4)*1e-12/1.6e-19;
    cout << "Gain=" << Gain << endl;
    sprintf(buff, "MGain=%0.2e", MGain);
    TLatex *l = DrawMyLatex(buff, 0.3, 0.6);
    l->Draw();
    sprintf(buff, "Gain=%0.2e", Gain);
    l = DrawMyLatex(buff, 0.3, 0.2);
    l->Draw();
    sprintf(buff, "%s/%sMPE2charge.png", path,name);
    c1->SaveAs(buff);
}
void drawSPE(int CanvasNum=1,double rangefac=10)
{
    //if (!hq->GetEntries())
        gethist();
    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    c1->SetLogy();
    RANGE qrange = {-1, 12};
    TH1F *hqfit = (TH1F *)SPSfit(hq, 8, qrange, rangefac);
    DrawMyHist(hqfit, "", "",1,3);

    TF1 *fq = hqfit->GetFunction("myFun");
    double Gain = (fq->GetParameter(4) - fq->GetParameter(2)) * 1e-12 / 1.6e-19;
    //double Gain=fq->GetParameter(4)*1e-12/1.6e-19;
    cout << "Gain=" << Gain << endl;
    sprintf(buff, "Gain=%0.2e", Gain);
    TLatex *l = DrawMyLatex(buff, 0.3, 0.6);
    l->Draw();
    sprintf(buff, "%s/%scharge.png", path,name);
    c1->SaveAs(buff);
}
void drawSPE2(int CanvasNum=1,double rangefac=5)
{
    if (!hq2->GetEntries())
        gethist();
    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    c1->SetLogy();
    RANGE qrange = {-1, 12};
    TH1F *hqfit = (TH1F *)SPSfit(hq2, 4, qrange, rangefac);
    DrawMyHist(hqfit, "", "",1,3);

    TF1 *fq = hqfit->GetFunction("myFun");
    double Gain = (fq->GetParameter(4) - fq->GetParameter(2)) * 1e-12 / 1.6e-19;
    //double Gain=fq->GetParameter(4)*1e-12/1.6e-19;
    cout << "Gain=" << Gain << endl;
    sprintf(buff, "Gain=%0.2e", Gain);
    TLatex *l = DrawMyLatex(buff, 0.3, 0.6);
    l->Draw();
    sprintf(buff, "%s/%scharge_1.png", path,name);
    c1->SaveAs(buff);
}
void drawTR(int CanvasNum=1,double fac=0.2)
{
    if (!ht->GetEntries())
    {

        gethist();
        cout << "Get hist ......" << endl;
    }
    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    ht->Draw();
    double tL = ht->GetBinCenter(ht->GetMaximumBin() - 800);
    double tR = ht->GetBinCenter(ht->GetMaximumBin() + 800);
    //cout<<"maximumbin: "<<ht->GetMaximumBin()<<endl;
    cout<<"set tL & tR: "<< tL <<"\t"<< tR<<endl;
    TH1F *htfit = (TH1F *)twogausfit(ht, fac, 4, 16, tL, tR);
    DrawMyHist(htfit, "", "",1,3);
    TF1 *f = (TF1 *)htfit->GetFunction("fit2");
    double sigma = f->GetParameter(2) * 1e3;
    sprintf(buff, "#sigma=%.0fps", sigma);
    TLatex *l = DrawMyLatex(buff, 0.2, 0.5);
    l->Draw();
    sprintf(buff, "%s/%sTR.png", path,name);
    c1->SaveAs(buff);
}
void drawrise(int CanvasNum=1)
{
    if (!hr->GetEntries())
    {

        gethist();
        cout << "Get hist ......" << endl;
    }
    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    hr->Draw();
    double tL = hr->GetBinCenter(hr->GetMaximumBin() - 600);
    double tR = hr->GetBinCenter(hr->GetMaximumBin() + 600);
    //cout<<"maximumbin: "<<ht->GetMaximumBin()<<endl;
    //cout<< tL <<"\t"<< tR<<endl;
    TH1F *hrfit = (TH1F *)gausfit(hr, 2.5, 1.8, 4, tL, tR);
    DrawMyHist(hrfit, "", "",1,3);
    TF1 *f = (TF1 *)hrfit->GetFunction("fitU");
    double mean = f->GetParameter(1) * 1e3;
    sprintf(buff, "Risetime=%.0fps", mean);
    TLatex *l = DrawMyLatex(buff, 0.55, 0.3);
    l->Draw();
    sprintf(buff, "%s/%sRisetime.png", path,name);
    c1->SaveAs(buff);
}

void drawblrms(int CanvasNum=1)
{
    if (!hblrms->GetEntries())
    {

        gethist();
        cout << "Get hist ......" << endl;
    }
    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    hblrms->Draw();
    double tL = hblrms->GetBinCenter(hblrms->GetMaximumBin() - 600);
    double tR = hblrms->GetBinCenter(hblrms->GetMaximumBin() + 600);
    //cout<<"maximumbin: "<<ht->GetMaximumBin()<<endl;
    //cout<< tL <<"\t"<< tR<<endl;
    TH1F *hblrmsfit = (TH1F *)gausfit(hblrms, 2.5, 1.8, 2, tL, tR);
    DrawMyHist(hblrmsfit, "", "",1,3);
    TF1 *f = (TF1 *)hblrmsfit->GetFunction("fitU");
    double mean = f->GetParameter(1) * 1e3;
    sprintf(buff, "BaselineRMS=%.2fmV", mean);
    TLatex *l = DrawMyLatex(buff, 0.55, 0.3);
    l->Draw();
    sprintf(buff, "%s/%sblrms.png", path,name);
    c1->SaveAs(buff);
}
void drawctratio(int CanvasNum=1)
{
        if (!hctratio->GetEntries())
    {

        gethist();
        cout << "Get hist ......" << endl;
    }
    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    hctratio->Draw();
    double tL = hctratio->GetBinCenter(hctratio->GetMaximumBin() - 100);
    double tR = hctratio->GetBinCenter(hctratio->GetMaximumBin() + 100);
    cout<<"maximumbin: "<<ht->GetMaximumBin()<<endl;
    cout<< tL <<"\t"<< tR<<endl;
    TH1F *hctratiofit = (TH1F *)gausfit(hctratio, 2.5, 1.2, 1, tL, tR);
    //return;
    DrawMyHist(hctratiofit, "", "",1,3);
    hctratiofit->GetXaxis()->SetNdivisions(505);
    TF1 *f = (TF1 *)hctratiofit->GetFunction("fitU");
    double mean = f->GetParameter(1);
    double sigma = f->GetParameter(2);
    sprintf(buff, "CrosstalkRatio=%.1f%%", mean*100);
    TLatex *l = DrawMyLatex(buff, 0.55, 0.3);
    l->Draw();

    Drawxline(mean+5*sigma,3,7,2);
    double purity = hctratio->Integral(0,hctratio->FindBin(mean+5*sigma))/hctratio->Integral();
    
    sprintf(buff, "SignalPurity=%.2f%%", purity*100);
    TLatex *l2 = DrawMyLatex(buff, 0.55, 0.45);
    l2->Draw();
    sprintf(buff, "%s/%sCrosstalkratio.png", path,name);
    c1->SaveAs(buff);

}

void drawall()
{
    drawctratio(1);
    drawblrms(2);
    drawrise(3);
    drawTR(4);
    drawSPE(5);
}