#include "TH1F.h"
#include "Include/DrawMyClass.h"
#define MyClass_cxx
//#include "MyClass3ch.h"
#include "MyClass.h"
//#include "MyClassnew4ch.h"

//char path[1024] = "/mnt/d/Experiment/labtest/XGS_MCP-PMT/12-2";
//char path[1024] = "/mnt/f/XOPtest/4Anode/Rb";
//char path[1024] = "/mnt/f/R10754/KT0881/A7A8";
//char path[1024] = "/mnt/f/R10754/KT0881/A7light_A4crosstalk";
//char path[1024] = "/mnt/f/R10754/KT0881/DarkNoise";
//char path[1024] = "/data2/R710liziwei/lab_data/R10754DATA/HVSCAN";
//char path[1024] = "/mnt/d/ExpDATA/labtest/R10754/KT0881_HVSCAN";
//char path[1024] = "/mnt/d/ExpDATA/labtest/R10754/code/root_results";
//char path[1024] = "/mnt/f/wypmttest/pmt";
//char path[1024] = "/mnt/f/R10754/Gain";
char path[1024] = "/mnt/f/R10754/batch/MT0197/V1742";
//char path[1024] = "/mnt/f/R10754/batch/KT0881/V1742";
//char path[1024] = "/mnt/f/R10754/batch/KT0890/V1742";
//char path[1024] = "/mnt/c/Users/liziwei/OneDrive/work/g4code/data";
//char path[1024] = "/mnt/f/R10754/KT0881/newbase-crosstalkstudy";
//char path[1024] = "/mnt/f/R10754/KT0881-A4-HVScan";
//char path[1024] = "/mnt/e/R10754DATA2";
//char path[1024] = "/mnt/e/R10754DATA2/AMP_v2";

//char path[1024] = "/mnt/f/XOPtest/XOP-1MCP-SingleAnode2";
//char path[1024] = "/mnt/f/LPZ/Sr90-EJ228";
//char name[1024] = "201901-A3";
//char name[1024] = "1244-3100v";
//char name[1024] = "HV2100-A4-DN-close470pF";
//char name[1024] = "2000V-th30mV";
//char name[1024] = "HV2000-A4-Light4";
//char name[1024] = "HV2200-D3.3";
char name[1024] = "calibrate-HV2200-A1A9";
//char name[1024] = "KT0881-HV1900";
//char name[1024] = "KT0890-HV1800";
//char name[1024] = "KT0890-HV1800-OSC-A1A9";

//char name[1024] =       "HV1900-ADL5545-220nH";
//char name[1024] = "HV1900-D>3.3";
char timename[1024] = "HV1900-D3.3_timeinfo";
char datname[1024] = "HV1900-D3.3_timeinfo_3.9mVTOT";
//char name[1024] = "HV2200-D3.3";
//char timename[1024] = "HV2200-D3.3_timeinfo";
//char datname[1024] = "HV2200-D3.3_timeinfo_50mVTOT";

//char name[1024] =       "XOP-1MCP-SA2-1200-3";
//char name[1024] =       "HV2000-A4-noamp-D";
//char timename[1024] =   "HV2000-A4-noamp-D_timeinfo";
//char datname[1024] =    "HV2000-A4-noamp-D_timeinfo_10mVTOT";
//char name[1024] = "HV2200-ADL5545-220nH-10k_1";
//char name[1024] = "HV2200-TRF37D73-220nH-10k";
//char name[1024] = "newbase-HV2000-A4-noamp";
//char name[1024] = "XOP-1MCP-1500-20mV";
//char name[1024] = "newbase-HV2100-SPE-half-D3V-copper-A4A7-";
//TH1F *hq = new TH1F("hq", ";charge (pC);Counts", 70e3, -0.1, 100);
int HV=1800;
char TubeNo[1024]="KT0890";

TH1F *hq = new TH1F("hq", ";charge (pC);Counts", 1e3, -0.1, 2);
TH1F *hq2 = new TH1F("hq2", ";charge (pC);Counts", 15e3, -1, 8);

TH1F *ha = new TH1F("ha", ";Amp(V);Counts", 30e3, -1e-3, 3000e-3);
TH1F *hr = new TH1F("hr", ";risetime (ns);Counts", 2e3, 0.001, 2);
TH1F *hw = new TH1F("hw", ";FWHM (ns);Counts", 1e3, 0.001, 2);
TH1F *ht = new TH1F("ht", ";time (ns);Counts", 105e3, -5, 100); // 1ps/bin
//TH1F *ht = new TH1F("ht", ";time (ns);Counts", 105e3, -5, 30); // 1ps/bin

TH1F *hbl = new TH1F("hbl", ";baseline (V);Counts", 20.e3, -100e-3, 100e-3);
TH1F *hblrms = new TH1F("hblrms", ";baselineRMS (V);Counts", 2e3, 0, 20e-3);
TH1F *hctratio = new TH1F("hctratio", ";Crosstalk ratio;Counts", 11e3, -1, 10);
TH1F *hringratio = new TH1F("hringratio", ";Ringing ratio;Counts", 11e3, -1, 10);
TH1F *hinvertringratio = new TH1F("hinvertringratio", ";invertRinging ratio;Counts", 11e3, -1, 10);
/*
TH1F *hq = new TH1F("hq", ";charge (pC);Counts", 2e3, -1, 20);

//#ifdef GATE
TH1F* hqgate[3];

//#endif

TH1F *ha = new TH1F("ha", ";Amp(V);Counts", 1e3, -1, 1);
TH1F *hr = new TH1F("hr", ";risetime (ns);Counts", 1e3, 0, 1);
TH1F *ht = new TH1F("ht", ";time (ns);Counts", 50e3, 0, 50); // 1ps/bin
TH1F *hbl = new TH1F("hbl", ";baseline (V);Counts", 400, -10e-3, 10e-3);
TH1F *hblrms = new TH1F("hblrms", ";baselineRMS (V);Counts", 2e3, 0, 20e-3);
TH1F *hctratio = new TH1F("hctratio", ";Crosstalk ratio;Counts", 11e3, -1, 10);
 */
TTree *t1 = new TTree();
bool fexist;
double CN = 0;
double risethmin = -11110.1;  //Unit: ns.
double risethmax = 111110.35; //Unit: ns.
double chargethmin = 0.06;     //Unit: pC.
//chargethmax =0.55; //Unit: pC.
double chargethmax = 30; //Unit: pC.
Color_t clr[] = {1, 2, kGreen + 3, 4, 6, 7, kOrange, kViolet + 2};

void setallname(const int theHV=1800,const char* theTubeNo="KT0890"){
    HV=theHV;
    sprintf(TubeNo, "%s", theTubeNo);
    //sprintf(path, "/mnt/f/R10754/batch/%s/V1742", TubeNo);
    sprintf(path, "/mnt/f/R10754/batch/%s/OSC", TubeNo);
    cout << "your path is: " << path << endl;
    //sprintf(name, "%s-HV%d", TubeNo,HV);
    sprintf(name, "%s-HV%d-OSC-A1A9", TubeNo,HV);
    cout << "your file name is: " << name << endl;
}

bool gethist(double chargethmax = 1e4)
{
    TGaxis::SetMaxDigits(3);

    double fcharge;
    double fchargegate[3];
    double frise;
    double fbaseline;
    double fbaselinerms;
    double famplitude;
    double fampnerbor;
    double finvring;
    double finvampnerbor;
    double ftime;
    double freftime;
    double fctx;
    double fwidth;
    hq->Reset();
    hq2->Reset();
    ha->Reset();
    hr->Reset();
    hw->Reset();
    ht->Reset();
    hbl->Reset();
    hblrms->Reset();
    hctratio->Reset();
    hringratio->Reset();
    hinvertringratio->Reset();

    //double chargethmin = 300; //Unit: V1742 bin.
    //double chargethmin = 0.052; //Unit: pC.
    //double chargethmin = 0.5; //Unit: pC.

    //double chargethmax = 0.04; //Unit: pC.
    double blrmsth = 2;
    sprintf(buff, "%s/%s.root", path, name);
    if (gSystem->AccessPathName(buff))
    {
        cout << "Error!! The File " << buff << " doesn't exist" << endl;
        return 0;
    }
    MyClass t(buff);
    t1 = (TTree *)t.fChain;
    int N = t1->GetEntries();
    for (int i = 0; i < N; i++)
    {
        t1->GetEntry(i);
        fcharge = t.MCP2_all_charge[0];
        //fcharge2 =    t.MCP2_all_charge[1];
        frise = t.MCP2_rise_time[0];
        fwidth = t.MCP2_width;
        fbaseline = t.MCP2_baseline_level;
        fbaselinerms = t.MCP2_baseline_rms;
        famplitude = t.MCP2_global_maximum_y;
        ftime = t.MCP2_CFDtime[3];

        fampnerbor = t.MCP4_global_maximum_y;
        finvring = t.MCP4_secondinvertpeak_y;
        finvampnerbor = t.MCP4_invert_maximum_y;
        fctx = t.MCP4_global_maximum_x;
        //fampnerbor = t.MCP4_secondinvertpeak_y;
        freftime = t.TR1_CFDtime[6];

        if (fbaselinerms < blrmsth)
        //if (1)
        {
            if (frise > risethmin && frise < risethmax)
                hq->Fill(fcharge);

            ha->Fill(famplitude);
            //cout<<Q[1]<<endl;
            hbl->Fill(fbaseline);
            hblrms->Fill(fbaselinerms);
            if (fcharge > chargethmin && fcharge < chargethmax)
            {
                hr->Fill(frise);
                hw->Fill(fwidth);
                //if(t.MCP1_invert_maximum_y/t.MCP1_global_maximum_y>-0.2)
                if (frise > risethmin && frise < risethmax)
                {
                    ht->Fill(ftime - freftime);

                    //if (1)
                    //if (fctx > 40 && fctx < 43)
                    if (fctx > 23 && fctx < 27)
                    {
                        //if(fctx>37.5&&fctx<40&&fampnerbor<4e-3){
                        //if(fctx>63&&fctx<70){
                        //if(1){

                        hringratio->Fill(fampnerbor / famplitude);
                        hinvertringratio->Fill(-1 * finvring / famplitude);
                        hctratio->Fill(-1 * finvampnerbor / famplitude);
                    }
                }
            }
        }
    }
    cout << " The File " << buff << " build hists complete" << endl;
    return 1;
}

void drawbaselinestablity(int start = 0, int step = 2e4)
{
    setgStyle();
    //sprintf(name, "%s", rootname);
    double fbaseline;
    double fbaselinerms;
    //hbl->Reset();
    //hblrms->Reset();
    sprintf(buff, "%s/%s.root", path, name);
    if (gSystem->AccessPathName(buff))
    {
        cout << "Error!! The File " << buff << " doesn't exist" << endl;
        return;
    }
    MyClass t(buff);
    t1 = (TTree *)t.fChain;
    int N = t1->GetEntries();
    //int step = 2e4;
    //int start = 0;
    vector<double> vx;    // 0.66min/1e4 waveform
    vector<double> vxerr; // 0.66min/1e4 waveform
    vector<double> vbl;
    vector<double> vblerr;
    vector<double> vblrms;
    vector<double> vblrmserr;
    while (start < N)
    {

        //for (int i = start; i < start + step && i < N; i++)
        for (int i = start; i < start + step && i < N; i++)
        {
            t1->GetEntry(i);
            //fcharge2 =    t.MCP2_all_charge[1];
            fbaseline = t.MCP2_baseline_level;
            fbaselinerms = t.MCP2_baseline_rms;
            hbl->Fill(fbaseline);
            hblrms->Fill(fbaselinerms);
        }
        TH1 *hblfit = gausfit(hbl, 2e-3, 3, 3, 1, -15e-3, 15e-3);
        TF1 *f1 = (TF1 *)hblfit->GetFunction("fitU");
        double mean1 = f1->GetParameter(1) * 1e3;
        double sigma1 = f1->GetParameter(2) * 1e3;

        TH1 *hblrmsfit = gausfit(hblrms, 1e-3, 3, 3, 1, 0e-3, 10e-3);
        TF1 *f2 = (TF1 *)hblrmsfit->GetFunction("fitU");
        double mean2 = f2->GetParameter(1) * 1e3;
        double sigma2 = f2->GetParameter(2) * 1e3;

        start = start + step;
        vx.push_back(start / 1e4 * 0.66);
        vxerr.push_back(0);
        vbl.push_back(mean1);
        vblerr.push_back(sigma1);
        vblrms.push_back(mean2);
        vblrmserr.push_back(sigma2);
        //return;
    }
    TGraphErrors *g1 = new TGraphErrors(vx.size(), &vx[0], &vblrms[0], &vxerr[0], &vblrmserr[0]);
    TGraphErrors *g2 = new TGraphErrors(vx.size(), &vx[0], &vbl[0], &vxerr[0], &vblerr[0]);
    TCanvas *c = cdC(0);
    DrawMyGraph(g1, "Time (min)", "baselineRMS (mV)");
    g1->Draw("AP");
    g1->GetYaxis()->SetRangeUser(0, 2);
    sprintf(buff, "%s/%sblRMSstability_MCP2.png", path, name);
    c->SaveAs(buff);

    c = cdC(1);
    DrawMyGraph(g2, "Time (min)", "baseline (mV)");
    g2->Draw("AP");
    g2->GetYaxis()->SetRangeUser(-2, 0);
    sprintf(buff, "%s/%sblstability_MCP2.png", path, name);
    c->SaveAs(buff);
}
void setrootname(const char *rootname = "201901-A3")
{
    sprintf(name, "%s", rootname);
    cout << "your file name is: " << name << endl;
    gethist();
    cout << "Get hist ... >> " << name << endl;
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

void MCPfit(double c, double alpha, double theta)
{
    TF1 *f1 = new TF1("f", "[0]*ROOT::Math::gamma_pdf(x,[1],[2])", 0, 100);
    //f1->SetNDF(10e4);
    f1->SetNpx(10e3);
    f1->SetParameter(0, c);
    f1->SetParameter(1, alpha);
    f1->SetParameter(2, theta);
    f1->Draw();
}

double mcpfun(double *x, double *par)
{
    Double_t val = 0.;
    Double_t amp = par[0];
    Double_t lambda = par[1];
    Double_t ped = par[2];
    Double_t pedSigma = par[3];
    Double_t gamma = par[4]; //shape parameter
    Double_t x0 = par[5];    // location parameter
    Double_t beta = par[6];  // scale parameter
    Double_t bw = par[7];

    //TF1 *myGaus = new TF1("myGaus","[0]*exp(-pow((x-[1])/[2],2)/2.)/sqrt(2.*TMath::Pi()*pow([2],2))",0,4000);
    //myGaus->SetParameters(TMath::PoissonI(0,lambda), ped, pedSigma);
    //val += myGaus->Eval(x[0]);
    val += TMath::PoissonI(0, lambda) * TMath::Gaus(x[0], ped, pedSigma, kTRUE);

    for (int i = 1; i < 4; i++)
    {
        //myGaus->SetParameters(TMath::PoissonI(i,lambda), ped+i*peak, sigma*sqrt(i));// CHANGED!!
        //val += myGaus->Eval(x[0]);
        val += TMath::PoissonI(i, lambda) * TMath::GammaDist(x[0], gamma * sqrt(i), ped + i * x0, beta);
    }

    //delete myGaus;
    return amp * val * bw;
}
TH1 *SPSMCPfit(TH1 *h, int rbq, RANGE u, double leftfac, double rightfac)

{
    TH1 *hqdc = (TH1 *)h->Clone();
    hqdc->Draw();
    hqdc->Rebin(rbq);
    hqdc->SetLineColor(1);
    TF1 *myGaus = new TF1("myGaus", "gaus", u.L, u.R);
    TF1 *myGamma = new TF1("myGamma", "[0]*TMath::GammaDist(x,[1],[2],[3])", 0, u.R);
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
    //hqdc->GetXaxis()->SetRangeUser(pedMean + leftfac * pedSigma, u.R);
    hqdc->GetXaxis()->SetRangeUser(pedMean + leftfac * pedSigma, pedMean + rightfac * pedSigma);
    ibin = hqdc->GetMaximumBin();
    double mean = hqdc->GetBinCenter(ibin) - pedMean;
    //double sigma = hqdc->GetStdDev()/10;
    double sigma = hqdc->GetStdDev();
    double mean2 = hqdc->GetMean() - pedMean;
    if (mean2 / mean > 5)
        mean = mean2;

    myGamma->SetParameter(1, mean * mean / sigma);
    myGamma->SetParameter(2, 0);
    myGamma->SetParameter(3, sigma / mean);
    hqdc->Fit(myGamma, "R", "", pedMean + leftfac * pedSigma, pedMean + rightfac * pedSigma);
    double alpha = myGamma->GetParameter(1);
    double beta = myGamma->GetParameter(3);
    Drawxline(pedMean + leftfac * pedSigma);
    Drawxline(pedMean + rightfac * pedSigma);
    cout << " init. par.: SPEmean = " << alpha * beta << "; SPEsigma = " << alpha * beta * beta << endl;
    //return hqdc;

    hqdc->GetXaxis()->SetRangeUser(u.L, u.R);

    TF1 *myFun = new TF1("myFun", "mcpfun", u.L, u.R, 8);
    //TF1 *myFun = new TF1("myFun", myfunc, &LZWfunc::mcpfun, u.L, u.R, 8);
    myFun->SetParNames("N", "#lambda", "#mu_{ped}", "#sigma_{ped}", "x0", "#Gamma", "#beta", "BW");

    //Int_t ibin =  h->GetMaximumBin();
    myFun->SetParameters(hqdc->GetEntries(), 0.01, pedMean, pedSigma, alpha, 0, beta, hqdc->GetBinWidth(5));
    myFun->FixParameter(0, hqdc->GetEntries());                       //fix total yield
    myFun->FixParameter(7, hqdc->GetBinWidth(5));                     //fix bin width
    myFun->SetParLimits(2, pedMean - 0.0001, pedMean + 0.0001);       //fix pedestal mean
    myFun->SetParLimits(3, pedSigma - 0.000001, pedSigma + 0.000001); //fix pedestal sigma
    //myFun->SetParLimits(4, pedMean+fac*pedSigma, 5*mean);//>10 for 1400V

    //myFun->FixParameter(4, mean); //>10 for 1400V
    myFun->SetParLimits(4, 1, 2 * alpha);                                   //>10 for 1400V
    myFun->SetParLimits(5, pedMean - 3 * pedSigma, pedMean + 3 * pedSigma); //>10 for 1400V
    myFun->SetParLimits(6, 0.5 * beta, 1. * beta);
    //if(pedSigma>5) myFun->SetRange(pedMean+12*pedSigma,  u.R);
    //else
    //myFun->SetRange(pedMean+fac*pedSigma,  mean+1*sigma);
    //myFun->SetRange(pedMean + leftfac * pedSigma, pedMean + mean + 1 * sigma);
    cout << " first fitting..." << endl;
    hqdc->Fit(myFun, "R");
    return hqdc;
    mean = myFun->GetParameter(4); //return;
    sigma = myFun->GetParameter(5);
    //myFun->SetRange(pedMean+mean-1.3*TMath::Abs(sigma),  u.R);
    myFun->SetRange(pedMean + leftfac * pedSigma, pedMean + 3 * mean);
    //myFun->SetParLimits(4,mean/2., u.R);
    myFun->SetParameter(4, mean);
    myFun->SetParameter(5, sigma);
    myFun->SetParLimits(4, 0.9 * (mean - pedMean), 1.2 * (mean - pedMean));
    myFun->SetParLimits(5, 0.5 * sigma, 1. * sigma);

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
        hqdc->GetXaxis()->SetRangeUser(xmin, pedMean + 10 * mean + 3 * sigma);
    else if (par[1] > 1)
        hqdc->GetXaxis()->SetRangeUser(xmin, pedMean + 7 * mean + 3 * sigma);
    else if (par[1] > 0.1)
        hqdc->GetXaxis()->SetRangeUser(xmin, pedMean + 5 * mean + 3 * sigma);
    else if (par[1] > 0.05)
        hqdc->GetXaxis()->SetRangeUser(xmin, pedMean + 4 * mean + 3 * sigma);
    else
        hqdc->GetXaxis()->SetRangeUser(xmin, pedMean + 4 * mean + 3 * sigma);

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
        //peakFun[j] = new TF1(buff, "[0]*exp(-pow((x-[1])/[2],2)/2.)/sqrt(2.*TMath::Pi()*pow([2],2))", u.L, u.R);
        peakFun[j] = new TF1(buff, "[0]*TMath::GammaDist(x, [2], [1], 1)", par[2] + j * par[4], u.R);
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
TH1 *SPSfit(TH1 *h, int rbq, RANGE u, double leftfac, double rightfac)

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
    //hqdc->GetXaxis()->SetRangeUser(pedMean + leftfac * pedSigma, u.R);
    hqdc->GetXaxis()->SetRangeUser(pedMean + leftfac * pedSigma, pedMean + rightfac * pedSigma);
    ibin = hqdc->GetMaximumBin();
    double mean = hqdc->GetBinCenter(ibin) - pedMean;
    //double sigma = hqdc->GetStdDev()/10;
    double sigma = hqdc->GetStdDev();
    double mean2 = hqdc->GetMean() - pedMean;
    if (mean2 / mean > 5)
        mean = mean2;

    myGaus->SetParameter(1, mean);
    myGaus->SetParameter(2, sigma);
    hqdc->Fit(myGaus, "R", "", pedMean + leftfac * pedSigma, pedMean + rightfac * pedSigma);
    mean = myGaus->GetParameter(1);
    sigma = myGaus->GetParameter(2);
    Drawxline(pedMean + leftfac * pedSigma);
    Drawxline(pedMean + rightfac * pedSigma);
    cout << " init. par.: SPEmean = " << mean << "; SPEsigma = " << sigma << endl;
    //return hqdc;

    hqdc->GetXaxis()->SetRangeUser(u.L, u.R);

    TF1 *myFun = new TF1("myFun", "pmtfun", u.L, u.R, 7);
    //TF1 *myFun = new TF1("myFun", myfunc, &LZWfunc::mcpfun, u.L, u.R, 8);
    myFun->SetParNames("N", "#lambda", "#mu_{ped}", "#sigma_{ped}", "#mu", "#sigma", "BW");

    //Int_t ibin =  h->GetMaximumBin();
    myFun->SetParameters(hqdc->GetEntries(), 0.01, pedMean, pedSigma, mean, 0.8 * sigma, hqdc->GetBinWidth(5));
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
    myFun->SetRange(pedMean + leftfac * pedSigma, pedMean + mean + 1 * sigma);
    cout << " first fitting..." << endl;
    hqdc->Fit(myFun, "R");
    //return hqdc;
    mean = myFun->GetParameter(4); //return;
    sigma = myFun->GetParameter(5);
    //myFun->SetRange(pedMean+mean-1.3*TMath::Abs(sigma),  u.R);
    myFun->SetRange(pedMean + leftfac * pedSigma, pedMean + 3 * mean);
    //myFun->SetParLimits(4,mean/2., u.R);
    myFun->SetParameter(4, mean);
    myFun->SetParameter(5, sigma);
    myFun->SetParLimits(4, 0.9 * (mean - pedMean), 1.2 * (mean - pedMean));
    myFun->SetParLimits(5, 0.5 * sigma, 1. * sigma);

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
        hqdc->GetXaxis()->SetRangeUser(xmin, pedMean + 10 * mean + 3 * sigma);
    else if (par[1] > 1)
        hqdc->GetXaxis()->SetRangeUser(xmin, pedMean + 7 * mean + 3 * sigma);
    else if (par[1] > 0.1)
        hqdc->GetXaxis()->SetRangeUser(xmin, pedMean + 5 * mean + 3 * sigma);
    else if (par[1] > 0.05)
        hqdc->GetXaxis()->SetRangeUser(xmin, pedMean + 4 * mean + 3 * sigma);
    else
        hqdc->GetXaxis()->SetRangeUser(xmin, pedMean + 4 * mean + 3 * sigma);

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

/*
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
 */
void drawMPE(int CanvasNum = CN++, double Gain = 5.54e5, int reb = 40)
{
    if (!hq->GetEntries())
        gethist();
    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    c1->SetLogy();
    //hq->Draw();
    //DrawMyHist(hq, "", "", 1, 3);
    //TF1 *fq;

    TH1F *hqfit = (TH1F *)gausfit(hq, 20, 1.8, 1.5, reb, 0.5, 150);
    TF1 *fq = (TF1 *)hqfit->GetFunction("fitU");
    double MGain = fq->GetParameter(1) * 1e-12 / 1.6e-19;
    double NPE = MGain / Gain;
    //double Gain=fq->GetParameter(4)*1e-12/1.6e-19;
    cout << "NPE=" << NPE << endl;
    sprintf(buff, "MGain=%0.2e", MGain);
    TLatex *l = DrawMyLatex(buff, 0.3, 0.6);
    l->Draw();
    sprintf(buff, "NPE=%0.2f", NPE);
    l = DrawMyLatex(buff, 0.3, 0.2);
    l->Draw();

    sprintf(buff, "%s/%sMPEcharge.png", path, name);
    c1->SaveAs(buff);
}
void drawMPE2(int CanvasNum = CN++, double NPE = 20.80, double leftrange = 2, double rightrange = 1.2, double sigma = 1)
{
    if (!hq->GetEntries())
        gethist();
    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    c1->SetLogy();
    DrawMyHist(hq, "", "", 1, 3);
    //TF1 *fq;
    TH1F *hqfit = (TH1F *)gausfit(hq, 1, leftrange, rightrange, 40, -1, 8);
    TF1 *fq = (TF1 *)hqfit->GetFunction("fitU");
    double MGain = fq->GetParameter(1) * 1e-12 / 1.6e-19;
    double Gain = MGain / NPE;
    //double Gain=fq->GetParameter(4)*1e-12/1.6e-19;
    cout << "Gain=" << Gain << endl;
    sprintf(buff, "MGain=%0.2e", MGain);
    TLatex *l = DrawMyLatex(buff, 0.3, 0.6);
    l->Draw();
    sprintf(buff, "Gain=%0.2e", Gain);
    l = DrawMyLatex(buff, 0.3, 0.2);
    l->Draw();
    sprintf(buff, "%s/%sMPE2charge.png", path, name);
    c1->SaveAs(buff);
}
void drawAMP(int CanvasNum = CN++, int reb = 20, double leftfac = 10, double rightfac = 100)
{
    if (!ha->GetEntries())
        gethist();
    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    c1->SetLogy();
    //    ha->Draw();
    //    DrawMyHist(ha, "", "", 1, 3);

    RANGE arange = {-1e-3, 8000e-3};
    TH1F *hafit = (TH1F *)SPSfit(ha, reb, arange, leftfac, rightfac);
    DrawMyHist(hafit, "", "", 1, 3);
    TF1 *fq = hafit->GetFunction("myFun");
    double SPEAmp = fq->GetParameter(4) - fq->GetParameter(2);
    //double Gain=fq->GetParameter(4)*1e-12/1.6e-19;
    cout << "SPEAmp=" << SPEAmp << endl;
    sprintf(buff, "SPEAmp=%0.2fmV", SPEAmp * 1e3);
    TLatex *l = DrawMyLatex(buff, 0.3, 0.75);
    l->Draw();

    sprintf(buff, "%s/%samplitude.png", path, name);
    c1->SaveAs(buff);
}

double drawSPE(int CanvasNum = CN++, int reb = 2, double leftfac = 20, double rightfac = 200)
{

    if (!hq->GetEntries())
    {
        fexist = gethist();
        if (!fexist)
            return 0;
    }
    setgStyle();
    //return;
    TCanvas *c1 = cdC(CanvasNum);
    c1->SetLogy();
    RANGE qrange = {-0.051, 80};
    TH1F *hqfit = (TH1F *)SPSfit(hq, reb, qrange, leftfac, rightfac);
    //TH1F *hqfit = (TH1F *)SPSMCPfit(hq, reb, qrange, leftfac, rightfac);
    DrawMyHist(hqfit, "", "", 1, 3);
    TF1 *fq = hqfit->GetFunction("myFun");
    double Gain = (fq->GetParameter(4) - fq->GetParameter(2)) * 1e-12 / 1.6e-19;
    //double Gain=fq->GetParameter(4)*1e-12/1.6e-19;
    cout << "Gain=" << Gain << endl;
    sprintf(buff, "Gain=%0.2e", Gain);
    TLatex *l = DrawMyLatex(buff, 0.3, 0.75);
    l->Draw();
    sprintf(buff, "%s/%sC%dcharge.png", path, name,CanvasNum);
    c1->SaveAs(buff);
    return Gain;
}
void drawSPE2(int CanvasNum = CN++, double rangefac = 5)
{
    if (!hq->GetEntries())
        gethist();
    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    c1->SetLogy();
    RANGE qrange = {-0.051, 12};
    TH1F *hqfit = (TH1F *)SPSfit(hq2, 4, qrange, rangefac, rangefac + 10);
    DrawMyHist(hqfit, "", "", 1, 3);

    TF1 *fq = hqfit->GetFunction("myFun");
    double Gain = (fq->GetParameter(4) - fq->GetParameter(2)) * 1e-12 / 1.6e-19;
    //double Gain=fq->GetParameter(4)*1e-12/1.6e-19;
    cout << "Gain=" << Gain << endl;
    sprintf(buff, "Gain=%0.2e", Gain);
    TLatex *l = DrawMyLatex(buff, 0.3, 0.6);
    l->Draw();
    sprintf(buff, "%s/%scharge_1.png", path, name);
    c1->SaveAs(buff);
}
void drawTR(int CanvasNum = CN++, double tL = 18, double tR = 23)
{
    if (!ht->GetEntries())
    {

        gethist();
        cout << "Get hist ......" << endl;
    }
    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    ht->Draw();
    //double tL = ht->GetBinCenter(ht->GetMaximumBin() - 200);
    //double tR = ht->GetBinCenter(ht->GetMaximumBin() + 200);
    //cout<<"maximumbin: "<<ht->GetMaximumBin()<<endl;
    cout << "set tL & tR: " << tL << "\t" << tR << endl;
    bool fitfailed = 1;
    int counter = 0;
    double fac[10] = {0.5, 0.55, 0.49, 0.45, 0.4};
    TF1 *fitT;
    TF1 *fittr;
    TF1 *fitbg;
    TH1F *htfit;
    double TRratio;
    double BGratio;

    while (fitfailed)
    {
        //htfit = (TH1F *)twogausfit(ht, 0.01, 4, 6, 8, tL, tR);
        htfit = (TH1F *)twogausfit(ht, fac[counter], 4, 6, 8, tL, tR);
        if (!htfit)
        {
            cout << "The twogausfit failed!" << endl;
            sprintf(buff, "%s/%sTR.png", path, name);
            c1->SaveAs(buff);
            return;
        }
        fitT = (TF1 *)htfit->GetFunction("fit2");
        if (fitT->GetChisquare() / fitT->GetNDF() < 2 || counter >= 5)
            fitfailed = 0;
        counter++;
    }
    fittr = new TF1("fittr", "gaus", tL, tR);
    fittr->SetParameter(0, fitT->GetParameter(0));
    fittr->SetParameter(1, fitT->GetParameter(1));
    fittr->SetParameter(2, fitT->GetParameter(2));
    fitbg = new TF1("fitbg", "gaus", tL, tR);
    fitbg->SetParameter(0, fitT->GetParameter(3));
    fitbg->SetParameter(1, fitT->GetParameter(4));
    fitbg->SetParameter(2, fitT->GetParameter(5));
    TRratio = fittr->Integral(tL, tR) / (fittr->Integral(tL, tR) + fitbg->Integral(tL, tR));
    BGratio = fitbg->Integral(tL, tR) / (fittr->Integral(tL, tR) + fitbg->Integral(tL, tR));
    //htfit->GetXaxis()->SetRangeUser(30,32.53);
    //TF1 *f = (TF1 *)htfit->GetFunction("fit2");
    double sigma = fitT->GetParameter(2) * 1e3;
    //sprintf(buff, "#sigma=%.0fps", sigma);
    sprintf(buff, "#sigma=%.0fps,%.2f%%", sigma, TRratio * 100);
    TLatex *l = DrawMyLatex(buff, 0.2, 0.5);
    l->Draw();
    DrawMyHist(htfit, "", "", 1, 3);
    htfit->SetNdivisions(505);
    sprintf(buff, "%s/%sTR.png", path, name);
    c1->SaveAs(buff);
}
void drawrise(int CanvasNum = CN++, float leftfac = 2., float rightfac = 2.0)
{
    if (!hr->GetEntries())
    {

        fexist = gethist();
        cout << "Get hist ......" << endl;
        if (!fexist)
            return;
        if (!hr)
        {
            cout << "the rise hist is NULL " << endl;
            return;
        }
    }
    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    hr->Draw();
    //double tL = 0.01;
    //double tR = 0.3;
    double tL = hr->GetBinCenter(hr->GetMaximumBin() - 2e2);
    double tR = hr->GetBinCenter(hr->GetMaximumBin() + 2e2);
    cout << "maximumbin: " << hr->GetMaximumBin() << endl;
    cout << tL << "\t" << tR << endl;
    TH1F *hrfit = (TH1F *)gausfit(hr, 0.1, leftfac, rightfac, 3, tL, tR);
    if (!hrfit)
    {
        cout << "the risefit hist is NULL " << endl;
        return;
    }
    DrawMyHist(hrfit, "", "", 1, 3);
    TF1 *f = (TF1 *)hrfit->GetFunction("fitU");
    double mean = f->GetParameter(1) * 1e3;
    sprintf(buff, "Risetime=%.0fps", mean);
    TLatex *l = DrawMyLatex(buff, 0.55, 0.3);
    l->Draw();
    sprintf(buff, "%s/%sRisetime.png", path, name);
    c1->SaveAs(buff);
}
void drawwidth(int CanvasNum = CN++, float leftfac = 2.5, float rightfac = 1.0)
{
    if (!hw->GetEntries())
    {

        fexist = gethist();
        cout << "Get hist ......" << endl;
        if (!fexist)
            return;
        if (!hw)
        {
            cout << "the rise hist is NULL " << endl;
            return;
        }
    }

    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    hw->Draw();
    double tL = hw->GetBinCenter(hw->GetMaximumBin() - 300);
    double tR = hw->GetBinCenter(hw->GetMaximumBin() + 300);
    cout << "maximumbin: " << ht->GetMaximumBin() << endl;
    cout << tL << "\t" << tR << endl;
    TH1F *hwfit = (TH1F *)gausfit(hw, 0.1, leftfac, rightfac, 4, tL, tR);
    if (!hwfit)
    {
        cout << "the risefit hist is NULL " << endl;
        return;
    }
    DrawMyHist(hwfit, "", "", 1, 3);
    TF1 *f = (TF1 *)hwfit->GetFunction("fitU");
    double mean = f->GetParameter(1) * 1e3;
    sprintf(buff, "FWHM=%.0fps", mean);
    TLatex *l = DrawMyLatex(buff, 0.55, 0.3);
    l->Draw();
    sprintf(buff, "%s/%swidth.png", path, name);
    c1->SaveAs(buff);
}
void drawbl(int CanvasNum = CN++, int rb = 10)
{
    if (!hbl->GetEntries())
    {

        gethist();
        cout << "Get hist ......" << endl;
    }
    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    hbl->Draw();
    double tL = hbl->GetBinCenter(hbl->GetMaximumBin() - 600 * rb);
    double tR = hbl->GetBinCenter(hbl->GetMaximumBin() + 600 * rb);
    //cout<<"maximumbin: "<<ht->GetMaximumBin()<<endl;
    //cout<< tL <<"\t"<< tR<<endl;
    TH1F *hblfit = (TH1F *)gausfit(hbl, 0.1, 2.5, 1.8, rb, tL, tR);
    DrawMyHist(hblfit, "", "", 1, 3);
    TF1 *f = (TF1 *)hblfit->GetFunction("fitU");
    double mean = f->GetParameter(1) * 1e3;
    sprintf(buff, "Baseline=%.2fmV", mean);
    TLatex *l = DrawMyLatex(buff, 0.55, 0.3);
    l->Draw();
    sprintf(buff, "%s/%sbl.png", path, name);
    c1->SaveAs(buff);
}

void drawblrms(int CanvasNum = CN++)
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
    TH1F *hblrmsfit = (TH1F *)gausfit(hblrms, 0.1, 2.5, 1.8, 2, tL, tR);
    DrawMyHist(hblrmsfit, "", "", 1, 3);
    TF1 *f = (TF1 *)hblrmsfit->GetFunction("fitU");
    double mean = f->GetParameter(1) * 1e3;
    sprintf(buff, "BaselineRMS=%.2fmV", mean);
    TLatex *l = DrawMyLatex(buff, 0.55, 0.3);
    l->Draw();
    sprintf(buff, "%s/%sblrms.png", path, name);
    c1->SaveAs(buff);
}
void drawctratio(int CanvasNum = CN++)
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
    double tR = hctratio->GetBinCenter(hctratio->GetMaximumBin() + 150);
    cout << "maximumbin: " << hctratio->GetMaximumBin() << endl;
    cout << "maximumbin pos: " << hctratio->GetBinCenter(hctratio->GetMaximumBin()) << endl;
    cout << tL << "\t" << tR << endl;
    TH1F *hctratiofit = (TH1F *)gausfit(hctratio, 0.001, 2, 1.2, 1, tL, tR);
    //return;
    DrawMyHist(hctratiofit, "", "", 1, 3);
    hctratiofit->GetXaxis()->SetNdivisions(505);
    TF1 *f = (TF1 *)hctratiofit->GetFunction("fitU");
    double mean = f->GetParameter(1);
    double sigma = f->GetParameter(2);
    sprintf(buff, "CrosstalkRatio=%.1f%%", mean * 100);
    //hctratiofit->GetXaxis()->SetRangeUser(-0.1, 1);
    TLatex *l = DrawMyLatex(buff, 0.55, 0.3);
    l->Draw();

    //Drawxline(mean + 5 * sigma, 3, 7, 2);
    //double purity = hctratio->Integral(0, hctratio->FindBin(mean + 5 * sigma)) / hctratio->Integral();
    Drawxline(0.2, 3, 7, 2);
    double purity = hctratio->Integral(0, hctratio->FindBin(0.2)) / hctratio->Integral();

    sprintf(buff, "SignalPurity=%.2f%%", purity * 100);
    TLatex *l2 = DrawMyLatex(buff, 0.55, 0.45);
    l2->Draw();
    sprintf(buff, "%s/%sCrosstalkratio.png", path, name);
    c1->SaveAs(buff);
}
void drawringratio(int CanvasNum = CN++)
{
    if (!hringratio->GetEntries())
    {

        gethist();
        cout << "Get hist ......" << endl;
    }
    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    hringratio->Draw();
    double tL = hringratio->GetBinCenter(hringratio->GetMaximumBin() - 100);
    double tR = hringratio->GetBinCenter(hringratio->GetMaximumBin() + 150);
    cout << "maximumbin: " << hringratio->GetMaximumBin() << endl;
    cout << "maximumbin pos: " << hringratio->GetBinCenter(hringratio->GetMaximumBin()) << endl;
    cout << tL << "\t" << tR << endl;
    TH1F *hringratiofit = (TH1F *)gausfit(hringratio, 0.001, 2, 1.2, 1, tL, tR);
    //return;
    DrawMyHist(hringratiofit, "", "", 1, 3);
    hringratiofit->GetXaxis()->SetNdivisions(505);
    TF1 *f = (TF1 *)hringratiofit->GetFunction("fitU");
    double mean = f->GetParameter(1);
    double sigma = f->GetParameter(2);
    sprintf(buff, "RingingRatio=%.1f%%", mean * 100);
    //hringratiofit->GetXaxis()->SetRangeUser(-0.1, 1);
    TLatex *l = DrawMyLatex(buff, 0.55, 0.3);
    l->Draw();

    //Drawxline(mean + 5 * sigma, 3, 7, 2);
    //double purity = hringratio->Integral(0, hringratio->FindBin(mean + 5 * sigma)) / hringratio->Integral();
    Drawxline(0.2, 3, 7, 2);
    double purity = hringratio->Integral(0, hringratio->FindBin(0.2)) / hringratio->Integral();

    sprintf(buff, "SignalPurity=%.2f%%", purity * 100);
    TLatex *l2 = DrawMyLatex(buff, 0.55, 0.45);
    l2->Draw();
    sprintf(buff, "%s/%sRingingratio.png", path, name);
    c1->SaveAs(buff);
}

void drawinvringratio(int CanvasNum = CN++)
{
    if (!hinvertringratio->GetEntries())
    {

        gethist();
        cout << "Get hist ......" << endl;
    }
    setgStyle();
    TCanvas *c1 = cdC(CanvasNum);
    hinvertringratio->Draw();
    double tL = hinvertringratio->GetBinCenter(hinvertringratio->GetMaximumBin() - 100);
    double tR = hinvertringratio->GetBinCenter(hinvertringratio->GetMaximumBin() + 150);
    cout << "maximumbin: " << hinvertringratio->GetMaximumBin() << endl;
    cout << "maximumbin pos: " << hinvertringratio->GetBinCenter(hinvertringratio->GetMaximumBin()) << endl;
    cout << tL << "\t" << tR << endl;
    TH1F *hinvertringratiofit = (TH1F *)gausfit(hinvertringratio, 0.001, 1.2, 1.5, 1, tL, tR);
    //return;
    DrawMyHist(hinvertringratiofit, "", "", 1, 3);
    hinvertringratiofit->GetXaxis()->SetNdivisions(505);
    TF1 *f = (TF1 *)hinvertringratiofit->GetFunction("fitU");
    double mean = f->GetParameter(1);
    double sigma = f->GetParameter(2);
    sprintf(buff, "RingingRatio=%.1f%%", mean * 100);
    //hinvertringratiofit->GetXaxis()->SetRangeUser(-0.1, 1);
    TLatex *l = DrawMyLatex(buff, 0.55, 0.3);
    l->Draw();

    //Drawxline(mean + 5 * sigma, 3, 7, 2);
    //double purity = hinvertringratio->Integral(0, hinvertringratio->FindBin(mean + 5 * sigma)) / hinvertringratio->Integral();
    Drawxline(0.2, 3, 7, 2);
    double purity = hinvertringratio->Integral(0, hinvertringratio->FindBin(0.2)) / hinvertringratio->Integral();

    sprintf(buff, "SignalPurity=%.2f%%", purity * 100);
    TLatex *l2 = DrawMyLatex(buff, 0.55, 0.45);
    l2->Draw();
    sprintf(buff, "%s/%sInvRingingratio.png", path, name);
    c1->SaveAs(buff);
}
void drawavewaveform(int CanvasNum = CN++, float left = 2, float right = 12)
{
    int CHid1 = 2;
    int CHid2 = 4;
    char buff[1024];
    //const char *path = "/mnt/f/R10754/KT0881/A7A8";
    //const char *path = "/mnt/f/R10754/KT0881/A4light_A7crosstalk";
    sprintf(buff, "%s/%saverage.root", path, name);
    if (gSystem->AccessPathName(buff))
    {
        cout << "Error!! The File " << buff << " doesn't exist" << endl;
        return;
    }
    TFile *f1 = new TFile(buff, "read");
    sprintf(buff, "CH%daverage_waveform", CHid1);
    TGraph *g1 = (TGraph *)f1->Get(buff);
    sprintf(buff, "CH%daverage_waveform", CHid2);
    TGraph *g2 = (TGraph *)f1->Get(buff);
    //DrawMyGraph(g1, "Time (ns)", "Amplitude (V)", 0.1, 20, kGreen+2);
    g1->SetLineWidth(3);
    g2->SetLineWidth(3);
    g1->GetXaxis()->SetTitle("Time (ns)");
    g1->GetYaxis()->SetTitle("Amplitude (V)");
    g1->GetXaxis()->SetTitleFont(42);
    g1->GetYaxis()->SetTitleFont(42);
    g1->GetXaxis()->SetTitleSize(0.07);
    g1->GetYaxis()->SetTitleSize(0.07);
    g1->GetXaxis()->SetRangeUser(left, right);
    g2->GetXaxis()->SetRangeUser(left, right);
    //g2->GetXaxis()->SetRangeUser(0, 60);
    TCanvas *c1 = cdC(CanvasNum++);
    g1->Draw();
    sprintf(buff, "%s/%sWFsignal.png", path, name);
    gPad->SaveAs(buff);

    c1 = cdC(CanvasNum++);
    g2->Draw();
    sprintf(buff, "%s/%sWFct.png", path, name);
    gPad->SaveAs(buff);

    g2->SetLineColor(2);
    c1 = cdC(CanvasNum++);
    g1->Draw();
    g2->Draw("same");
    sprintf(buff, "%s/%sWFtogether.png", path, name);
    gPad->SaveAs(buff);
}
void drawDN(const char *path = "/mnt/f/R10754/KT0881/DarkNoiseScan")
{
    int HV[] = {1900, 1950, 2000, 2100, 2200};
    const int n = sizeof(HV) / sizeof(HV[0]);
    double min = 0;
    double max = 0;
    double th[100];
    double DN[100];
    char str[1024];
    char buff[1024];
    ifstream input;
    TGraph *g1;

    TCanvas *c1;
    for (int i = 1; i < 2; i++)
    {
        sprintf(buff, "%s/%d.txt", path, HV[i]);
        input.open(buff);
        if (!input)
        {
            cout << "file can't be found: \n"
                 << buff << endl;
            return 0;
        }
        int k = 0;
        while (!input.eof())
        {

            cout << "k=" << k << endl;
            input >> th[k] >> DN[k];
            cout << th[k] << "\t" << DN[k] << endl;
            k++;
        }
        max = DN[0] * 1.2;
        input.close();

        g1 = new TGraph(k, th, DN);
        c1 = cdC(i);
        SetMyPad(gPad, 0.15, 0.05, 0.1, 0.14);
        //DrawMyPad(gPad, "Threshold (mV)", "DarkNoiseRate (Hz) ", 0, 100, 0, 100, 0, 0);
        DrawMyGraph(g1, "Threshold (mV)", "DarkNoiseRate (Hz) ", 1.5, 20, kGreen + 2);
        g1->Draw("AP");
        g1->GetYaxis()->SetRangeUser(min, max);
        sprintf(buff, "%s/HV%d_DN.png", path, HV[i]);
        //sprintf(buff, "UV-6-HVscan.png");
        c1->SaveAs(buff);
    }
}
void drawHV()
{
    setgStyle();
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(1112);

    int N = 14;
    double gain[N];

    double x[N];

    char str[1024];
    char buff[1024];
    sprintf(buff, "%s/HVScan.txt", path);
    ifstream input;
    input.open(buff);
    if (!input)
    {
        cout << "file can't be found: \n"
             << buff << endl;
        return 0;
    }
    for (int i = 0; i < N; i++)
    {
        cout << "enter the file: " << buff << endl;
        input >> x[i] >> gain[i];
    }

    //double GainSyserr = 0.0257; //unit:pC
    double GainSyserr = 0.0; //unit:pC

    /*    
    for (int i = 0; i < n; i++)
    {
        yerr[i] = sqrt(yerr[i] * yerr[i] + GainSyserr * GainSyserr);
        yerr[i] = yerr[i] * 1.e-12 / 1.6e-19;
    }
*/
    TGraphErrors *g1 = new TGraphErrors(N, x, gain, 0, 0);

    TCanvas *c1;
    c1 = cdC(1);
    c1->SetLogy();
    SetMyPad(gPad, 0.15, 0.15, 0.1, 0.14);
    DrawMyPad(gPad, "Work voltage (kV)", "Gain ", 1800, 2400, 2e4, 1e8, 0, 0);
    g1->Draw("Psame");

    //mydraw.Graph(g1,"NPE","TimeRes (ps)",1.5,20,4);

    DrawMyGraph(g1, "Work voltage (kV)", "Gain ", 1.5, 20, kGreen + 2);

    TF1 *fhv = new TF1("fhv", HVfun, 1800, 2500, 3);
    fhv->SetParNames("cons", "#delta", "#alpha");
    //fhv->SetParLimits(0,1,1e7);
    //fhv->SetParLimits(1,1,10);
    fhv->SetParameter(0, -10);
    fhv->FixParameter(2, 40);
    g1->Fit(fhv, "", "", 1850, 2350);
    fhv->Draw("same");
    //g1->GetXaxis()->SetRangeUser(1800, 2900);
    //g1->GetYaxis()->SetRangeUser(1e4, 1e7);
    gPad->Update();
    gPad->Modified();

    sprintf(buff, "%s/_HVscan.png", path);
    //sprintf(buff, "UV-6-HVscan.png");
    c1->SaveAs(buff);
}

bool GetGain(int N = 4)
{
    TFile *sfile;
    sprintf(buff, "%s/%sgain_5ns.root", path, name);
    sfile = new TFile(buff, "recreate");

    double fcharge[4];
    TH1F *hq[N];

    for (int i = 0; i < N; i++)
    {
        sprintf(buff, "hq_Light%d", i);
        hq[i] = new TH1F(buff, ";charge (pC);Counts", 10e3, -1, 10);
        sprintf(buff, "%s/HV2000-A4-Light%d.root", path, i);
        if (gSystem->AccessPathName(buff))
        {
            cout << "Error!! The File " << buff << " doesn't exist" << endl;
            return 0;
        }
        cout << "====>>  Start to open the root file : " << buff << endl;
        TFile *f1 = new TFile(buff, "read");
        TTree *t2 = (TTree *)f1->Get("Pico");
        t2->SetMakeClass(1);

        t2->SetBranchAddress("MCP2_all_charge", fcharge);
        int NEvents = t2->GetEntries();
        for (int j = 0; j < NEvents; j++)
        {
            t2->GetEntry(j);
            hq[i]->Fill(fcharge[1]);
        }

        f1->Close();
        sfile->cd();
        //sfile->Write();
        sfile->WriteTObject(hq[i]);
    }
    return 1;
}

double gainfun(double *x, double *par)
{
    double val = 0;
    double lambda = par[0];
    double D0 = par[1];
    double D1 = par[2];
    //double D2 = par[3];
    //val = TMath::Poisson(0, x[0] * lambda) * D0 + TMath::Poisson(1, x[0] * lambda) * D1 + TMath::Poisson(2, x[0] * lambda) * D2;
    val = TMath::Poisson(0, x[0] * lambda) * D0 + TMath::Poisson(1, x[0] * lambda) * D1;
    return val;
}
Double_t pmtFun4(Double_t *x, Double_t *par)
{
    Double_t val = 0.;
    Double_t lama = par[0];
    Double_t f1 = par[1];
    Double_t f2 = par[2];
    Double_t N = par[3];
    Double_t Num0 = par[4];
    val += Num0 * TMath::Exp(-x[0] * lama) + N * x[0] * lama * f1 * TMath::Exp(-x[0] * lama) + 0.5 * N * x[0] * x[0] * lama * lama * f2 * TMath::Exp(-x[0] * lama);
    //val += Num0*TMath::Exp(-x[0]*lama)+N*x[0]*lama*f1*TMath::Exp(-x[0]*lama)+0.5*N*x[0]*lama*lama*f2*TMath::Exp(-x[0]*lama);
    return val;
}
void fitgain(int N = 4)
{
    setgStyle();
    TGaxis::SetMaxDigits(3);
    TH1F *hq[N];
    TCanvas *c;
    TLegend *leg = DrawMyLeg(0.6, 0.7, 0.85, 0.9);
    sprintf(buff, "%s/%sgain_5ns.root", path, name);
    if (gSystem->AccessPathName(buff))
    {
        fexist = GetGain();
        if (!fexist)
        {

            cout << "Error!! The source File of " << buff << " doesn't exist" << endl;
            return;
        }
    }
    cout << "====>>  Start to open the root file : " << buff << endl;
    TFile *f1 = new TFile(buff, "read");
    c = cdC(0);
    int histID[] = {0, 1, 2, 4, 6, 7};
    double lambda_initial[] = {0, 1.355, 0.263, 1.209, 1.5, 0.822};
    //double lambda_initial[]={0,1.355,0.263,0.815,1.209,4.5,1.5,0.822};
    N = sizeof(histID) / sizeof(histID[0]);
    for (int i = 0; i < N; i++)
    {
        sprintf(buff, "hq_Light%d", histID[i]);
        hq[i] = (TH1F *)f1->Get(buff);
        DrawMyHist(hq[i], "Charge (pC)", "Counts", clr[i], 3);
        leg->AddEntry(hq[i], buff, "l");
        if (i == 0)
        {

            hq[i]->Draw();
            hq[i]->GetXaxis()->SetRangeUser(-0.5, 4);
        }
        else
            hq[i]->Draw("same");
    }
    leg->Draw();
    c->SetLogy();
    sprintf(buff, "%s/%sgain.png", path, name);
    c->SaveAs(buff);

    double x[N - 1];
    double y[N - 1];
    int baseindex = 3;
    double D0;
    double thL = -0.05;
    double thR = 0.25;
    double th0 = 0;
    double th1 = 0;
    double thegain;
    double thegainerr;
    double lambdaerr;
    int sampleN = 100;
    TGraphErrors *gfitresult[N - 1];
    TGraphErrors *gfitlambda[N - 1];
    c = cdC(1);
    leg->Clear();
    //for(int k =0; k<1;k++)
    for (int k = 0; k < N - 1; k++)
    {
        gfitresult[k] = new TGraphErrors();
        gfitlambda[k] = new TGraphErrors();
        baseindex = k + 1;
        for (int j = 0; j < sampleN; j++)
        //for (int j = 0; j < 5; j++)
        {
            th0 = thL + (thR - thL) / sampleN * j;
            //th1 = thL + (thR - thL) / sampleN * (j+1);
            //D0 = hq[0]->Integral(hq[0]->FindBin(th0), hq[0]->FindBin(th1)) / hq[0]->Integral();
            D0 = hq[0]->Integral(0, hq[0]->FindBin(th0)) / hq[0]->Integral();
            cout << "D0= " << D0 << endl;
            for (int i = 0; i < N - 1; i++)
            {

                x[i] = (hq[i + 1]->GetMean() - hq[0]->GetMean()) / (hq[baseindex]->GetMean() - hq[0]->GetMean());
                //y[i] = hq[i + 1]->Integral(hq[i + 1]->FindBin(th0), hq[i + 1]->FindBin(th1)) / hq[i + 1]->Integral();
                //y[i] = hq[i + 1]->Integral(0, hq[i + 1]->FindBin(th0)) / hq[i + 1]->Integral();
                y[i] = hq[i + 1]->Integral(0, hq[i + 1]->FindBin(th0)) / hq[i + 1]->Integral() * hq[0]->Integral();
                //cout << "Point " << i << " ," << x[i] << "\t" << y[i] << endl;
            }
            TGraph *ggain = new TGraph(N - 1, x, y);
            DrawMyGraph(ggain, "#alpha", "cdf", 1, 8, 1);
            ggain->Draw("AP");
            ggain->SetMarkerStyle(8);
            //TF1 *fit1 = new TF1("fit1", "gainfun", 0, 5, 3);
            //fit1->SetParNames("#lambda", "D0", "D1", "D2");
            TF1 *fit1 = new TF1("fit1", "pmtFun4", 0, 5, 5);
            fit1->SetParNames("#lambda", "f1", "f2", "N", "N0");
            //fit1->SetParameter(0, 1.2);
            fit1->SetParameter(0, lambda_initial[baseindex]);
            fit1->SetParLimits(0, lambda_initial[baseindex] - 0.5, lambda_initial[baseindex] + 0.5);
            //fit1->FixParameter(1, D0);
            //fit1->FixParameter(2, 0);
            //fit1->SetParameter(2, 0.01);
            //fit1->SetParLimits(2, 0,0.1);
            //fit1->FixParameter(1, 0);
            fit1->SetParameter(1, 0.01);
            fit1->SetParLimits(1, 0, 0.1);
            fit1->FixParameter(2, 0);
            fit1->FixParameter(3, hq[0]->Integral());
            fit1->FixParameter(4, D0 * hq[0]->Integral());
            //fit1->SetParLimits(3, 0,0.1);
            ggain->Fit(fit1, "q");
            thegain = (hq[baseindex]->GetMean() - hq[0]->GetMean()) / fit1->GetParameter(0) * 1e-12 / 1.6e-19;
            //thegainerr=hq[baseindex]->GetMean()/fit1->GetParError(0)*1e-12/1.6e-19;
            thegainerr = 0;
            lambdaerr = fit1->GetParError(0);
            if (thegain > 100e6 || thegain < 0)
                thegain = 0;
            gfitresult[k]->SetPoint(j, th0, thegain);
            gfitresult[k]->SetPointError(j, 0, thegainerr);
            gfitlambda[k]->SetPoint(j, th0, fit1->GetParameter(0));
            gfitlambda[k]->SetPointError(j, 0, lambdaerr);
            cout << "Point " << j << " ," << th0 << "\t" << fit1->GetParameter(0) << ", " << fit1->GetParError(0) << endl;
            cout << "The fit Number:" << j << endl;
        }
        DrawMyGraph(gfitresult[k], "th (pC)", "Gain_{fit}", 1, 24, clr[k + 1]);
        DrawMyGraph(gfitlambda[k], "th (pC)", "#lambda_{fit}", 1, 24, clr[k + 1]);
        sprintf(buff, "#lambda of hq%d", k + 1);
        leg->AddEntry(gfitresult[k], buff, "lp");
    }
    c = cdC(2);
    //gfitresult->GetYaxis()->SetRangeUser(0.6,1.4e6);
    //N=2;
    //for(int j=1; j<2;j++){
    bool flag = 1;
    for (int j = 0; j < N - 1; j++)
    {

        if (flag)
        {
            gfitresult[j]->Draw("AP");
            flag = 0;
        }
        else
            gfitresult[j]->Draw("Psame");
        //if(gfitresult[1]) gfitresult[1]->Draw("Psame");
        //if(gfitresult[2]) gfitresult[2]->Draw("Psame");
    }
    cout << "program check" << endl;
    leg->Draw();
    sprintf(buff, "%s/%sgainvsth.png", path, name);
    c->SaveAs(buff);

    c = cdC(3);
    //gfitresult->GetYaxis()->SetRangeUser(0.6,1.4e6);
    flag = 1;
    for (int j = 0; j < N - 1; j++)
    {

        if (flag)
        {
            gfitlambda[j]->Draw("AP");
            flag = 0;
        }
        else
            gfitlambda[j]->Draw("Psame");
    }
    //gfitlambda[0]->Draw("AP");
    //if(gfitlambda[1]) gfitlambda[1]->Draw("Psame");
    //if(gfitlambda[2]) gfitlambda[2]->Draw("Psame");
    leg->Draw();
    sprintf(buff, "%s/%slambdavsth.png", path, name);
    c->SaveAs(buff);
}

bool GetV1742Gain(int rb=1,int left =20,int right=100, int N = 16,int startID=0)
{

    //setgStyle();
    sprintf(buff, "%s/%sV1742gain.dat", path, name);
    ofstream op;
    op.open(buff,ios::trunc);
    double fcharge[16][4];
    TH1F *hqarray[16];
    double fgain[16];
    double UL=-0.1;
    double UR=5;
    int Ubin=(UR-UL)/1e-3;
    
    sprintf(buff, "%s/%s.root", path, name);
    if (gSystem->AccessPathName(buff))
    {
        cout << "Error!! The File " << buff << " doesn't exist" << endl;
        return 0;
    }
    cout << "====>>  Start to open the root file : " << buff << endl;

    TFile *f1 = new TFile(buff, "read");
    TTree *t2 = (TTree *)f1->Get("Pico");
    t2->SetMakeClass(1);
    for (int i = startID; i < startID+N; i++)
    {
        sprintf(buff, "hqarray%d", i);
        hqarray[i] = new TH1F(buff, ";charge (pC);Counts", Ubin, UL, UR);
        sprintf(buff, "MCP%d_all_charge", i);
        t2->SetBranchAddress(buff, fcharge[i]);
    }
    int NEvents = t2->GetEntries();
    for (int i = 0; i < NEvents; i++)
    {
        t2->GetEntry(i);
        for (int j = startID; j < startID+N; j++)
        {
            hqarray[j]->Fill(fcharge[j][0]);
        }
    }
    for (int j = startID; j < startID+N; j++)
    {
        hq = hqarray[j];
        fgain[j]=drawSPE(j,rb,left,right);
        op<<fgain[j]<<endl;
    }
    
    //GainUniformtext->Draw("colz TEXT SAME");
    op.close();
    //f1->Close();
    return 1;
}
void drawBias(int N=16,int startID=0){
    double fgain[16];
    double mean;
    double var;
    TH2D* GainUniform= new TH2D("GainUniform","",4,1,5,4,1,5);
    ifstream input;
    sprintf(buff, "%s/%sgain.dat", path, name);
    if (gSystem->AccessPathName(buff))
    {
        cout << "Error!! The File " << buff << " doesn't exist" << endl;
        bool flag=GetV1742Gain();
        if(!flag) {
            cout << "Error!! Can't get the gain" << endl;
            return 0;
        }
    }
    sprintf(buff, "%s/%sgain.dat", path, name);
    input.open(buff);
    for(int j=0; j<N;j++)
    {
        input>>fgain[j];
        GainUniform->SetBinContent(j/4+1,j%4+1,URound(fgain[j]/1e6,2));
    }
    mean = TMath::Mean(N,&fgain[startID]);
    var = TMath::RMS(N,&fgain[startID]);
    //TH2D* GainUniformtext =(TH2D*) GainUniform->Clone();
    TCanvas *cg = cdC(100);
    TGaxis::SetMaxDigits(3);

    cg->SetGrid(0);
    DrawMy2dHist(GainUniform, "CHID X", "CHID Y");
    GainUniform->SetStats(0);
    sprintf(buff, "HV=%d,Gain(#times10^{6}),Bias=%.1f%%", HV,var/mean*100);
    gStyle->SetOptTitle(1);
    GainUniform->SetTitle(buff);
    GainUniform->Draw("colz");
    //TExec *ex1 = new TExec("ex1", "Pal();");
    //TExec *ex1 = new TExec("ex1", "gStyle->SetPalette(55, 0, 0.45);");
    //ex1->Draw();
    GainUniform->SetMarkerSize(2.4);
    GainUniform->SetMarkerColor(1);
    GainUniform->Draw("Colz TEXT SAME");
    sprintf(buff, "%s/%sgain.png", path, name);
    cg->SaveAs(buff);
}

bool GetOSCGain(int rb=1,int left =20,int right=100,int N=2)
{
    sprintf(buff, "%s/%sOSCgain.dat", path, name);
    ofstream op;
    op.open(buff,ios::trunc);
    //setgStyle();
    int CHID[]={2,4};
    int AnodeID[]={1,9};
    double fcharge[4][4];
    TH1F *hqarray[4];
    double fgain[4];
    double UL=-0.1;
    double UR=5;
    int Ubin=(UR-UL)/1e-3;
    double mean;
    double var;
    sprintf(buff, "%s/%s.root", path, name);
    if (gSystem->AccessPathName(buff))
    {
        cout << "Error!! The File " << buff << " doesn't exist" << endl;
        return 0;
    }
    cout << "====>>  Start to open the root file : " << buff << endl;

    TFile *f1 = new TFile(buff, "read");
    TTree *t2 = (TTree *)f1->Get("Pico");
    t2->SetMakeClass(1);
    for (int i = 0; i < N; i++)
    {
        sprintf(buff, "hqarray%d", AnodeID[i]);
        hqarray[i] = new TH1F(buff, ";charge (pC);Counts", Ubin, UL, UR);
        sprintf(buff, "MCP%d_all_charge", CHID[i]);
        t2->SetBranchAddress(buff, fcharge[CHID[i]-1]);
    }
    int NEvents = t2->GetEntries();
    for (int i = 0; i < NEvents; i++)
    {
        t2->GetEntry(i);
        for (int j = 0; j < N; j++)
        {
            hqarray[j]->Fill(fcharge[CHID[j]-1][0]);
        }
    }
    for (int j = 0; j < N; j++)
    {
        hq = hqarray[j];
        fgain[CHID[j]-1]=drawSPE(AnodeID[j],rb,left,right);
        op<<fgain[CHID[j]-1]<<endl;
    }
    //TH2D* GainUniformtext =(TH2D*) GainUniform->Clone();
   op.close();
    return 1;
}

bool Preprocesstimeinformation()
{
    double fLEDth[14];
    double fTOT[14];
    double fLEDtime[14];
    double freftime[8];
    double frise[4];
    double fcharge[4];
    double fwidth;
    double famp;

    double ftimediff[14];

    TFile *sfile;

    sprintf(timename, "%s_timeinfo", name);
    sprintf(buff, "%s/%s.root", path, timename);

    sfile = new TFile(buff, "recreate");
    TTree *t3 = new TTree("timeinfo", "");
    t3->Branch("LEDthrd", fLEDth, "LEDthrd[14]/D");
    t3->Branch("LEDTOT", fTOT, "LEDTOT[14]/D");
    t3->Branch("width", &fwidth, "width/D");
    t3->Branch("charge", &fcharge[0], "charge/D");
    t3->Branch("timediff", ftimediff, "timediff[14]/D");

    sprintf(buff, "%s/%s.root", path, name);
    if (gSystem->AccessPathName(buff))
    {
        cout << "Error!! The File " << buff << " doesn't exist" << endl;
        return 0;
    }
    cout << "====>>  Start to open the root file : " << buff << endl;

    TFile *f1 = new TFile(buff, "read");
    TTree *t2 = (TTree *)f1->Get("Pico");
    t2->SetMakeClass(1);
    t2->SetBranchAddress("MCP2_LEDthrd", fLEDth);
    t2->SetBranchAddress("MCP2_LEDtime", fLEDtime);
    t2->SetBranchAddress("MCP2_TOT", fTOT);
    t2->SetBranchAddress("MCP2_width", &fwidth);
    t2->SetBranchAddress("MCP2_rise_time", frise);
    t2->SetBranchAddress("MCP2_all_charge", fcharge);
    t2->SetBranchAddress("MCP2_global_maximum_y", &famp);
    t2->SetBranchAddress("TR1_CFDtime", freftime);

    //t3->Branch("reftime",&freftime[6],"reftime/D");
    //t3->Branch("rise",&frise[0],"rise/D");
    //t3->Branch("charge",&fcharge[0],"charge/D");

    int NEvents = t2->GetEntries();
    for (int i = 0; i < NEvents; i++)
    {
        t2->GetEntry(i);

        if (fcharge[0] > chargethmin && fcharge[0] < chargethmax &&
            frise[0] > risethmin && frise[0] < risethmax)
        {
            for (int k = 0; k < 14; k++)
            {
                ftimediff[k] = fLEDtime[k] - freftime[6];
            }
            t3->Fill();
        }

        //fLEDth[k] = t.MCP2_LEDthrd[k];
        //fLEDtime[k] = t.MCP2_LEDtime[k];
        //fTOT[k] = t.MCP2_TOT[k];
        //frise = t.MCP2_rise_time[0];
        //fcharge = t.MCP2_all_charge[0];
        //freftime = t.TR1_CFDtime[6];
        //cout<<fLEDth[k]<<"\t"<<fTOT[k]<<"\t"<<fLEDtime[k]<<"\t"<<frise[0]<<"\t"<<fcharge[0]<<"\t"<<freftime[6]<<endl;
        //if(fcharge[0]>5) cout<<fTOT[k]<<endl;
    }
    f1->Close();
    sfile->cd();
    sfile->WriteTObject(t3);
    sfile->Close();
    return 1;
}

void AnalysisTime(int start = 0)
{

    setgStyle();
    bool fexist;
    bool fitfailed;
    const int ThNum = 14;
    const int CorNo = 3;

    double fLEDth[14];
    double fLEDTOT[14];
    double ftimediff[14];
    double fwidth;
    double fcharge;

    int iter = 5;
    double fac[10] = {0.5, 0.55, 0.49, 0.45, 0.4}; //the ratio of base to signal
    double rangeL = 3;                             // the fit range factor
    double rangeR = 5;                             // the fit range factor
    //int rbt = 8;
    int rbt = 8;
    int rbU = 4;
    int counter = 0;
    //double UL = 4, UR = 60;
    double UL = -0.05, UR = 3;
    double tL = -1, tR = 100;

    double TAcor = 0;
    double T[400000] = {0};

    ofstream output;

    //sprintf(datname, "%s_19mVTOT", timename);
    sprintf(buff, "%s/%s.dat", path, datname);
    output.open(buff, ios::trunc);

    //MyClass t(buff);
    //MyClass* t = new MyClass(buff);
    //TTree *t2 = (TTree *)t.fChain;
    sprintf(buff, "%s/%s.root", path, timename);
    if (gSystem->AccessPathName(buff))
    {
        fexist = Preprocesstimeinformation();
        if (!fexist)
        {

            cout << "Error!! The source File of " << buff << " doesn't exist" << endl;
            return;
        }
    }
    cout << "====>>  Start to open the root file : " << buff << endl;
    TFile *f1 = new TFile(buff, "read");
    TTree *t2 = (TTree *)f1->Get("timeinfo");
    t2->SetMakeClass(1);
    t2->SetBranchAddress("LEDthrd", fLEDth);
    t2->SetBranchAddress("LEDTOT", fLEDTOT);
    t2->SetBranchAddress("timediff", ftimediff);
    t2->SetBranchAddress("width", &fwidth);
    t2->SetBranchAddress("charge", &fcharge);
    int NEvents = t2->GetEntries();
    TCanvas *cCor;
    cCor = cdC(0);
    int CNum = 0;
    TH1D *htime = new TH1D("htime", "Time Resolution;T (ps);Counts", 101e3, tL, tR);
    TH1D *hTOT = new TH1D("hTOT", "TOT;T (ps);Counts", 1e3, UL, UR);

    TH2D *hATorigin[2];
    TH2D *hAT[2];

    for (int h = 0; h < 1; h++)
    {
        sprintf(buff, "hAT%d", h);
        hATorigin[h] = new TH2D(buff, "hAT;TOT (ps); Timediff (ps)", 1e3, UL, UR, 101e3, tL, tR);
    }
    hAT[0] = (TH2D *)hATorigin[0]->Clone();
    TF1 *fitAT = new TF1("fitAT", "0", UL, UR);
    TF1 *fitT;
    TF1 *fittr;
    TF1 *fitbg;
    TF1 *fitTOT;
    TH1F *htfit;
    TH1F *hTOTfit;
    double TOTmean;
    double TOTsigma;
    double TRratio;
    double BGratio;
    //for (int k = start; k < start + 1; k++)
    for (int k = 0; k < ThNum; k++)
    {
        if (k != 0)
            hTOT->Reset();
        for (int s = 0; s < iter; s++)
        {
            //cout << "progress check iter1: " << s << endl;
            if (s != 0)
            {
                if (TOTmean - 3 * TOTsigma <= 0)

                    fitAT = profilefit(hAT[0], rbU, rbt * 4, tL, tR, 0.1, TOTmean + 3 * TOTsigma, buff);
                else
                    fitAT = profilefit(hAT[0], rbU, rbt * 4, tL, tR, TOTmean - 3 * TOTsigma, TOTmean + 3 * TOTsigma, buff);
                if (!fitAT)
                {
                    cout << " the profilefit is failed! " << endl;
                    return;
                }
                htime->Reset();
                hAT[0]->Reset();
            }

            for (int i = 0; i < NEvents; i++)
            {
                t2->GetEntry(i);
                //fLEDth[k] = t.MCP2_LEDthrd[k];
                //fLEDtime[k] = t.MCP2_LEDtime[k];
                //fLEDTOT[k] = t.MCP2_TOT[k];
                //frise = t.MCP2_rise_time[0];
                //fcharge = t.MCP2_all_charge[0];
                //freftime = t.TR1_CFDtime[6];
                //cout<<fLEDth[k]<<"\t"<<fLEDTOT[k]<<"\t"<<fLEDtime[k]<<"\t"<<frise[0]<<"\t"<<fcharge[0]<<"\t"<<freftime[6]<<endl;
                //if(fcharge[0]>5) cout<<fLEDTOT[k]<<endl;
                if (ftimediff[k] > 0)
                {
                    if (s == 0)
                    {
                        hTOT->Fill(fLEDTOT[CorNo]);
                        T[i] = ftimediff[k];
                        //tL = 23.5;
                        //tR = 25;
                        tL = 18;
                        tR = 23;
                        //tL = 43.5;
                        //tR = 46;
                    }

                    //if (s == 0 && fcharge > 0.05) cout << "T[i]=" << T[i] << "\t,freftime" << freftime << endl;
                    else
                    {

                        TAcor = fitAT->Eval(fLEDTOT[CorNo]);
                        T[i] = T[i] - TAcor;
                        tL = -1;
                        tR = 1;
                    }
                    htime->Fill(T[i]);
                    hAT[0]->Fill(fLEDTOT[CorNo], T[i]);

                    //if (s == 1) cout << "TAcor=" << TAcor<< ",T[i]=" << T[i] << "\t,fLEDTOT[k]" << fLEDTOT[k] << endl;
                }
                // cout<<"process check"<<endl;
            }

            //cCor = cdC(CNum++);
            //ht->Draw();

            //cCor = cdC(CNum++);
            //hAT[0]->Draw();

            cCor->cd();
            if (s == 0)
            {
                hTOTfit = (TH1F *)gausfit(hTOT, 0.1, 3, 3, rbU, UL, UR);
                if (!hTOTfit)
                {
                    cout << "the TOTfit hist is NULL " << endl;
                    return;
                }
                DrawMyHist(hTOTfit, "", "", 1, 3);
                fitTOT = (TF1 *)hTOTfit->GetFunction("fitU");
                TOTmean = fitTOT->GetParameter(1);
                TOTsigma = fitTOT->GetParameter(2);
                cout << "TOTmean=" << TOTmean << ",\tTOTsigma=" << TOTsigma << endl;
                sprintf(buff, "%s/%s_TH%.2fmV_TOT.png", path, name, fLEDth[k] * 1e3);
                cCor->SaveAs(buff);
            }
            cCor->Clear();
            //if (s == 1)    return;
            fitfailed = 1;
            counter = 0;
            while (fitfailed)
            {

                htfit = (TH1F *)twogausfit(htime, fac[counter], rangeL, rangeR, rbt, tL, tR);
                //htfit->Draw();
                if (!htfit)
                {
                    cout << "The twogausfit failed!" << endl;
                    sprintf(buff, "%s/%s_TH%.2fmV_TR_cor%d.png", path, name, fLEDth[k] * 1e3, s);
                    cCor->SaveAs(buff);
                    return;
                }
                fitT = (TF1 *)htfit->GetFunction("fit2");
                if (fitT->GetChisquare() / fitT->GetNDF() < 2 || counter >= 5)
                    fitfailed = 0;
                counter++;
            }
            fittr = new TF1("fittr", "gaus", tL, tR);
            fittr->SetParameter(0, fitT->GetParameter(0));
            fittr->SetParameter(1, fitT->GetParameter(1));
            fittr->SetParameter(2, fitT->GetParameter(2));
            fitbg = new TF1("fitbg", "gaus", tL, tR);
            fitbg->SetParameter(0, fitT->GetParameter(3));
            fitbg->SetParameter(1, fitT->GetParameter(4));
            fitbg->SetParameter(2, fitT->GetParameter(5));
            TRratio = fittr->Integral(tL, tR) / (fittr->Integral(tL, tR) + fitbg->Integral(tL, tR));
            BGratio = fitbg->Integral(tL, tR) / (fittr->Integral(tL, tR) + fitbg->Integral(tL, tR));
            double sigma = fitT->GetParameter(2) * 1e3;

            DrawMyHist(htfit, "", "", 1, 3);
            htfit->SetNdivisions(505);
            sprintf(buff, "#sigma=%.0fps,%.2f%%", sigma, TRratio * 100);
            TLatex *l = DrawMyLatex(buff, 0.2, 0.5);
            l->Draw();
            if (!fitT)
            {
                fitT = new TF1("fitT", "landau", tL, tR);
                TFitResultPtr failed = htime->Fit("RQ");
                if (failed)
                {
                    fitT = 0;
                    cout << "gaus and landau fit both failed !" << endl;
                }
            }

            sprintf(buff, "%s/%s_TH%.2fmV_TR_cor%d.png", path, name, fLEDth[k] * 1e3, s);
            cCor->SaveAs(buff);
            if (fitT)
                output << k << "\t" << fLEDth[k] << "\t" << s << "\t" << fitT->GetParameter(2) << "\t" << fitT->GetParError(2) << endl;
            else
                output << fLEDth[k] << "\t" << s << "\t"
                       << "0000"
                       << "\t"
                       << "0000" << endl;
            //sprintf(buff, "%s/%s_Th%.2fmV_At_pfx_cor%d", path, name, fLEDth[k]*1e3, s);
            //if(h==1) {hAT[h]->Draw("colz");return;}

            //cout << "progress check iter2: " << s << endl;
        }
        //cout << "progress check iter3: " << endl;
    }
    //cout << "progress check iter4: " << endl;
    //return;
}
int drawtimeinfo(double left = 35, double right = 55)
{

    setgStyle();

    //**
    //** set your parameters **//
    //

    int N = 14;
    int iterN = 5;
    int THID;
    int iter;
    double fttserr[iterN][N];
    double ftts[iterN][N];
    double ttserr[2][N];
    double tts[2][N];

    double x[N];
    //****************************
    //

    ifstream input;

    char str[1024];
    char buff[1024];
    sprintf(buff, "%s/%s.dat", path, datname);
    //sprintf(buff,"%s%d.dat",name,k);
    input.open(buff);

    if (!input)
    {
        cout << "file can't be found: \n"
             << buff << endl;
        return 0;
    }
    int ct = 0;
    for (int i = 0; i < N; i++)
    {
        /*for(i=0;i<charN;i++)
        {
            input>>str;
            cout<<str<<endl;
        }
        */
        cout << "enter the file: " << buff << endl;
        //for( j=0;j<4;j++){
        for (int k = 0; k < iterN; k++)
        {
            input >> THID >> x[ct] >> iter >> ftts[k][ct] >> fttserr[k][ct];
            //input >> x[ct] >> iter >> ftts[k][ct] >> fttserr[k][ct];
            cout << x[ct] << "\t" << iter << "\t" << ftts[k][ct] << "\t" << fttserr[k][ct] << endl;
            x[ct] = x[ct] * 1e3;
            ftts[k][ct] = TMath::Abs(ftts[k][ct]) * 1e3;
            fttserr[k][ct] = fttserr[k][ct] * 1e3;
            //tr[j][k]=hitsigma*1e3;
        }
        tts[0][ct] = ftts[0][ct];
        ttserr[0][ct] = fttserr[0][ct];
        tts[1][ct] = 9999;
        for (int j = 1; j < iterN - 1; j++)
        {
            if (tts[1][ct] > ftts[j][ct])
            {

                tts[1][ct] = ftts[j][ct];
                ttserr[1][ct] = fttserr[j][ct];
            }
        }

        ct++;
        //}
        //if(j==4&&k==10)
        //break;
    }
    input.close();
    TGraphErrors *g1 = new TGraphErrors(N, x, tts[0], 0, ttserr[0]);
    TGraphErrors *g2 = new TGraphErrors(N, x, tts[1], 0, ttserr[1]);

    TCanvas *c1;
    c1 = cdC(0);

    g1->Draw("AP");

    //mydraw.Graph(g1,"NPE","TimeRes (ps)",1.5,20,4);
    DrawMyGraph(g1, "Fixed Threshold (mV)", "TTS (ps)", 1.5, 20, 1, 1, 2);
    DrawMyGraph(g2, "Fixed Threshold (mV)", "TTS (ps)", 1.5, 20, 2, 2, 2);
    g1->GetYaxis()->SetRangeUser(left, right);
    g1->GetYaxis()->SetNdivisions(505);
    g1->GetXaxis()->SetNdivisions(505);
    g2->Draw("sameP");

    TLegend *leg = DrawMyLeg(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(g1, "TTS", "lp");
    leg->AddEntry(g2, "TTS after TAcorrection", "lp");
    leg->Draw();

    sprintf(buff, "%s/%sTTS.png", path, datname);
    c1->SaveAs(buff);
    return 0;
}
void drawall()
{
    gethist();
    //drawSPE();
    //drawTR();
    //drawrise();
    //drawwidth();
    drawctratio();
    drawringratio();
    drawinvringratio();
    drawavewaveform();
}