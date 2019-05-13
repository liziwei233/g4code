#include "../Include/LZWfunc.h"
#include <iostream>

LZWfunc::LZWfunc()
{
    //*
    //* if calculate average time resolution, set ADD =1;
    ADD = 0;
    name = "defaultname";
    t.L = -1e3;
    t.R = 1e3;
    npk = 10;
}

LZWfunc::~LZWfunc()
{
    cout << "The destructor is called" << endl;
}

// #define __PEAKS_C_FIT_AREAS__ 1 /* fit peaks' areas */
Double_t LZWfunc::fpeaks(Double_t *x, Double_t *par)
{
    //Double_t result = par[0] + par[1]*x[0];
    Double_t result = par[0] * TMath::Gaus(x[0], par[1], par[2]);

    for (Int_t p = 0; p < npk; p++)
    {
        Double_t norm = par[3 * p + 3]; // "height" or "area"
        Double_t mean = par[3 * p + 4];
        Double_t sigma = par[3 * p + 5];
#if defined(__PEAKS_C_FIT_AREAS__)
        norm /= sigma * (TMath::Sqrt(TMath::TwoPi())); // "area"
#endif                                                 /* defined(__PEAKS_C_FIT_AREAS__) */
        result += norm * TMath::Gaus(x[0], mean, sigma);
    }
    return result;
}

TF1 *LZWfunc::fpeaksfit(TH1 *ha, int npeaks, double res, double sigma, double thrd)
{
    npk = TMath::Abs(npeaks);
    TCanvas *cfp1 = new TCanvas("cfp1", "cfp1", 800, 600);

    cfp1->cd();
    TH1 *h = (TH1 *)ha->Clone("h");
    TH1 *h2 = (TH1 *)h->Clone("h2");
    h->Draw();
    RANGE hx;
    hx.L = h->GetXaxis()->GetXmin();
    hx.R = h->GetXaxis()->GetXmax();
    cout << "hx.L" << hx.L << ",   hx.R" << hx.R << endl;
    double par[300] = {0};

    double *xpeaks;

    TSpectrum *s = new TSpectrum(2 * npk, res);
    int nfound = s->Search(h, sigma, "", thrd);
    npk = 0;
    printf("Found %d candidate peaks to fit\n", nfound);
    //ha->Draw();
    //return;
    TH1 *hb = s->Background(h, 20, "same");
    if (hb)
        cfp1->Update();
    //if (np <npeaks) return;

    //estimate linear background using a fitting method
    //cfp->cd(1);
    hb->Draw();
    TF1 *fbk = new TF1("fbk", "gaus", hx.L, hx.R);
    fbk->SetLineColor(3);
    hb->Fit("fbk", "R");
    // Loop on all found peaks. Eliminate peaks at the background level
    par[0] = fbk->GetParameter(0);
    par[1] = fbk->GetParameter(1);
    par[2] = fbk->GetParameter(2);
    //fbk->Draw("same");
    gPad->Clear();
    h->Draw();
    h->Rebin();
    fbk->Draw("same");
    TF1 *fp[nfound];

    xpeaks = s->GetPositionX();

    for (int p = 0; p < nfound; p++)
    {
        Double_t xp = xpeaks[p];
        Int_t bin = h->GetXaxis()->FindBin(xp);
        Double_t yp = h->GetBinContent(bin);
        if (yp - TMath::Sqrt(yp) < fbk->Eval(xp))
            continue;
        sprintf(buff, "fp%d", p);

        fp[p] = new TF1(buff, "gaus", hx.L, hx.R);
        fp[p]->SetLineColor(3);
        h->Fit(fp[p], "+q", "", xp - res * sigma, xp + res * sigma);
        //fp[p]->Draw("same");
        if (fp[p]->GetParameter(1) > hx.L && fp[p]->GetParameter(1) < hx.R && fp[p]->GetParameter(0) > 0)
        {
            par[3 * npk + 3] = fp[p]->GetParameter(0); // "height"
            par[3 * npk + 4] = fp[p]->GetParameter(1); // "mean"
            par[3 * npk + 5] = fp[p]->GetParameter(2); // "sigma"
#if defined(__PEAKS_C_FIT_AREAS__)
            par[3 * npk + 3] *= par[3 * npk + 5] * (TMath::Sqrt(TMath::TwoPi())); // "area"
#endif                                                                            /* defined(__PEAKS_C_FIT_AREAS__) */
            npk++;
        }
    }

    printf("Found %d useful peaks to fit\n", npk);
    sprintf(buff, "%s_1.png", name.c_str());
    cfp1->SaveAs(buff);
    //printf("Now fitting: Be patient\n");
    //return fp[1];
    TCanvas *cfp2 = new TCanvas("cfp2", "cfp2", 800, 600);
    cfp2->cd();
    h2->Draw();
    h2->Rebin();
    cfp2->Update();
    //TF1 *fit = new TF1("fit",fpeaks,0,1000,2+3*npeaks);
    LZWfunc *myfunc = new LZWfunc();
    TF1 *fit = new TF1("fit", myfunc, &LZWfunc::fpeaks, hx.L, hx.R, 3 + 3 * npk);
    // We may have more than the default 25 parameters
    TVirtualFitter::Fitter(h2, 10 + 3 * npk);
    fit->SetParameters(par);
    fit->SetNpx(1000);
    //TPaveStats* stat = (TPaveStats*)h2->FindObject("stats");
    //h->SetStats(0);
    h2->Fit("fit");
    for (int i = 0; i < npk; i++)

    {
        gPar[i].h = fit->GetParameter(3 * i + 3);
        gPar[i].m = fit->GetParameter(3 * i + 4);
        gPar[i].s = fit->GetParameter(3 * i + 5);
    }
    sort(gPar, gPar + npk);

    TPaveStats *stat = (TPaveStats *)h2->GetListOfFunctions()->FindObject("stats");
    stat->SetOptFit(1100);
    sprintf(buff, "%s_2.png", name.c_str());
    cfp2->SaveAs(buff);
    //fit->Draw("same");
    //gPad->Update();
    return fit;
}
double LZWfunc::mcpfun(double *x, double *par){
  Double_t val = 0.;
  Double_t amp = par[0];
  Double_t alpha = par[1];
  Double_t lambda = par[2];
  Double_t ped = par[3];
  Double_t pedSigma = par[4];
  Double_t peak = par[5];
  Double_t  sigma = par[6];
  Double_t bw = par[7];

  //TF1 *myGaus = new TF1("myGaus","[0]*exp(-pow((x-[1])/[2],2)/2.)/sqrt(2.*TMath::Pi()*pow([2],2))",0,4000);
  //myGaus->SetParameters(TMath::PoissonI(0,lambda), ped, pedSigma);
  //val += myGaus->Eval(x[0]);
  
  //
  //** dark noise
  //
  val += alpha*TMath::Gaus(x[0], ped + peak, sigma, kTRUE);
  
  //
  //** pedstal
  //
  val += TMath::PoissonI(0,lambda)*TMath::Gaus(x[0], ped , pedSigma, kTRUE);
  
  
  //
  //** photonelectron
  //
  for(int i=1;i<8;i++) {
    //myGaus->SetParameters(TMath::PoissonI(i,lambda), ped+i*peak, sigma*sqrt(i));// CHANGED!!
    //val += myGaus->Eval(x[0]);
    val += TMath::PoissonI(i,lambda)*TMath::Gaus(x[0], ped + i*peak, sigma*sqrt(i), kTRUE);
  }

  //delete myGaus;
  return amp*val*bw;
}
double LZWfunc::pmtfun(double *x, double *par){
  Double_t val = 0.;
  Double_t amp = par[0];
  Double_t lambda = par[1];
  Double_t ped = par[2];
  Double_t pedSigma = par[3];
  Double_t peak = par[4];
  Double_t  sigma = par[5];
  Double_t bw = par[6];

  //TF1 *myGaus = new TF1("myGaus","[0]*exp(-pow((x-[1])/[2],2)/2.)/sqrt(2.*TMath::Pi()*pow([2],2))",0,4000);
  //myGaus->SetParameters(TMath::PoissonI(0,lambda), ped, pedSigma);
  //val += myGaus->Eval(x[0]);
  val += TMath::PoissonI(0,lambda)*TMath::Gaus(x[0], ped , pedSigma, kTRUE);

  for(int i=1;i<4;i++) {
    //myGaus->SetParameters(TMath::PoissonI(i,lambda), ped+i*peak, sigma*sqrt(i));// CHANGED!!
    //val += myGaus->Eval(x[0]);
    val += TMath::PoissonI(i,lambda)*TMath::Gaus(x[0], ped + i*peak, sigma*sqrt(i), kTRUE);
  }

  //delete myGaus;
  return amp*val*bw;
}
double LZWfunc::HVfun(double *x, double *par){
    double val=0.;
    double A=par[1];
    double alpha=par[2];
    double C=par[0];
    //val = A*TMath::Power(x[0],beta);
    val = TMath::Exp(C+x[0]*A*alpha);
    return val;
}
TF1* LZWfunc::SPSfit(TH1* h,int rbq,RANGE u,double fac)
{
    TH1 *hqdc = (TH1 *)h->Clone();
    hqdc->Draw();
    hqdc->Rebin(rbq);
    hqdc->SetLineColor(1);
    TF1 *myGaus = new TF1("myGaus","gaus",u.L,u.R);
    int ibin =  hqdc->GetMaximumBin();
    double pedMean = hqdc->GetBinCenter(ibin);
    //
    //* try to fit pedstal
    hqdc->Fit(myGaus,"","",u.L, pedMean+pedMean);//???????????
    pedMean = myGaus->GetParameter(1);
    double pedSigma = myGaus->GetParameter(2);
    //
    //* fit pedstal formly
    hqdc->Fit(myGaus,"","",pedMean-10*pedSigma, pedMean+fac*pedSigma);
    pedMean = myGaus->GetParameter(1);
    pedSigma = myGaus->GetParameter(2);
    cout<<" init. par.: pedmean = "<<pedMean<<"; pedsigma = "<<pedSigma<<endl;
    //return myGaus;
    /*
    hqdc->Fit(myGaus,"","",pedMean-10*pedSigma, pedMean+3*pedSigma);
    pedMean = myGaus->GetParameter(1);
    pedSigma = myGaus->GetParameter(2);
    */

    //TF1 *myFun = new TF1("myFun",pmtfun,u.L,u.R,7);
    //if(pedSigma>5e-3) hqdc->GetXaxis()->SetRangeUser(pedMean+12*pedSigma,u.R);
    //else  
    
    //
    //* find the position of SPE
    hqdc->GetXaxis()->SetRangeUser(pedMean+fac*pedSigma,  u.R);
    ibin = hqdc->GetMaximumBin();
    double mean = hqdc->GetBinCenter(ibin)-pedMean;
    double sigma = hqdc->GetStdDev()/10;
    double mean2 = hqdc->GetMean()-pedMean;
    if(mean2/mean>5) mean = mean2;

    myGaus->SetParameter(1,mean);
    myGaus->SetParameter(2,sigma);
    hqdc->Fit(myGaus,"","",pedMean+fac*pedSigma, pedMean+fac*pedSigma+1*mean);
    mean=myGaus->GetParameter(1);
    sigma=myGaus->GetParameter(2);
    cout<<" init. par.: mean = "<<mean<<"; sigma = "<<sigma<<endl;
   //return myGaus;
    
    hqdc->GetXaxis()->SetRangeUser(u.L,  u.R);

    LZWfunc *myfunc = new LZWfunc();
    TF1 *myFun = new TF1("myFun", myfunc, &LZWfunc::pmtfun, u.L, u.R, 7);
    //TF1 *myFun = new TF1("myFun", myfunc, &LZWfunc::mcpfun, u.L, u.R, 8);
    myFun->SetParNames("N","#lambda","#mu_{ped}","#sigma_{ped}","#mu","#sigma", "BW");

     //Int_t ibin =  h->GetMaximumBin();
  myFun->SetParameters(hqdc->GetEntries(),0.1, pedMean, pedSigma, mean, sigma, hqdc->GetBinWidth(5));
  myFun->FixParameter(0, hqdc->GetEntries());//fix total yield
  myFun->FixParameter(6, hqdc->GetBinWidth(5));//fix bin width
  myFun->SetParLimits(2, pedMean-0.0001,pedMean+0.0001);//fix pedestal mean
  myFun->SetParLimits(3, pedSigma-0.000001,pedSigma+0.000001);//fix pedestal sigma
  myFun->SetParLimits(4, pedMean+fac*pedSigma, 2*mean);//>10 for 1400V
  myFun->SetParLimits(5, 0, 2*sigma);
  //if(pedSigma>5) myFun->SetRange(pedMean+12*pedSigma,  u.R);
  //else 
  myFun->SetRange(pedMean+fac*pedSigma,  mean+1*sigma);
  cout<<" first fitting..."<<endl;
  hqdc->Fit(myFun,"R");
  //return myFun;
  mean=myFun->GetParameter(4);//return;
  sigma = myFun->GetParameter(5);
  //myFun->SetRange(pedMean+mean-1.3*TMath::Abs(sigma),  u.R);
  myFun->SetRange(pedMean+mean-1.3*TMath::Abs(sigma),  pedMean+5*mean);
  //myFun->SetParLimits(4,mean/2., u.R);
  //myFun->SetParameter(4,mean);
  //myFun->SetParameter(5,sigma);
  myFun->SetParLimits(4,mean/2., u.R);
  
  //myFun->SetRange(u.L,  u.R);
  //cout<<" second fitting..."<<endl;
  //myFun->SetParLimits(4,mean/2., u.R);
  hqdc->Fit(myFun,"R");
  //return myFun;
  cout<<" fitting done"<<endl;
  mean = myFun->GetParameter(4);
  sigma = myFun->GetParameter(5);
  Double_t par[10], parErr[10];
  for(int j=0;j<7;j++) {
    par[j] = myFun->GetParameter(j);
    parErr[j] = myFun->GetParError(j);
  }

  Double_t xmin = 0, xmax = 0;
  if((pedMean-15*pedSigma)>0) xmin = pedMean-12*pedSigma;
  if(par[1]>2) hqdc->GetXaxis()->SetRangeUser(xmin,pedMean+10*mean);
  else if(par[1]>1) hqdc->GetXaxis()->SetRangeUser(xmin,pedMean+7*mean);
  else if(par[1]>0.1) hqdc->GetXaxis()->SetRangeUser(xmin,pedMean+5*mean);
  else if(par[1]>0.05) hqdc->GetXaxis()->SetRangeUser(xmin,pedMean+4*mean);
  else hqdc->GetXaxis()->SetRangeUser(xmin,pedMean+4*mean);

  hqdc->SetLineWidth(1.5);
  myGaus->SetLineWidth(2);
  myGaus->SetLineColor(4);
  myGaus->Draw("same");

  myFun->SetNpx(1000);
  //myFun->SetLineWidth(1);
  myFun->SetLineColor(1);
  //myFun->SetLineStyle(7);
  myFun->GetRange(xmin, xmax);
  myFun->SetRange(0, xmax);
  myFun->DrawCopy("same");
  myFun->SetLineStyle(1);
  myFun->SetLineWidth(3);
  myFun->SetRange(pedMean-4*pedSigma,u.R);
  //myFun->Draw("same");

  TF1 *peakFun[10];
  for(int j=1;j<4;j++) {
    sprintf(buff,"peakFun%d",j);
    peakFun[j] = new TF1(buff,"[0]*exp(-pow((x-[1])/[2],2)/2.)/sqrt(2.*TMath::Pi()*pow([2],2))",u.L,u.R);
    //peakFun[j] = new TF1(buff,"gaus",u.L,u.R);
    //if(j==0) peakFun[j]->SetParameters(par[1]*par[0]*par[7],par[3]+par[5], par[6]);
    //else 
    peakFun[j]->SetParameters(TMath::PoissonI(j, par[1])*par[0]*par[6], par[2]+j*par[4], par[5]*sqrt(1.0*j));//CHANGED!!
    //if(j==0) peakFun[j]->SetLineColor(6);
    //else 
    peakFun[j]->SetLineColor(6);
    peakFun[j]->SetLineWidth(2);
    peakFun[j]->Draw("same");
  }

  myFun->Draw("same");
  return myFun;
}
TF1* LZWfunc::mcpSPfit(TH1* h,int rbq,RANGE u,double fac)
{
    TH1 *hqdc = (TH1 *)h->Clone();
    hqdc->Draw();
    hqdc->Rebin(rbq);
    hqdc->SetLineColor(1);
    TF1 *pedGaus = new TF1("pedGaus","gaus",u.L,u.R);
    int ibin =  hqdc->GetMaximumBin();
    double pedMean = hqdc->GetBinCenter(ibin);
    //
    //* try to fit pedstal
    hqdc->Fit(pedGaus,"","",0-pedMean, pedMean+1*pedMean);//???????????
    pedMean = pedGaus->GetParameter(1);
    double pedSigma = pedGaus->GetParameter(2);
    //
    //* fit pedstal formly
    hqdc->Fit(pedGaus,"","",pedMean-10*pedSigma, pedMean+1*pedSigma);
    pedMean = pedGaus->GetParameter(1);
    pedSigma = pedGaus->GetParameter(2);
    //return pedGaus;
    /*
    hqdc->Fit(myGaus,"","",pedMean-10*pedSigma, pedMean+3*pedSigma);
    pedMean = myGaus->GetParameter(1);
    pedSigma = myGaus->GetParameter(2);
    */

    //TF1 *myFun = new TF1("myFun",pmtfun,u.L,u.R,7);
    //if(pedSigma>5e-3) hqdc->GetXaxis()->SetRangeUser(pedMean+12*pedSigma,u.R);
    //else  
    
    TF1 *SPEGaus = new TF1("SPEGaus","gaus",u.L,u.R);
    //
    //* find the position of SPE
    hqdc->GetXaxis()->SetRangeUser(pedMean+fac*pedSigma,  u.R);
    ibin = hqdc->GetMaximumBin();
    double mean = hqdc->GetBinCenter(ibin)-pedMean;
    double sigma = hqdc->GetStdDev()/10;
    double mean2 = hqdc->GetMean()-pedMean;
    if(mean2/mean>5) mean = mean2;

    SPEGaus->SetParameter(1,mean);
    SPEGaus->SetParameter(2,sigma);
    hqdc->Fit(SPEGaus,"","",pedMean+fac*pedSigma, pedMean+fac*pedSigma+1*mean);
    cout<<" fit range.: L = "<<pedMean+fac*pedSigma<<"; R = "<<pedMean+fac*pedSigma+2*mean<<endl;
    mean=SPEGaus->GetParameter(1);
    sigma=SPEGaus->GetParameter(2);
    cout<<" init. par.: mean = "<<mean<<"; sigma = "<<sigma<<endl;
    //return SPEGaus;
  
    hqdc->GetXaxis()->SetRangeUser(pedMean-2*fac*pedSigma,  pedMean+5*mean);
    TF1*  myGaus = new TF1("myGaus", "gaus(0)+gaus(3)", u.L, u.R);
    myGaus->SetParNames("C_{ped}", "#mu_{ped}", "#sigma_{ped}", "C_{PE}", "#mu_{PE}", "#sigma_{PE}");
    //myGaus->SetParameter(0, pedGaus->GetParameter(0));
    myGaus->FixParameter(0, pedGaus->GetParameter(0));
    myGaus->FixParameter(1, pedMean);
    myGaus->FixParameter(2, pedSigma);
    //myGaus->SetParameter(3, SPEGaus->GetParameter(0));
    myGaus->FixParameter(3, SPEGaus->GetParameter(0));
    myGaus->FixParameter(4, mean);
    myGaus->FixParameter(5, sigma);
    hqdc->Fit(myGaus);
   //return SPEGaus;
  

  
  
  hqdc->SetLineWidth(1.5);
  SPEGaus->SetLineWidth(2);
  SPEGaus->SetLineColor(2);
  SPEGaus->Draw("same");
  pedGaus->SetLineWidth(2);
  pedGaus->SetLineColor(4);
  pedGaus->Draw("same");

  myGaus->SetLineWidth(2);
  myGaus->SetLineColor(1);
  myGaus->Draw("same");

  
  return myGaus;
}

TF1 *LZWfunc::gausfit(TH1 *h, int rbU, double fac, RANGE U)
{

    double mean = 0;
    double sigma = 0;
    double max=0;
    TH1 *hU = (TH1 *)h->Clone();
    hU->Draw();
    hU->Rebin(rbU);
    hU->GetXaxis()->SetRangeUser(U.L, U.R);
    TF1 *fitU = new TF1("fitU", "gaus", U.L, U.R);
    max = hU->GetBinCenter(hU->GetMaximumBin());
    fitU->SetParLimits(0,0, hU->GetBinContent(hU->GetMaximumBin()));
    fitU->SetParameter(1, max);
    //cout << mean << "\t" << sigma << endl;
    hU->Fit(fitU, "R");
    mean = fitU->GetParameter(1);
    sigma = fitU->GetParameter(2);

    cout << mean << "\t" << sigma << endl;

    TFitResultPtr failed = hU->Fit(fitU, "", "", mean - fac * sigma, mean + fac * sigma);
    //failed =1 means fit failed
    if (failed)
        return fitU = 0;

    else
    {

        if (U.L < mean - 8 * sigma)
            U.L = mean - 8 * sigma;
        if (U.R > mean + 8 * sigma)
            U.R = mean + 8 * sigma;

        hU->GetXaxis()->SetRangeUser(U.L, U.R);

        return fitU;
    }
}

TF1 *LZWfunc::gausfit(TH1 *h, int rbU, double fac, RANGE *U)
{

    double mean = 0;
    double sigma = 0;
    TH1 *hU = (TH1 *)h->Clone();
    hU->Draw();
    hU->Rebin(rbU);
    hU->GetXaxis()->SetRangeUser((*U).L, (*U).R);
    TF1 *fitU = new TF1("fitU", "gaus", (*U).L, (*U).R);
    mean = hU->GetBinCenter(hU->GetMaximumBin());
    fitU->SetParameter(1, mean);
    cout << mean << "\t" << sigma << endl;
    hU->Fit(fitU, "R");
    mean = fitU->GetParameter(1);
    sigma = fitU->GetParameter(2);

    cout << mean << "\t" << sigma << endl;

    TFitResultPtr failed = hU->Fit(fitU, "", "", mean - fac * sigma, mean + fac * sigma);
    //failed =1 means fit failed
    if (failed)
        return fitU = 0;

    else
    {

        if ((*U).L < mean - 8 * sigma)
            (*U).L = mean - 8 * sigma;
        if ((*U).R > mean + 8 * sigma)
            (*U).R = mean + 8 * sigma;

        hU->GetXaxis()->SetRangeUser((*U).L, (*U).R);

        return fitU;
    }
}

TF1 *LZWfunc::twoguasfit(TH1 *ht, double fac, int rbt, RANGE *t)
{
    //First fit for ensuring the rangement of histgram;
    TH1 *h = (TH1 *)ht->Clone();
    h->Rebin(rbt);
    double mean = h->GetBinCenter(h->GetMaximumBin());
    double sigma = h->GetRMS();
    TF1 *fit = new TF1("fit", "gaus", (*t).L, (*t).R);
    h->GetXaxis()->SetRangeUser((*t).L, (*t).R);
    fit->SetParameter(1, mean);
    //fit->SetParameter(2,sigma);
    h->Fit(fit);
    mean = fit->GetParameter(1);
    sigma = TMath::Abs(fit->GetParameter(2));
    if ((*t).L < mean - 5 * sigma || sigma > 1)
    {
        (*t).L = mean - 5 * sigma;
        (*t).R = mean + 5 * sigma;
    }
    cout << h->GetName() << "\t" << (*t).L << "\t" << (*t).R << endl;

    h->GetXaxis()->SetRangeUser((*t).L, (*t).R);

    TF1 *fit2 = new TF1("fit2", "gaus(0)+gaus(3)", (*t).L, (*t).R);
    fit2->SetParNames("C_{TR}", "#mu_{TR}", "#sigma_{TR}", "C_{bkgnd}", "#mu_{bkgnd}", "#sigma_{bkgnd}");
    fit2->SetParameter(1, mean);
    fit2->SetParameter(2, sigma);
    fit2->SetParLimits(3, 0, fit->GetParameter(0) * fac);
    fit2->SetParameter(4, mean);
    fit2->SetParameter(5, 1.5 * sigma);

    h->Fit(fit2, "", "", mean - 5 * sigma, mean + 5 * sigma);
    TF1 *fit_tr = new TF1("fit_tr", "gaus", (*t).L, (*t).R);
    fit_tr->SetParameter(0, fit2->GetParameter(0));
    fit_tr->SetParameter(1, fit2->GetParameter(1));
    fit_tr->SetParameter(2, fit2->GetParameter(2));
    TF1 *fit_bg = new TF1("fit_bg", "gaus", (*t).L, (*t).R);
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

    if ((*t).L < mean - 10 * sigma)
    {
        (*t).L = mean - 10 * sigma;
    }
    if ((*t).R > mean + 10 * sigma)
    {
        (*t).R = mean + 10 * sigma;
    }

    cout << h->GetName() << "\t" << (*t).L << "\t" << (*t).R << endl;
    h->GetXaxis()->SetRangeUser((*t).L, (*t).R);
    return fit2;
}
TF1 *LZWfunc::twoguasfit(TH1 *ht, double fac, int rbt, RANGE t)
{
    //First fit for ensuring the rangement of histgram;
    TH1 *h = (TH1 *)ht->Clone();
    h->Rebin(rbt);
    double mean = h->GetBinCenter(h->GetMaximumBin());
    double sigma = h->GetRMS();
    TF1 *fit = new TF1("fit", "gaus", t.L, t.R);
    h->GetXaxis()->SetRangeUser(t.L, t.R);
    fit->SetParameter(1, mean);
    //fit->SetParameter(2,sigma);
    h->Fit(fit);
    mean = fit->GetParameter(1);
    sigma = TMath::Abs(fit->GetParameter(2));
    if (t.L < mean - 5 * sigma || sigma > 1)
    {
        t.L = mean - 5 * sigma;
        t.R = mean + 5 * sigma;
    }
    cout << h->GetName() << "\t" << t.L << "\t" << t.R << endl;

    h->GetXaxis()->SetRangeUser(t.L, t.R);

    TF1 *fit2 = new TF1("fit2", "gaus(0)+gaus(3)", t.L, t.R);
    fit2->SetParNames("C_{TR}", "#mu_{TR}", "#sigma_{TR}", "C_{bkgnd}", "#mu_{bkgnd}", "#sigma_{bkgnd}");
    fit2->SetParameter(1, mean);
    fit2->SetParameter(2, sigma);
    fit2->SetParLimits(3, 0, fit->GetParameter(0) * fac);
    fit2->SetParameter(4, mean);
    fit2->SetParameter(5, 1.5 * sigma);

    h->Fit(fit2, "", "", mean - 5 * sigma, mean + 5 * sigma);
    TF1 *fit_tr = new TF1("fit_tr", "gaus", t.L, t.R);
    fit_tr->SetParameter(0, fit2->GetParameter(0));
    fit_tr->SetParameter(1, fit2->GetParameter(1));
    fit_tr->SetParameter(2, fit2->GetParameter(2));
    TF1 *fit_bg = new TF1("fit_bg", "gaus", t.L, t.R);
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

    if (t.L < mean - 10 * sigma)
    {
        t.L = mean - 10 * sigma;
    }
    if (t.R > mean + 10 * sigma)
    {
        t.R = mean + 10 * sigma;
    }

    cout << h->GetName() << "\t" << t.L << "\t" << t.R << endl;
    h->GetXaxis()->SetRangeUser(t.L, t.R);
    return fit2;
}

TF1 *LZWfunc::profilefit(TH2 *Rt, double rbU, double rbt, RANGE t, RANGE U, char *name)
{

    TCanvas *cpfx = new TCanvas("cpfx", "cpfx", 1600, 600);
    cpfx->Clear();
    cpfx->Divide(2, 1);
    cpfx->cd(1);
    TH2 *Qt = (TH2 *)Rt->Clone();
    //TH2* Qt = (TH2*) Rt->Clone("tmp");

    Qt->Draw("colz");

    Qt->RebinX(rbU);
    Qt->RebinY(rbt);
    Qt->GetYaxis()->SetRangeUser(t.L, t.R);
    Qt->GetXaxis()->SetRangeUser(U.L, U.R);
    Qt->ProfileX();

    cpfx->cd(2);
    TH1 *Qpfx = Qt->ProfileX();
    //Qpfx->Reset();

    //Qpfx=(TH1*)gDirectory->Get("Qt_pfx");
    Qpfx->Draw();
    Qpfx->GetYaxis()->SetRangeUser(t.L, t.R);
    Qpfx->GetXaxis()->SetRangeUser(U.L, U.R);

    TF1 *fitQt = new TF1("fitQt", "[0]+[1]/TMath::Sqrt(abs(x))+[2]/abs(x)+[3]/abs(x)/TMath::Sqrt(abs(x))+[4]/abs(x)/abs(x)", U.L, U.R);

    Qpfx->Fit(fitQt, "R");

    sprintf(buff, "%s.png", name);
    cpfx->SaveAs(buff);
    Qpfx->Reset();
    //delete Qt;
    //delete Qpfx;
    //delete c5;
    return fitQt;
}

void LZWfunc::CH3Correction(TTree *t1, vector<EVENT *> ch, double *p, vector<charRANGE> range, vector<CUT> cut, OPTION opt, string name)
{
    //*
    //** The number of channels : 3
    //** Discrimination: multi
    //*

    int rbU = opt.rbu;
    int rbt = opt.rbt;
    int iter = opt.iter;
    double fac = opt.fac;
    RANGE t = range.at(0).t;
    vector<RANGE> U;
    U.push_back(range.at(0).y);
    U.push_back(range.at(1).y);
    U.push_back(range.at(2).y);
    //(A+B)/2-MCP
    //write for beamtest data of SiPM prototype
    //with complicated cut

    //Correction QT
    ofstream output;
    sprintf(buff, "%s.dat", name.c_str());
    output.open(buff, ios::trunc);

    TCanvas *cCor = new TCanvas("cCor", "cCor", 800, 600);
    TH1D *ht = new TH1D("ht", "Time Resolution;T (ps);Counts", 2e3, t.L, t.R);

    TF1 *fitT;

    RANGE initial_t;
    initial_t.L = t.L;
    initial_t.R = t.R;

    double TAcor = 0;
    double T[1000000] = {0};
    const int chN = 3;
    double* Q[chN];
    for(int i=0;i<chN;i++){
#ifndef G4_FLAG
                    Q[i]=&(*ch.at(i)).Q;
#else
                    Q[i]=&(*ch.at(i)).Amp;
#endif
    }
    TH2D *hAT[chN];
    TF1 *fitAT = new TF1("fitAT", "0", U.at(0).L, U.at(0).R);

    for (int h = 0; h < chN; h++)
    {
        sprintf(buff, "hAT%d", h);
        hAT[h] = new TH2D(buff, buff, 200, U.at(h).L, U.at(h).R, 4e3, t.L, t.R);
    }

    double L = 20;
    for (int j = 0; j < L; j++)
    {
        //multi CFD fraction dicrimination
        for (int s = 0; s < iter; s++)
        {

            int N = t1->GetEntries();
            cout << "Total entries is :" << N << endl;

            //correction of A1
            for (int h = 0; h < chN; h++)
            {
                t.L = initial_t.L;
                t.R = initial_t.R;

                for (int i = 0; i < N; i++)
                {
                    t1->GetEntry(i);


                    if (h - 1 < 0)
                        TAcor = fitAT->Eval(*Q[chN - 1]);
                    else
                        TAcor = fitAT->Eval(*Q[h-1]);
                    if (abs(p[j]) < 1e-12)
                        return fitAT;
                    else
                    {
                        int il = 0; // if all boolean conditions are true , il==chN.
                        for (il = 0; il < chN; il++)
                        {
                            if (!ifstat(*ch.at(il), cut.at(il)))
                                break;
                        }
                        if (il == chN)
                        {

                            if (s == 0 && h == 0)
                                T[i] = (*ch.at(2)).CFD[j] - ((*ch.at(0)).CFD[j] + (*ch.at(1)).CFD[j]) / 2;

                            else
                                T[i] = T[i] - TAcor;

                            ht->Fill(T[i]);
                            hAT[h]->Fill(*Q[h], T[i]);
                        }
                    }
                }

                cCor->cd();
                cCor->Clear();
                ht->Draw();
                //return;
                //t.L=-5;
                //t.R=5;

                fitT = gausfit(ht, rbt, fac, &t);
                if (!fitT)
                {
                    fitT = new TF1("fitT", "landau", t.L, t.R);
                    TFitResultPtr failed = ht->Fit("R");
                    if (failed)
                    {
                        fitT = 0;
                        cout << "gaus and landau fit both failed !" << endl;
                    }
                }
                sprintf(buff, "%s_ch%d_CFDfrac%g_TR_cor%d.png", name.c_str(), h, p[j], s);
                cCor->SaveAs(buff);
                if (!fitT)
                    output << name.c_str() << "\t" << h << "\t" << p[j] << "\t" << s << "\t" << fitT->GetParameter(2) << "\t" << fitT->GetParError(2) << endl;
                else
                    output << name.c_str() << "\t" << h << "\t" << p[j] << "\t" << s << "\t"
                           << "0000"
                           << "\t"
                           << "0000" << endl;

                sprintf(buff, "%s_ch%d_CFDfrac%gAt_pfx_cor%d", name.c_str(), h, p[j], s);

                fitAT = profilefit(hAT[h], rbU, rbt * 8, t, U.at(h), buff);

                ht->Reset();
                hAT[h]->Reset();
            }
        }

        fitAT = new TF1("fitAT", "0", U.at(0).L, U.at(0).R);
    }
    output.close();
}

TF1 *LZWfunc::CH3Correction(TTree *t1, vector<EVENT *> ch, vector<charRANGE> range, vector<CUT> cut, OPTION opt, string name)
{
    //*
    //** The number of channels : 3
    //** Discrimination: single
    //*

    int rbU = opt.rbu;
    int rbt = opt.rbt;
    int iter = opt.iter;
    double fac = opt.fac;
    RANGE t = range.at(0).t;
    vector<RANGE> U;
    U.push_back(range.at(0).y);
    U.push_back(range.at(1).y);
    U.push_back(range.at(2).y);
    //(A+B)/2-MCP
    //write for beamtest data of SiPM prototype
    //with complicated cut

    //Correction QT
    ofstream output;
    sprintf(buff, "%s.dat", name.c_str());
    output.open(buff, ios::trunc);

    TCanvas *cCor = new TCanvas("cCor", "cCor", 800, 600);
    TH1D *ht = new TH1D("ht", "Time Resolution;T (ps);Counts", 2e3, t.L, t.R);

    TF1 *fitT;

    RANGE initial_t;
    initial_t.L = t.L;
    initial_t.R = t.R;

    double TAcor = 0;
    double T[1000000] = {0};
    const int chN = 3;
    double* Q[chN];
    for(int i=0;i<chN;i++){
#ifndef G4_FLAG
                    Q[i]=&(*ch.at(i)).Q;
#else
                    Q[i]=&(*ch.at(i)).Amp;
#endif
    }

    TH2D *hAT[chN];
    TF1 *fitAT = new TF1("fitAT", "0", U.at(0).L, U.at(0).R);

    for (int h = 0; h < chN; h++)
    {
        sprintf(buff, "hAT%d", h);
        hAT[h] = new TH2D(buff, buff, 200, U.at(h).L, U.at(h).R, 4e3, t.L, t.R);
    }

    for (int s = 0; s < iter; s++)
    {

        int N = t1->GetEntries();
        cout << "Total entries is :" << N << endl;

        //correction of A1
        for (int h = 0; h < chN; h++)
        {
            t.L = initial_t.L;
            t.R = initial_t.R;

            for (int i = 0; i < N; i++)
            {
                t1->GetEntry(i);

                if (h - 1 < 0)
                    TAcor = fitAT->Eval(*Q[chN-1]);
                else
                    TAcor = fitAT->Eval(*Q[h-1]);

                int il = 0; // if all boolean conditions are true , il==chN.
                for (il = 0; il < chN; il++)
                {
                    if (!ifstat(*ch.at(il), cut.at(il)))
                        break;
                }
                if (il == chN)
                {

                    if (s == 0 && h == 0)
                        T[i] = (*ch.at(2)).time - ((*ch.at(0)).time + (*ch.at(1)).time) / 2;

                    else
                        T[i] = T[i] - TAcor;

                    ht->Fill(T[i]);
                    hAT[h]->Fill(*Q[h], T[i]);
                }
            }

            cCor->cd();
            cCor->Clear();
            ht->Draw();
            //return;
            //t.L=-5;
            //t.R=5;

            fitT = gausfit(ht, rbt, fac, &t);
            if (!fitT)
            {
                fitT = new TF1("fitT", "landau", t.L, t.R);
                TFitResultPtr failed = ht->Fit("R");
                if (failed)
                {
                    fitT = 0;
                    cout << "gaus and landau fit both failed !" << endl;
                }
            }
            sprintf(buff, "%s_ch%d_TR_cor%d.png", name.c_str(), h, s);
            cCor->SaveAs(buff);
            if (!fitT)
                output << name.c_str() << "\t" << h << "\t" << s << "\t" << fitT->GetParameter(2) << "\t" << fitT->GetParError(2) << endl;
            else
                output << name.c_str() << "\t" << h << "\t" << s << "\t"
                       << "0000"
                       << "\t"
                       << "0000" << endl;

            sprintf(buff, "%s_ch%d_At_pfx_cor%d", name.c_str(), h, s);

            fitAT = profilefit(hAT[h], rbU, rbt * 8, t, U.at(h), buff);

            ht->Reset();
            hAT[h]->Reset();
        }
    }

    fitAT = new TF1("fitAT", "0", U.at(0).L, U.at(0).R);

    output.close();
    return fitT;
}

TF1 *LZWfunc::CH2Correction(TTree *t1, vector<EVENT *> ch, vector<charRANGE> range, vector<CUT> cut, OPTION opt, string name)
{
    //*
    //** The number of channels : 2
    //** Discrimination: single
    //*

    //pass your parameter
    RANGE t = range.at(0).t;
    vector<RANGE> U;
    U.push_back(range.at(0).y);
    U.push_back(range.at(1).y);

    RANGE initial_t;
    initial_t.L = t.L;
    initial_t.R = t.R;

    int iter = opt.iter;
    double fac = opt.fac;
    int rbt = opt.rbt;
    int rbU = opt.rbu;

    ofstream output;
    sprintf(buff, "%s.dat", name.c_str());
    output.open(buff, ios::trunc);
    //cout<<buff<<" has been built"<<endl;

    TCanvas *cCor = new TCanvas("cCor", "cCor", 800, 600);
    TH1D *ht = new TH1D("ht", "Time Resolution;T (ps);Counts", 1e3, t.L, t.R);

    TF1 *fitT;

    double TAcor = 0;

    double T[1000000] = {0};

    //EVENT** ch[chN];
    const int chN = 2;
    double* Q[chN];
    for(int i=0;i<chN;i++){
#ifndef G4_FLAG
                    //Q[i]=&(*ch.at(i)).Q;
                    Q[i]=&(*ch.at(i)).charge[0];
#else
                    Q[i]=&(*ch.at(i)).Amp;
#endif
    }
    TH2D *hAT[chN];
    TF1 *fitAT = new TF1("fitAT", "0", U.at(0).L, U.at(0).R);
    for (int h = 0; h < chN; h++)
    {
        sprintf(buff, "hAT%d", h);
        hAT[h] = new TH2D(buff, buff, 200, U.at(h).L, U.at(h).R, 4e3, t.L, t.R);
    }

    //ch[0]=&A;
    //ch[1]=&B;
    //ch[2]=&MCP;
    //double L=20;
    //cout<<"the number of values of fraction of CFD: "<<L<<endl;
    //for(int j=0;j<L;j++){
    //multi CFD fraction dicrimination

    for (int s = 0; s < iter; s++)
    {

        int N = t1->GetEntries();
        cout << "Total entries is :" << N << endl;

        //correction of A1
        for (int h = 0; h < chN; h++)
        {
            t.L = initial_t.L;
            t.R = initial_t.R;
            for (int i = 0; i < N; i++)
            {
                t1->GetEntry(i);

                //Q1=(*A).Q;
                //Q2=(*B).Q;
                if (h - 1 < 0)
                    TAcor = fitAT->Eval(*Q[chN-1]);
                else
                    TAcor = fitAT->Eval(*Q[h-1]);
                int il = 0; // if all boolean conditions are true , il==chN.
                for (il = 0; il < chN; il++)
                {
                    if (!ifstat(*ch.at(il), cut.at(il)))
                        break;
                }
                if (il == chN)
                {

                    if (s == 0 && h == 0)
                    {
                      if(ADD){
                        // average the time of all channels ;
                        for (int iT = 0; iT < chN; iT++)
                        {
                            T[i] += (*ch.at(iT)).time;
                        }
                        T[i] = T[i] / chN;
                        }
                      else
                        T[i] = (*ch.at(0)).time - (*ch.at(1)).time;

                    }
                    else
                        T[i] = T[i] - TAcor;

                    ht->Fill(T[i]);
                    hAT[h]->Fill(*Q[h], T[i]);
                }
            }

            cCor->cd();
            cCor->Clear();
            ht->Draw();
            //return;
            //tRL=-5;
            //tRR=5;

            //twoguasfit(ht,&tL,&tR,fac,rbt);
            fitT = gausfit(ht, rbt, fac, &t);
            //	sprintf(buff,"%s_TR_ch%d_cor%d.png",name,h,s);
            sprintf(buff, "%s_ch%d_TR_cor%d.png", name.c_str(), h, s);
            cCor->SaveAs(buff);

            //sprintf(buff,"%s_At_ch%d_pfx_cor%d",name,h,s);

            //*
            //* fit your histgram by gaus-function or landau function.
            fitT = gausfit(ht, rbt, fac, &t);
            if (!fitT)
            {
                fitT = new TF1("fitT", "landau", t.L, t.R);
                TFitResultPtr failed = ht->Fit("R");
                if (failed)
                {
                    fitT = 0;
                    cout << "gaus and landau fit both failed !" << endl;
                }
            }
            sprintf(buff, "%s_ch%d_TR_cor%d.png", name.c_str(), h, s);
            cCor->SaveAs(buff);

            //*
            //* output your fit parameter if fit successed.
            if (!fitT)
                output << name.c_str() << "\t" << h << "\t" << s << "\t" << fitT->GetParameter(2) << "\t" << fitT->GetParError(2) << endl;
            else
                output << name.c_str() << "\t" << h << "\t" << s << "\t"
                       << "0000"
                       << "\t"
                       << "0000" << endl;

            sprintf(buff, "%s_ch%d_At_pfx_cor%d", name.c_str(), h, s);
            fitAT = profilefit(hAT[h], rbU, rbt * 8, t, U.at(h), buff);

            ht->Reset();
            hAT[h]->Reset();
        }
    }
    fitAT = new TF1("fitAT", "0", U.at(0).L, U.at(0).R);

    //}

    output.close();
    return fitT;
}

void LZWfunc::CH2Correction(TTree *t1, vector<EVENT *> ch, double *p, vector<charRANGE> range, vector<CUT> cut, OPTION opt, string name)
{
    //*
    //** The number of channels : 2
    //** Discrimination: multi
    //*

    int iter = opt.iter;
    double fac = opt.fac;
    int rbt = opt.rbt;
    int rbU = opt.rbu;
    RANGE t = range.at(0).t;
    vector<RANGE> U;
    U.push_back(range.at(0).q);
    U.push_back(range.at(1).q);

    RANGE initial_t;
    initial_t.L = t.L;
    initial_t.R = t.R;

    ofstream output;
    sprintf(buff, "%s.dat", name.c_str());
    output.open(buff, ios::trunc);
    //cout<<buff<<" has been built"<<endl;

    TCanvas *cCor = new TCanvas("cCor", "cCor", 800, 600);
    TH1D *ht = new TH1D("ht", "Time Resolution;T (ps);Counts", 1e3, t.L, t.R);

    TF1 *fitT;

    double TAcor = 0;

    double T[1000000] = {0};
    const int chN = 2;
    double* Q[chN];
    for(int i=0;i<chN;i++){
#ifndef G4_FLAG
                    //Q[i]=&(*ch.at(i)).Q;
                    Q[i]=&(*ch.at(i)).charge[0];
#else
                    Q[i]=&(*ch.at(i)).Amp;
#endif
    }
    //EVENT** ch[chN];
    TH2D *hAT[chN];
    TF1 *fitAT = new TF1("fitAT", "0", U.at(0).L, U.at(0).R);
    for (int h = 0; h < chN; h++)
    {
        sprintf(buff, "hAT%d", h);
        hAT[h] = new TH2D(buff, buff, 200, U.at(h).L, U.at(h).R, 4e3, t.L, t.R);
    }

    //ch[0]=&A;
    //ch[1]=&B;
    //ch[2]=&MCP;
    double L = 20;
    //cout<<"the number of values of fraction of CFD: "<<L<<endl;
    for (int j = 0; j < L; j++)
    {
        //multi CFD fraction dicrimination

        for (int s = 0; s < iter; s++)
        {

            int N = t1->GetEntries();
            cout << "Total entries is :" << N << endl;

            //correction of A1
            for (int h = 0; h < chN; h++)
            {
                t.L = initial_t.L;
                t.R = initial_t.R;
                for (int i = 0; i < N; i++)
                {
                    t1->GetEntry(i);

                    //Q1=(*A).Q;
                    //Q2=(*B).Q;
                    if (h - 1 < 0)
                        TAcor = fitAT->Eval(*Q[chN-1]);
                    else
                        TAcor = fitAT->Eval(*Q[h-1]);

                    //TAcor=fitAT->Eval(Q2);
                    //if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*A).rise<4&&(*B).rise<4){
                    //if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*A).rise<4&&(*B).rise<4&&(*MCP).Q>0.3){
                    //if(1){
                    //if(1){
                    //if((*A).A>0.3&&(*B).A>0.3){
                    if (abs(p[j]) < 1e-12)
                    {
                        //cout<<"There are "<<j<<" thredholds !"<<endl;
                        return fitAT;
                    }
                    else
                    {
                        int il = 0; // if all boolean conditions are true , il==chN.
                        for (il = 0; il < chN; il++)
                        {
                            if (!ifstat(*ch.at(il), cut.at(il)))
                                break;
                        }
                        //cout<<"There are "<<il<<" conditions have been meeted !"<<endl;

                        if (il == chN)
                        {

                            if (s == 0 && h == 0)
                            {
                             if(ADD){
                                // average the time of all channels ;
                                for (int iT = 0; iT < chN; iT++)
                                {
                                    T[i] += (*ch.at(iT)).CFD[j];
                                }
                                T[i] = T[i] / chN;
                                }
                             else
                                T[i] = (*ch.at(0)).CFD[j] - (*ch.at(1)).CFD[j];

                            }

                            else
                                T[i] = T[i] - TAcor;
                            //cout<<T[i]<<"\t"<<*Q[h]<<endl;
                            ht->Fill(T[i]);
                            hAT[h]->Fill(*Q[h], T[i]);
                        }
                    }
                }

                cCor->cd();
                cCor->Clear();
                ht->Draw();
                //return;
                //tRL=-5;
                //tRR=5;

                //twoguasfit(ht,&tL,&tR,fac,rbt);
                fitT = gausfit(ht, rbt, fac, &t);
                //	sprintf(buff,"%s_TR_ch%d_cor%d.png",name,h,s);
                if (!fitT)
                {
                    fitT = new TF1("fitT", "landau", t.L, t.R);
                    TFitResultPtr failed = ht->Fit("R");
                    if (failed)
                    {
                        fitT = 0;
                        cout << "gaus and landau fit both failed !" << endl;
                    }
                }
                sprintf(buff, "%s_ch%d_CFDfrac%g_TR_cor%d.png", name.c_str(), h, p[j], s);
                cCor->SaveAs(buff);
                if (!fitT)
                    output << name.c_str() << "\t" << h << "\t" << p[j] << "\t" << s << "\t" << fitT->GetParameter(2) << "\t" << fitT->GetParError(2) << endl;
                else
                    output << name.c_str() << "\t" << h << "\t" << p[j] << "\t" << s << "\t"
                           << "0000"
                           << "\t"
                           << "0000" << endl;

                sprintf(buff, "%s_ch%d_CFDfrac%gAt_pfx_cor%d", name.c_str(), h, p[j], s);
                //sprintf(buff,"%s_At_ch%d_pfx_cor%d",name,h,s);

                fitAT = profilefit(hAT[h], rbU, rbt * 8, t, U.at(h), buff);

                ht->Reset();
                hAT[h]->Reset();
            }
        }
        fitAT = new TF1("fitAT", "0", U.at(0).L, U.at(0).R);
    }

    output.close();
}

TF1 *LZWfunc::CH1Correction(TTree *t1, EVENT *A, charRANGE range, CUT cut, OPTION opt, string name)
{
    //*
    //** The number of channels : 1
    //** Discrimination: signle
    //*

    int iter = opt.iter;
    double fac = opt.fac;
    int rbt = opt.rbt;
    int rbU = opt.rbu;

    RANGE t = range.t;
    RANGE U = range.y;

    double Q;
    RANGE initial_t;
    initial_t.L = -5;
    initial_t.R = 5;

    ofstream output;
    sprintf(buff, "%s.dat", name.c_str());
    output.open(buff, ios::trunc);
    TCanvas *cCor = new TCanvas("cCor", "cCor", 800, 600);

    TH1D *ht;
    TH2D *hAT;

    TH1D *htraw = new TH1D("htraw", "time Resolution;T (ps);Counts", 2e4, t.L, t.R);

    TH2D *hATraw = new TH2D("hATraw", "", 600, U.L, U.R, 2e4, t.L, t.R);

    TH1D *htcor = new TH1D("htcor", "time Resolution;T (ps);Counts", 2e4, initial_t.L, initial_t.R);

    TH2D *hATcor = new TH2D("hATcor", "", 600, U.L, U.R, 2e4, initial_t.L, initial_t.R);

    TF1 *fitAT = new TF1("fitAT", "0", U.L, U.R);
    TF1 *fitT;

    double TAcor = 0;
    double T[1000000] = {0};
    for (int s = 0; s < iter; s++)
    {

        int N = t1->GetEntries();
        cout << "Total entries is :" << N << endl;
        for (int i = 0; i < N; i++)
        {
            t1->GetEntry(i);
#ifndef G4_FLAG
            Q = (*A).Q;
#else 
            Q = (*A).Amp;
#endif
            TAcor = fitAT->Eval(Q);
            //if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3){
            //if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3&&(*A).rise<4&&(*B).rise<4){
            // else if((*A).A<0&&(*A).time>0){
            if (ifstat((*A), cut))
            {
                if (s == 0)
                //T[i]=(*A).time-(*MCP).time-TAcor;
                {
                    T[i] = (*A).time;
                    htraw->Fill(T[i]);
                    hATraw->Fill(Q, T[i]);
                }
                else
                {
                    T[i] = T[i] - TAcor;
                    htcor->Fill(T[i]);
                    hATcor->Fill(Q, T[i]);
                }
                //				cout<<"Q="<<Q<<" T="<<T[i]<<endl;
            }
        }

        //cout<<"<====process check====>"<<endl;

        cCor->cd();
        cCor->Clear();
        if (s == 0)
        {
            ht = (TH1D *)htraw->Clone();
            hAT = (TH2D *)hATraw->Clone();
        }
        else
        {
            ht = (TH1D *)htcor->Clone();
            hAT = (TH2D *)hATcor->Clone();
            cout << "progress check" << endl;
        }
        //return;
        //t.L=-5;
        //t.R=5;
        t.L = ht->GetMean() - 5 * ht->GetRMS();
        t.R = ht->GetMean() + 5 * ht->GetRMS();
        //cout<<"iter="<<s<<", MEAN="<<ht->GetMean()<<"\t RMS="<<ht->GetRMS()<<endl;
        //cout<<"iter="<<s<<", tL="<<tL<<"\t tR="<<tR<<endl;
        //if (s==1) return;
        //twoguasfit(ht,&tL,&tR,fac,rbt);
        fitT = gausfit(ht, rbt, fac, &t);
        if (!fitT)
        {
            fitT = new TF1("fitT", "landau", t.L, t.R);
            TFitResultPtr failed = ht->Fit("R");
            if (failed)
            {
                fitT = 0;
                cout << "gaus and landau fit both failed !" << endl;
            }
        }
        sprintf(buff, "%s_TR_cor%d.png", name.c_str(), s);
        cCor->SaveAs(buff);

        if (!fitT)
            output << name.c_str() << "\t" << s << "\t" << fitT->GetParameter(2) << "\t" << fitT->GetParError(2) << endl;
        else
            output << name.c_str() << "\t" << s << "\t"
                   << "0000"
                   << "\t"
                   << "0000" << endl;

        sprintf(buff, "%s_At_pfx_cor%d", name.c_str(), s);

        fitAT = profilefit(hAT, rbU, rbt * 2, t, U, buff);

        ht->Reset();
        hAT->Reset();
        htraw->Reset();
        hATraw->Reset();
        htcor->Reset();
        hATcor->Reset();
    }

    output.close();
    return fitT;
}

TF1 *LZWfunc::CH1Correction(TTree *t1, EVENT *A, double *t0, charRANGE range, CUT cut, OPTION opt, string name)
{
    //*
    //** The number of channels : 1 and reffered t0
    //** Discrimination: signle
    //*
    int iter = opt.iter;
    double fac = opt.fac;
    int rbt = opt.rbt;
    int rbU = opt.rbu;

    RANGE t = range.t;
    RANGE U = range.y;

    double Q;
    RANGE initial_t;
    initial_t.L = -5;
    initial_t.R = 5;

    ofstream output;
    sprintf(buff, "%s.dat", name.c_str());
    output.open(buff, ios::trunc);
    TCanvas *cCor = new TCanvas("cCor", "cCor", 800, 600);

    TH1D *ht;
    TH2D *hAT;

    TH1D *htraw = new TH1D("htraw", "time Resolution;T (ps);Counts", 2e4, t.L, t.R);

    TH2D *hATraw = new TH2D("hATraw", "", 600, U.L, U.R, 2e4, t.L, t.R);

    TH1D *htcor = new TH1D("htcor", "time Resolution;T (ps);Counts", 2e4, initial_t.L, initial_t.R);

    TH2D *hATcor = new TH2D("hATcor", "", 600, U.L, U.R, 2e4, initial_t.L, initial_t.R);

    TF1 *fitAT = new TF1("fitAT", "0", U.L, U.R);
    TF1 *fitT;

    double TAcor = 0;
    double T[1000000] = {0};
    for (int s = 0; s < iter; s++)
    {

        int N = t1->GetEntries();
        cout << "Total entries is :" << N << endl;
        for (int i = 0; i < N; i++)
        {
            t1->GetEntry(i);
#ifndef G4_FLAG
            Q = (*A).Q;
#else
            Q = (*A).Amp;
#endif 
            TAcor = fitAT->Eval(Q);
            //if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3){
            //if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3&&(*A).rise<4&&(*B).rise<4){
            // else if((*A).A<0&&(*A).time>0){
            if (ifstat((*A), cut))
            {
                if (s == 0)
                //T[i]=(*A).time-(*MCP).time-TAcor;
                {
                    T[i] = (*A).time - *t0;
                    htraw->Fill(T[i]);
                    hATraw->Fill(Q, T[i]);
                }
                else
                {
                    T[i] = T[i] - TAcor;
                    htcor->Fill(T[i]);
                    hATcor->Fill(Q, T[i]);
                }
                //				cout<<"Q="<<Q<<" T="<<T[i]<<endl;
            }
        }

        //cout<<"<====process check====>"<<endl;

        cCor->cd();
        cCor->Clear();
        if (s == 0)
        {
            ht = (TH1D *)htraw->Clone();
            hAT = (TH2D *)hATraw->Clone();
        }
        else
        {
            ht = (TH1D *)htcor->Clone();
            hAT = (TH2D *)hATcor->Clone();
            cout << "progress check" << endl;
        }
        //return;
        //t.L=-5;
        //t.R=5;
        t.L = ht->GetMean() - 5 * ht->GetRMS();
        t.R = ht->GetMean() + 5 * ht->GetRMS();
        //cout<<"iter="<<s<<", MEAN="<<ht->GetMean()<<"\t RMS="<<ht->GetRMS()<<endl;
        //cout<<"iter="<<s<<", tL="<<tL<<"\t tR="<<tR<<endl;
        //if (s==1) return;
        //twoguasfit(ht,&tL,&tR,fac,rbt);
        fitT = gausfit(ht, rbt, fac, &t);
        if (!fitT)
        {
            fitT = new TF1("fitT", "landau", t.L, t.R);
            TFitResultPtr failed = ht->Fit("R");
            if (failed)
            {
                fitT = 0;
                cout << "gaus and landau fit both failed !" << endl;
            }
        }
        sprintf(buff, "%s_TR_cor%d.png", name.c_str(), s);
        cCor->SaveAs(buff);

        if (!fitT)
            output << name.c_str() << "\t" << s << "\t" << fitT->GetParameter(2) << "\t" << fitT->GetParError(2) << endl;
        else
            output << name.c_str() << "\t" << s << "\t"
                   << "0000"
                   << "\t"
                   << "0000" << endl;

        sprintf(buff, "%s_At_pfx_cor%d", name.c_str(), s);

        fitAT = profilefit(hAT, rbU, rbt * 2, t, U, buff);

        ht->Reset();
        hAT->Reset();
        htraw->Reset();
        hATraw->Reset();
        htcor->Reset();
        hATcor->Reset();
    }

    output.close();
    return fitT;
}

void LZWfunc::drawcharacter(TTree *t3, int chN, string *chname, vector<charRANGE> chR)
{
    DrawMyfunc draw;
#ifndef G4_FLAG
    for (int i = 0; i < chN; i++)
    {
        TGaxis::SetMaxDigits(3);
        TH2F *hqy = new TH2F("hqy", "", 200, chR.at(i).y.L, chR.at(i).y.R, 200, chR.at(i).q.L, chR.at(i).q.R);
        draw.Hist(hqy, "Amp (V)", "Charge (pC)", 1.5);
        sprintf(buff, "MCP%d_all_charge:MCP%d_global_maximum_y>>hqy", i + 1, i + 1);
        t3->Draw(buff, "", "colz");
        sprintf(buff, "%s_Qvsy.png", chname[i].data());
        gPad->SaveAs(buff);

        TH2F *hqr = new TH2F("hqr", "", 200, chR.at(i).r.L, chR.at(i).r.R, 200, chR.at(i).q.L, chR.at(i).q.R);
        draw.Hist(hqr, "risetime (ns)", "Charge (pC)", 1.5);
        sprintf(buff, "MCP%d_all_charge:MCP%d_rise_time>>hqr", i + 1, i + 1);
        t3->Draw(buff, "", "colz");
        sprintf(buff, "%s_Qvsrise.png", chname[i].data());
        gPad->SaveAs(buff);

        TH2F *hqx = new TH2F("hqx", "", 200, chR.at(i).x.L, chR.at(i).x.R, 200, chR.at(i).q.L, chR.at(i).q.R);
        draw.Hist(hqx, "time (ns)", "Charge (pC)", 1.5);
        sprintf(buff, "MCP%d_all_charge:MCP%d_global_maximum_x>>hqx", i + 1, i + 1);
        t3->Draw(buff, "", "colz");
        sprintf(buff, "%s_Qvsx.png", chname[i].data());
        gPad->SaveAs(buff);

        TH1F *hbl = new TH1F("hbl", "", 200, chR.at(i).bl.L, chR.at(i).bl.R);
        draw.Hist(hbl, "baseline level(V)", "Counts", 1.5);
        sprintf(buff, "MCP%d_baseline_level>>hbl", i + 1);
        t3->Draw(buff, "", "colz");
        sprintf(buff, "%s_bl.png", chname[i].data());
        gPad->SaveAs(buff);

        TH1F *hblrms = new TH1F("hblrms", "", 200, chR.at(i).blrms.L, chR.at(i).blrms.R);
        draw.Hist(hblrms, "baseline RMS(V)", "Counts", 1.5);
        sprintf(buff, "MCP%d_baseline_rms>>hblrms", i + 1);
        t3->Draw(buff, "", "colz");
        sprintf(buff, "%s_blrms.png", chname[i].data());
        gPad->SaveAs(buff);
    }
#endif
    /*TH2F* hqy = new TH2F("hqy","",2e3,chR.at(i).y.L,chR.at(i).y.R,200,chR.at(i).q.L,chR.at(i).q.R);
	TH2F* hqr = new TH2F("hqr","",2e3,,10);
	TH1F* hrise = new TH1F("hy","",2e3,-10,10);
	TH1F* hy = new TH1F("hy","",2e3,-10,10);
	TH1F* hy = new TH1F("hy","",2e3,-10,10);
	TH1F* hy = new TH1F("hy","",2e3,-10,10);*/
}

bool LZWfunc::ifstat(EVENT ch, CUT cut)
{
#ifndef G4_FLAG
/*cout<<ch.Amp <<"\t"<<cut.y.L<<"\t"<< cut.y.R <<endl;
cout<<ch.charge[0] <<"\t"<<cut.q.L<<"\t"<< cut.q.R <<endl;
cout<<ch.rise <<"\t" <<cut.r.L <<"\t"<< cut.r.R <<endl;
cout<<ch.bl <<"\t"<<cut.bl.L <<"\t"<< cut.bl.R <<endl;
cout<<ch.blrms <<"\t" <<cut.blrms.L <<"\t"<< cut.blrms.R <<endl;
cout<<ch.x <<"\t" <<cut.x.L <<"\t"<< cut.x.R <<endl;
*/
    if (ch.Amp > cut.y.L &&
        ch.Amp < cut.y.R &&
        ch.charge[0] > cut.q.L &&
        ch.charge[0] < cut.q.R &&
        ch.rise > cut.r.L &&
        ch.rise < cut.r.R &&
        ch.bl > cut.bl.L &&
        ch.bl < cut.bl.R &&
        ch.blrms > cut.blrms.L &&
        ch.blrms < cut.blrms.R &&
        ch.x > cut.x.L &&
        ch.x < cut.x.R)
        return 1;
#else
    if (ch.Amp > cut.y.L &&
        ch.Amp < cut.y.R &&
        ch.time > cut.t.L &&
        ch.time < cut.t.R &&
        ch.npe > cut.npe.L &&
        ch.npe < cut.npe.R)
        return 1;
#endif
    else
        return 0;
}

string LZWfunc::getTime()
{
    time_t timep;
    time(&timep);
    char tmp[64];
    strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S", localtime(&timep));
    return tmp;
    /*
    how to use
    string time = getTime();
    output<< time << endl;
*/
}