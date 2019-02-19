#include "/mnt/c/Subsys/work/g4code/mycode/Include/LZWfunc.h"
#include <iostream>
using namespace std;

LZWfunc::LZWfunc()
{
    name="defaultname";
    t.L=-1e3;
    t.R=1e3;

    
}

LZWfunc::~LZWfunc()
{
    cout << "The destructor is called" << endl;
}

TF1* LZWfunc::gausfit(TH1 *h, int rbU, double fac, RANGE U)
{

    double mean = 0;
    double sigma = 0;
    TH1 *hU = (TH1 *)h->Clone();

    hU->Draw();
    hU->Rebin(rbU);
    TF1 *fitU = new TF1("fitU", "gaus", U.L, U.R);
    mean = hU->GetBinCenter(hU->GetMaximumBin());
    fitU->SetParameter(1, mean);
    hU->Fit(fitU, "R");
    mean = fitU->GetParameter(1);
    sigma = fitU->GetParameter(2);
    if (U.L < mean - 8 * sigma)
        U.L = mean - 8 * sigma;
    if (U.R > mean + 8 * sigma)
        U.R = mean + 8 * sigma;

    hU->GetXaxis()->SetRangeUser(U.L, U.R);
    hU->Fit(fitU, "", "", mean - fac * sigma, mean + fac * sigma);

    return fitU;
}
TF1* LZWfunc::gausfit(TH1 *h, int rbU, double fac, RANGE* U)
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
    cout<<mean<<"\t"<<sigma<<endl;
    hU->Fit(fitU, "R");
    mean = fitU->GetParameter(1);
    sigma = fitU->GetParameter(2);
    
    cout<<mean<<"\t"<<sigma<<endl;
    
    hU->Fit(fitU, "", "", mean - fac * sigma, mean + fac * sigma);
   
    if ((*U).L < mean - 8 * sigma)
        (*U).L = mean - 8 * sigma;
    if ((*U).R > mean + 8 * sigma)
        (*U).R = mean + 8 * sigma;

    hU->GetXaxis()->SetRangeUser((*U).L, (*U).R);

    return fitU;
}

TF1* LZWfunc::twoguasfit(TH1 *ht, double fac, int rbt, RANGE* t)
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
TF1* LZWfunc::twoguasfit(TH1 *ht, double fac, int rbt, RANGE t)
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

TF1* LZWfunc::profilefit(TH2 *Rt, double rbU, double rbt, RANGE t, RANGE U, char *name)
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

void LZWfunc::CH3Correction(TTree *t1, vector<EVENT *> ch, double* p, int rbU, int rbt, double fac, int iter)
{

    //(A+B)/2-MCP
    //write for beamtest data of SiPM prototype 
    //with complicated cut 

    //Correction QT
    sprintf(buff,"%s_ch3Correction.dat",name.c_str());
    output.open(buff,ios::trunc);

    TCanvas *cCor = new TCanvas("cCor", "cCor", 800, 600);
    TH1D *ht = new TH1D("ht", "Time Resolution;T (ps);Counts", 2e3, t.L, t.R);

	TF1 *fitT;

    
    RANGE initial_t;
    initial_t.L=t.L;
    initial_t.R=t.R;

    double TAcor = 0;
    double T[1000000] = {0};
    const int chN = 3;

    TH2D *hAT[chN];
    TF1 *fitAT= new TF1("fitAT", "0", U.at(0).L, U.at(0).R);

    for (int h=0; h<chN; h++ )
    {
        sprintf(buff,"hAT%d",h);
        hAT[h] = new TH2D(buff,buff, 200, U.at(h).L, U.at(h).R, 4e3, t.L, t.R);
    }     


    double L=20;
    for(int j=0;j<L;j++){
		//multi CFD fraction dicrimination
    for (int s = 0; s < iter; s++)
    {

        int N = t1->GetEntries();
        cout << "Total entries is :" << N << endl;

        //correction of A1
        for (int h = 0; h < chN; h++)
        {
            t.L=initial_t.L;
            t.R=initial_t.R;

            for (int i = 0; i < N; i++)
            {
                t1->GetEntry(i);

                if (h - 1 < 0)
				TAcor=fitAT->Eval((*ch.at(chN-1)).Q);
                else
				TAcor=fitAT->Eval((*ch.at(h-1)).Q);
                
                if(abs(p[j])<1e-12) return fitAT;
                else if(
                    (*ch.at(0)).Amp>cut.at(0).Amplow&&
                    (*ch.at(0)).Amp<cut.at(0).Ampup&&
                    (*ch.at(0)).Q>cut.at(0).Qlow&&
                    (*ch.at(0)).Q<cut.at(0).Qup&&
                    (*ch.at(0)).rise>cut.at(0).riselow&&
                    (*ch.at(0)).rise<cut.at(0).riseup&&
                    (*ch.at(1)).Amp>cut.at(1).Amplow&&
                    (*ch.at(1)).Amp<cut.at(1).Ampup&&
                    (*ch.at(1)).Q>cut.at(1).Qlow&&
                    (*ch.at(1)).Q<cut.at(1).Qup&&
                    (*ch.at(1)).rise>cut.at(1).riselow&&
                    (*ch.at(1)).rise<cut.at(1).riseup&&
                    (*ch.at(2)).Amp>cut.at(2).Amplow&&
                    (*ch.at(2)).Amp<cut.at(2).Ampup&&
                    (*ch.at(2)).Q>cut.at(2).Qlow&&
                    (*ch.at(2)).Q<cut.at(2).Qup&&
                    (*ch.at(2)).rise>cut.at(2).riselow&&
                    (*ch.at(2)).rise<cut.at(2).riseup
                ){
     
				if (s==0&&h==0) 
					T[i]=(*ch.at(2)).CFD[j]-((*ch.at(0)).CFD[j]+(*ch.at(1)).CFD[j])/2;
                            
                        else
                            T[i] = T[i] - TAcor;

                        ht->Fill(T[i]);
				        hAT[h]->Fill((*ch.at(h)).Q,T[i]);
                    
                }
            }

            cCor->cd();
            cCor->Clear();
            ht->Draw();
            //return;
            //t.L=-5;
            //t.R=5;

            fitT = gausfit(ht,rbt,fac,&t);
        sprintf(buff,"%s_ch%d_CFDfrac%g_TR_cor%d.png",name.c_str(),h,p[j],s);
            cCor->SaveAs(buff);
        output << name.c_str() << "\t" << h << "\t" <<p[j] << "\t" << s  << "\t" << fitT->GetParameter(2)  << "\t" << fitT->GetParError(2) << endl;

        sprintf(buff,"%s_ch%d_CFDfrac%gAt_pfx_cor%d",name.c_str(),h,p[j],s);

		fitAT=profilefit(hAT[h],rbU,rbt*8,t,U.at(h),buff);

            ht->Reset();
            hAT[h]->Reset();
        }
    }

fitAT= new TF1("fitAT", "0", U.at(0).L, U.at(0).R);
    }
        output.close();
}

TF1* LZWfunc::CH3Correction(TTree *t1, EVENT *A, EVENT *B, EVENT *MCP, CUT basecut, CUT selcut, int rbU, int rbt, char *name, double fac, int iter, RANGE t, RANGE U)
{

    //(A+B)/2-MCP
    //write for beamtest data of SiPM prototype 
    //with complicated cut 

    //Correction QT

    TCanvas *cCor = new TCanvas("cCor", "cCor", 800, 600);
    TH1D *ht = new TH1D("ht", "Time Resolution;T (ps);Counts", 2e3, t.L, t.R);

    TH2D *hAT = new TH2D("hAT", "", 200, U.L, U.R, 2e3, t.L, t.R);

    TF1 *fitAT = new TF1("fitAT", "0", U.L, U.R);
    RANGE initial_t;
    initial_t.L=t.L;
    initial_t.R=t.R;

    double TAcor = 0;
    double T[1000000] = {0};
    const int chN = 2;
    EVENT **ch[chN];

    ch[0] = &A;
    ch[1] = &B;
    //ch[2]=&MCP;
    for (int s = 0; s < iter; s++)
    {

        int N = t1->GetEntries();
        cout << "Total entries is :" << N << endl;

        //correction of A1
        for (int h = 0; h < chN; h++)
        {
            t.L=initial_t.L;
            t.R=initial_t.R;

            for (int i = 0; i < N; i++)
            {
                t1->GetEntry(i);

                if (h - 1 < 0)
                    TAcor = fitAT->Eval((**ch[chN - 1]).Q);
                else
                    TAcor = fitAT->Eval((**ch[h - 1]).Q);

                //TAcor=fitAT->Eval(Q2);
                if ((*A).Q > basecut.Qlow 
                && (*A).Amp < basecut.Ampup 
                && (*B).Q > basecut.Qlow 
                && (*B).Amp < basecut.Ampup 
                && (*A).rise < basecut.riseup 
                && (*B).rise < basecut.riseup
                && (*MCP).Q > 0.3 )
                {
                    //if(1){
                    if ((*A).Q < selcut.Qup && (*B).Q < selcut.Qup)
                    {
                        if (s == 0)
                            T[i] = ((*A).time + (*B).time) / 2 - (*MCP).time - TAcor;
                        else
                            T[i] = T[i] - TAcor;

                        ht->Fill(T[i]);
                        hAT->Fill((**ch[h]).Q, T[i]);
                    }
                }
            }

            cCor->cd();
            cCor->Clear();
            ht->Draw();
            //return;
            //t.L=-5;
            //t.R=5;

            twoguasfit(ht, fac, rbt, &t);
            sprintf(buff, "%s_TR_ch%d_cor%d.png", name, h, s);
            cCor->SaveAs(buff);

            sprintf(buff, "%s_At_ch%d_pfx_cor%d", name, h, s);

            fitAT = profilefit(hAT, rbU, rbt * 2, t, U, buff);

            ht->Reset();
            hAT->Reset();
        }
    }

    return fitAT;
}

TF1* LZWfunc::CH2Correction(TTree *t1, EVENT *A, EVENT *B, EVENT *MCP, CUT basecut, CUT selcut, int rbU, int rbt, char *name, double fac, int iter, RANGE t, RANGE U)
{
    //A-B
    //write for beamtest data of SiPM prototype 
    //with complicated cut

    //Correction QT

    TCanvas *cCor = new TCanvas("cCor", "cCor", 800, 600);
    TH1D *ht = new TH1D("ht", "Time Resolution;T (ps);Counts", 2e3, -10, 10);

    TH2D *hAT = new TH2D("hAT", "", 200, 0, 1.6e3, 2e3, -10, 10);

    TF1 *fitAT = new TF1("fitAT", "0", -1e4, 1e4);
    
    RANGE initial_t;
    initial_t.L=t.L;
    initial_t.R=t.R;
    U.L = 0;
    U.R = 1.6e3;

    double TAcor = 0;
    double T[1000000] = {0};
    const int chN = 2;
    EVENT **ch[chN];

    ch[0] = &A;
    ch[1] = &B;
    //ch[2]=&MCP;
    for (int s = 0; s < iter; s++)
    {

        int N = t1->GetEntries();
        cout << "Total entries is :" << N << endl;

        //correction of A1
        for (int h = 0; h < chN; h++)
        {
            t.L=initial_t.L;
            t.R=initial_t.R;
            

            for (int i = 0; i < N; i++)
            {
                t1->GetEntry(i);

                if (h - 1 < 0)
                    TAcor = fitAT->Eval((**ch[chN - 1]).Q);
                else
                    TAcor = fitAT->Eval((**ch[h - 1]).Q);

                //TAcor=fitAT->Eval(Q2);
                //if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*A).rise<4&&(*B).rise<4){
                if ((*A).Q > basecut.Qlow 
                && (*A).Amp < basecut.Ampup 
                && (*B).Q > basecut.Qlow 
                && (*B).Amp < basecut.Ampup 
                && (*A).rise < basecut.riseup 
                && (*B).rise < basecut.riseup
                && (*MCP).Q > 0.3 )
                {
                    //if(1){
                    //if(1){
                    if ((*A).Amp > selcut.Amplow 
                    && (*B).Amp > selcut.Amplow)
                    {
                        if (s == 0 && h == 0)
                            //T[i]=(*A).time-(*B).time-TAcor;
                            T[i] = (*A).time - (*B).time;
                        else
                            T[i] = T[i] - TAcor;

                        ht->Fill(T[i]);
                        hAT->Fill((**ch[h]).Q, T[i]);
                    }
                }
            }

            cCor->cd();
            cCor->Clear();
            ht->Draw();
            //return;
            //t.L=-5;
            //t.R=5;

            twoguasfit(ht,fac, rbt, &t);
            sprintf(buff, "%s_TR_ch%d_cor%d.png", name, h, s);
            cCor->SaveAs(buff);

            sprintf(buff, "%s_At_ch%d_pfx_cor%d", name, h, s);

            fitAT = profilefit(hAT, rbU, rbt * 2, t, U, buff);

            ht->Reset();
            hAT->Reset();
        }
    }

    return fitAT;
}

TF1* LZWfunc::CH2Correction(TTree *t1, EVENT *A, EVENT *B, double Uth, int rbU, int rbt, char *name, double fac, int iter, RANGE t, RANGE U)
{
    //A-B 
    //threshold loop

    //Correction QT
    RANGE initial_t;
    initial_t.L=t.L;
    initial_t.R=t.R;

    U.L = 0;
    U.R = 20;

    TCanvas *cCor = new TCanvas("cCor", "cCor", 800, 600);
    TH1D *ht = new TH1D("ht", "Time Resolution;T (ps);Counts", 10e3, t.L, t.R);

    TH2D *hAT = new TH2D("hAT", "", 200, U.L, U.R, 4e3, t.L, t.R);

    TF1 *fitAT = new TF1("fitAT", "0", U.L, U.R);
    
    TF1 *fitT;


    double TAcor = 0;
    double T[1000000] = {0};
    const int chN = 2;
    EVENT **ch[chN];

    ch[0] = &A;
    ch[1] = &B;
    //ch[2]=&MCP;
    for (int s = 0; s < iter; s++)
    {

        int N = t1->GetEntries();
        cout << "Total entries is :" << N << endl;

        //correction of A1
        for (int h = 0; h < chN; h++)
        {
            t.L=initial_t.L;
            t.R=initial_t.R;

            for (int i = 0; i < N; i++)
            {
                t1->GetEntry(i);

                if (h - 1 < 0)
                    TAcor = fitAT->Eval((**ch[chN - 1]).Q);
                else
                    TAcor = fitAT->Eval((**ch[h - 1]).Q);

                //TAcor=fitAT->Eval(Q2);
                //if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*A).rise<4&&(*B).rise<4){
                //if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*A).rise<4&&(*B).rise<4&&(*MCP).Q>0.3){
                if ((*A).time > 0 && (*B).time > 0)
                {
                    //if(1){
                    if ((*A).Q > Uth * 0.32 && (*B).Q > Uth * 0.16)
                    {
                        if (s == 0 && h == 0)
                            //T[i]=(*A).time-(*B).time-TAcor;
                            T[i] = (*A).time - (*B).time;
                        else
                            T[i] = T[i] - TAcor;

                        ht->Fill(T[i]);
                        hAT->Fill((**ch[h]).Q, T[i]);
                    }
                }
            }

            cCor->cd();
            cCor->Clear();
            ht->Draw();
            //return;
            //t.L=-5;
            //t.R=5;
            //cout<<"before tL vs tR= "<<tL<<"\t"<<tR<<endl;
            //twoguasfit(ht,&tL,&tR,fac,rbt);
            fitT = gausfit(ht, rbt, fac, &t);
            //cout<<"gausfit tL vs tR= "<<tL<<"\t"<<tR<<endl;
            sprintf(buff, "%s_Uth%g_TR_ch%d_cor%d.png", name,Uth, h, s);
            cCor->SaveAs(buff);
            output << name << "\t" << Uth << "\t" << h << "\t" << s << "\t" << fitT->GetParameter(2) / sqrt(2) << "\t" << fitT->GetParError(2) << endl;

            sprintf(buff, "%s_Uth%g_At_ch%d_pfx_cor%d", name, Uth, h, s);

            fitAT = profilefit(hAT, rbU, rbt * 2, t, U, buff);
            //cout<<"profile tL vs tR= "<<tL<<"\t"<<tR<<endl;

            ht->Reset();
            hAT->Reset();
        }
    }

    return fitAT;
}

TF1* LZWfunc::CH2Correction(TTree * t1, EVENT * A, EVENT * B, POSITION * pos, double range, int rbU, int rbt, char *name, double fac, int iter, RANGE t,RANGE U)
{
    //A-B
    //npe threshold loop
    //or position selection
    // write for SIMU

        //Correction QT
    RANGE initial_t;
    initial_t.L=t.L;
    initial_t.R=t.R;

        U.L = -3e3;
        U.R = 0e3;

        TCanvas *cCor = new TCanvas("cCor", "cCor", 800, 600);
        //sprintf(buff,"%s_ht",name);
        TH1D *ht = new TH1D("ht", "time Resolution;T (ps);Counts", 1e3, t.L, t.R);
        TH1D *hnpeA = new TH1D("hnpeA", ";npe;Counts", 2e2, 1, 2e2);
        TH1D *hnpeB = new TH1D("hnpeB", ";npe;Counts", 2e2, 1, 2e2);

        //sprintf(buff,"%s_hAT",name);
        TH2D *hAT = new TH2D("hAT", "", 600, U.L, U.R, 1e3, t.L, t.R);

        TF1 *fitAT = new TF1("fitAT", "0", U.L, U.R);
        TF1 *fitT;

        double TAcor = 0;
        const int chN = 2;
        EVENT **ch[chN];

        ch[0] = &A;
        ch[1] = &B;
        //ch[2]=&MCP;
        //for(int j=0;j<12;j++){

        double T[1000000] = {0};

        for (int s = 0; s < iter; s++)
        {

            int N = t1->GetEntries();
            cout << "Total entries is :" << N << endl;

            //correction of A1
            for (int h = 0; h < chN; h++)
            {
            t.L=initial_t.L;
            t.R=initial_t.R;

                for (int i = 0; i < N; i++)
                {
                    t1->GetEntry(i);

                    //Q1=(*A).Q;
                    //Q2=(*B).Q;
                    if (h - 1 < 0)
                        TAcor = fitAT->Eval((**ch[chN - 1]).Amp);
                    else
                        TAcor = fitAT->Eval((**ch[h - 1]).Amp);

                    //TAcor=fitAT->Eval(Q2);
                    //if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*A).rise<4&&(*B).rise<4){
                    //if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*A).rise<4&&(*B).rise<4&&(*MCP).Q>0.3){
                    //if(1){
                    //if(1){
                    //if((*A).A>0.3&&(*B).A>0.3){
                    //if(abs(p[j])<1e-3) return fitAT;
                    //else if((*A).A<-300&&(*A).time>0&&(*B).A<-300&&(*B).time>0){
                    //else if((*A).A<0&&(*A).time>0&&(*B).A<0&&(*B).time>0
                    if ((*A).Amp < 0 && (*A).time > 0 && (*B).Amp < 0 && (*B).time > 0
                        //&&sqrt((*pos).x*(*pos).x+(*pos).y*(*pos).y)>range-1
                        //&&sqrt((*pos).x*(*pos).x+(*pos).y*(*pos).y)<range
                        && (*A).npe > range && (*B).npe > range)
                    {

                        if (s == 0 && h == 0)
                        //T[i]=(*A).time-(*B).time-TAcor;
                        {
                            T[i] = (*A).time - (*B).time;
                            hnpeA->Fill((*A).npe);
                            hnpeB->Fill((*B).npe);
                        }
                        else
                            T[i] = T[i] - TAcor;

                        ht->Fill(T[i]);
                        hAT->Fill((**ch[h]).Amp, T[i]);
                    }
                }

                cCor->cd();
                cCor->Clear();
                ht->Draw();
                //return;
                //t.L=-5;
                //t.R=5;

                //twoguasfit(ht,&tL,&tR,fac,rbt);
                fitT = gausfit(ht, rbt, fac, &t);
                cout << "tL=" << t.L << ",tR=" << t.R << endl;
                //sprintf(buff,"%s_TR_ch%d_cor%d.png",name,h,s);
                sprintf(buff, "%s_TR_range%g_ch%d_cor%d.png", name, range, h, s);
                cCor->SaveAs(buff);
                output << name << "\t" << h << "\t" << s << "\t" << range << "\t" << fitT->GetParameter(2) / sqrt(2) << "\t" << fitT->GetParError(2) << "\t" << hnpeA->GetMean() << "\t" << hnpeB->GetMean() << endl;
                sprintf(buff, "%s_range%g_At_ch%d_pfx_cor%d", name, range, h, s);
                //sprintf(buff,"%s_At_ch%d_pfx_cor%d",name,h,s);

                fitAT = profilefit(hAT, rbU, rbt * 2, t, U, buff);

                ht->Reset();
                hAT->Reset();
            }
        }
        //fitAT=new TF1("fitAT","0",QL,QR);
        //}
        hnpeA->Delete();
        hnpeB->Delete();
        ht->Delete();
        hAT->Delete();

        return fitAT;
    }


//TF1* LZWfunc::CH2Correction(TTree* t1,EVENT* A,EVENT* B,double* p,int rbU,int rbt,char* name,double fac,int iter, RANGE t, RANGE U){
void LZWfunc::CH2Correction(TTree* t1,vector<EVENT*> ch, double* p,int rbU,int rbt,double fac,int iter)
{
	//Correction QT

	RANGE initial_t;
    initial_t.L=t.L;
    initial_t.R=t.R;


    sprintf(buff,"%s_ch2Correction.dat",name.c_str());
    output.open(buff,ios::trunc);
    //cout<<buff<<" has been built"<<endl;

	TCanvas *cCor = new TCanvas("cCor", "cCor", 800, 600);
    TH1D *ht = new TH1D("ht", "Time Resolution;T (ps);Counts", 1e3, t.L, t.R);

	TF1 *fitT;

	double TAcor=0;

	double T[1000000]={0};
	const int chN=2;
	//EVENT** ch[chN];
	TH2D *hAT[chN];
    TF1 *fitAT= new TF1("fitAT", "0", U.at(0).L, U.at(0).R);
    for (int h=0; h<chN; h++ )
    {
        sprintf(buff,"hAT%d",h);
        hAT[h] = new TH2D(buff,buff, 200, U.at(h).L, U.at(h).R, 4e3, t.L, t.R);
    }     

    

	//ch[0]=&A;
	//ch[1]=&B;
	//ch[2]=&MCP;
    double L=20;
    //cout<<"the number of values of fraction of CFD: "<<L<<endl;
    for(int j=0;j<L;j++){
		//multi CFD fraction dicrimination
        
	for(int s =0;s<iter;s++){



		int N = t1->GetEntries();
		cout<< "Total entries is :"<<N<<endl;

		//correction of A1
		for (int h= 0;h<chN;h++){
		t.L=initial_t.L;
        t.R=initial_t.R;
		for(int i=0; i<N; i++)
		{
			t1->GetEntry(i);

			//Q1=(*A).Q;
			//Q2=(*B).Q;
			if (h-1<0)
				TAcor=fitAT->Eval((*ch.at(chN-1)).Q);
			else
				TAcor=fitAT->Eval((*ch.at(h-1)).Q);
			
			//TAcor=fitAT->Eval(Q2);
			//if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*A).rise<4&&(*B).rise<4){
			//if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*A).rise<4&&(*B).rise<4&&(*MCP).Q>0.3){
				//if(1){
				//if(1){
				//if((*A).A>0.3&&(*B).A>0.3){
                if(abs(p[j])<1e-12) return fitAT;
                else if(
                    (*ch.at(0)).Amp>cut.at(0).Amplow&&
                    (*ch.at(0)).Amp<cut.at(0).Ampup&&
                    (*ch.at(0)).Q>cut.at(0).Qlow&&
                    (*ch.at(0)).Q<cut.at(0).Qup&&
                    (*ch.at(0)).rise>cut.at(0).riselow&&
                    (*ch.at(0)).rise<cut.at(0).riseup&&
                    (*ch.at(1)).Amp>cut.at(1).Amplow&&
                    (*ch.at(1)).Amp<cut.at(1).Ampup&&
                    (*ch.at(1)).Q>cut.at(1).Qlow&&
                    (*ch.at(1)).Q<cut.at(1).Qup&&
                    (*ch.at(1)).rise>cut.at(1).riselow&&
                    (*ch.at(1)).rise<cut.at(1).riseup
                ){
    
				if (s==0&&h==0) 
					//T[i]=(*A).p20-(*B).p20-TAcor;
					T[i]=(*ch.at(0)).CFD[j]-(*ch.at(1)).CFD[j];
				else 
					T[i]=T[i]-TAcor;


				ht->Fill(T[i]);		
				hAT[h]->Fill((*ch.at(h)).Q,T[i]);
			
			}

		}


		cCor->cd();
		cCor->Clear();
		ht->Draw();
		//return;
		//tRL=-5;
		//tRR=5;

		//twoguasfit(ht,&tL,&tR,fac,rbt);
        fitT = gausfit(ht,rbt,fac,&t);
	//	sprintf(buff,"%s_TR_ch%d_cor%d.png",name,h,s);
        sprintf(buff,"%s_ch%d_CFDfrac%g_TR_cor%d.png",name.c_str(),h,p[j],s);
		cCor->SaveAs(buff);
        output << name.c_str() << "\t" << h << "\t" <<p[j] << "\t" << s  << "\t" << fitT->GetParameter(2)  << "\t" << fitT->GetParError(2) << endl;

        sprintf(buff,"%s_ch%d_CFDfrac%gAt_pfx_cor%d",name.c_str(),h,p[j],s);
		//sprintf(buff,"%s_At_ch%d_pfx_cor%d",name,h,s);

		fitAT=profilefit(hAT[h],rbU,rbt*8,t,U.at(h),buff);

		ht->Reset();
		hAT[h]->Reset();
		}
	}
fitAT= new TF1("fitAT", "0", U.at(0).L, U.at(0).R);
    
    }

	
        output.close();

}

TF1* LZWfunc::CH1Correction(TTree * t1, EVENT * A, EVENT * B, double *p, int rbU, int rbt, char *name, double fac, int iter, RANGE t, RANGE U)
    {
        //A channel with trigger
        //write for SIMU

        //Correction QT

        TCanvas *cCor = new TCanvas("cCor", "cCor", 800, 600);

        TH1D *ht = new TH1D("ht", "time Resolution;T (ps);Counts", 2e4, t.L, t.R);

        TH2D *hAT = new TH2D("hAT", "", 600, U.L, U.R, 2e4, t.L, t.R);

        TF1 *fitAT = new TF1("fitAT", "0", U.L, U.R);


        double TAcor = 0;
        double Q = 0;
        for (int j = 0; j < 20; j++)
        {
            double T[1000000] = {0};
            for (int s = 0; s < iter; s++)
            {

                int N = t1->GetEntries();
                cout << "Total entries is :" << N << endl;
                for (int i = 0; i < N; i++)
                {
                    t1->GetEntry(i);

                    Q = (*A).Amp;
                    TAcor = fitAT->Eval(Q);
                    //if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3){
                    //if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3&&(*A).rise<4&&(*B).rise<4){
                    if (abs(p[j]) < 1e-3)
                        return fitAT;
                    // else if((*A).A<0&&(*A).time>0){
                    else if ((*A).Amp < -0 && (*A).time > 0 && (*B).Amp < -0 && (*B).time > 0)
                    {

                        //if((*A).Q<600){
                        //if(1){
                        if (s == 0)
                            //T[i]=(*A).time-(*MCP).time-TAcor;
                            T[i] = (*A).time;
                        else
                            T[i] = T[i] - TAcor;

                        ht->Fill(T[i]);
                        //				cout<<"Q="<<Q<<" T="<<T[i]<<endl;
                        hAT->Fill(Q, T[i]);
                    }
                }

                //cout<<"<====process check====>"<<endl;

                cCor->cd();
                cCor->Clear();
                ht->Draw();
                //return;
                //t.L=-5;
                //t.R=5;
                t.L = ht->GetMean() - 3 * ht->GetRMS();
                t.R = ht->GetMean() + 3 * ht->GetRMS();
                //cout<<"iter="<<s<<", MEAN="<<ht->GetMean()<<"\t RMS="<<ht->GetRMS()<<endl;
                //cout<<"iter="<<s<<", tL="<<tL<<"\t tR="<<tR<<endl;
                //if (s==1) return;
                //twoguasfit(ht,&tL,&tR,fac,rbt);
                gausfit(ht, rbt, fac, &t);

                sprintf(buff, "%s_th%g_TR_cor%d.png", name, p[j], s);
                cCor->SaveAs(buff);

                sprintf(buff, "%s_th%gAt_pfx_cor%d", name, p[j], s);

                fitAT = profilefit(hAT, rbU, rbt * 2, t, U, buff);

                ht->Reset();
                hAT->Reset();
            }
            fitAT = new TF1("fitAT", "0", U.L, U.R);
        }

        return fitAT;
    
    }
void LZWfunc::CH1Correction(TTree * t1, EVENT * A, double * t0, RANGE t,RANGE U,int rbU, int rbt,  double fac, int iter)
    {
        //A channel with trigger
        //write for SIMU

        //Correction QT
        double Q; 
        RANGE initial_t;
        initial_t.L=-5;
        initial_t.R=5;
        sprintf(buff,"%s_ch1Correction.dat",name.c_str());
        output.open(buff,ios::trunc);
        TCanvas *cCor = new TCanvas("cCor", "cCor", 800, 600);

        TH1D *ht; 
        TH2D *hAT;
        
        TH1D *htraw= new TH1D("htraw", "time Resolution;T (ps);Counts", 2e4, t.L, t.R);

        TH2D *hATraw= new TH2D("hATraw", "", 600, U.L, U.R, 2e4, t.L, t.R);

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

                    Q = (*A).Amp;
                    TAcor = fitAT->Eval(Q);
                    //if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3){
                    //if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3&&(*A).rise<4&&(*B).rise<4){
                    // else if((*A).A<0&&(*A).time>0){
                    if ((*A).x > charcut.x.L &&
                    (*A).x < charcut.x.R &&
                    (*A).Q > charcut.q.L &&
                    (*A).Q <  charcut.q.R&&
                    (*A).rise > charcut.r.L &&
                    (*A).rise < charcut.r.R &&
                    (*A).bl > charcut.bl.L &&
                    (*A).bl < charcut.bl.R &&
                    (*A).blrms > charcut.blrms.L &&
                    (*A).blrms < charcut.blrms.R 
                      )
                    {

                        if((*A).Amp > charcut.y.L &&
                    (*A).Amp < charcut.y.R ){
                        if (s == 0)
                            //T[i]=(*A).time-(*MCP).time-TAcor;
                            {T[i] = (*A).time-*t0;
                            htraw->Fill(T[i]);
                            hATraw->Fill(Q, T[i]);
                            }
                        else
                            {T[i] = T[i] - TAcor;
                            htcor->Fill(T[i]);
                            hATcor->Fill(Q, T[i]);
                            }
                        //				cout<<"Q="<<Q<<" T="<<T[i]<<endl;
                    }
                    }
                }

                //cout<<"<====process check====>"<<endl;

                cCor->cd();
                cCor->Clear();
                if(s==0){
                    ht=(TH1D *)htraw->Clone();
                    hAT=(TH2D *)hATraw->Clone();
                }
                else{
                    ht=(TH1D *)htcor->Clone();
                    hAT=(TH2D *)hATcor->Clone();
                    cout<<"progress check"<<endl;
                }
                //return;
                //t.L=-5;
                //t.R=5;
                t.L = ht->GetMean() - 3 * ht->GetRMS();
                t.R = ht->GetMean() + 3 * ht->GetRMS();
                //cout<<"iter="<<s<<", MEAN="<<ht->GetMean()<<"\t RMS="<<ht->GetRMS()<<endl;
                //cout<<"iter="<<s<<", tL="<<tL<<"\t tR="<<tR<<endl;
                //if (s==1) return;
                //twoguasfit(ht,&tL,&tR,fac,rbt);
                fitT=gausfit(ht, rbt, fac, &t);

                sprintf(buff, "%s_TR_cor%d.png", name.c_str(), s);
                cCor->SaveAs(buff);

                output << name.c_str() << "\t" << s  << "\t" << fitT->GetParameter(2)  << "\t" << fitT->GetParError(2) << endl;

                sprintf(buff, "%s_At_pfx_cor%d", name.c_str(), s);

                fitAT = profilefit(hAT, rbU, rbt * 2, t, U, buff);

                ht->Reset();
                hAT->Reset();
                htraw->Reset();
                hATraw->Reset();
                htcor->Reset();
                hATcor->Reset();
                }
                
            
       

        
    
    }
void LZWfunc::drawcharacter(TTree *t3, int chN, string* chname, vector<charRANGE> chR){
	DrawMyfunc draw;

	for (int i=0;i<chN;i++)
	{
	TGaxis::SetMaxDigits(3);
    TH2F *hqy = new TH2F("hqy","",200,chR.at(i).y.L,chR.at(i).y.R,200,chR.at(i).q.L,chR.at(i).q.R);
    draw.Hist(hqy,"Amp (V)","Charge (pC)",1.5);
	sprintf(buff,"MCP%d_all_charge:MCP%d_global_maximum_y>>hqy",i+1,i+1);
	t3->Draw(buff,"","colz");
	sprintf(buff,"%s_Qvsy.png",chname[i].data());
	gPad->SaveAs(buff);
	

	TH2F* hqr=new TH2F("hqr","",200,chR.at(i).r.L,chR.at(i).r.R,200,chR.at(i).q.L,chR.at(i).q.R);
    draw.Hist(hqr,"risetime (ns)","Charge (pC)",1.5);
	sprintf(buff,"MCP%d_all_charge:MCP%d_rise_time>>hqr",i+1,i+1);
	t3->Draw(buff,"","colz");
	sprintf(buff,"%s_Qvsrise.png",chname[i].data());
	gPad->SaveAs(buff);
	
	
	TH2F* hqx=new TH2F("hqx","",200,chR.at(i).x.L,chR.at(i).x.R,200,chR.at(i).q.L,chR.at(i).q.R);
    draw.Hist(hqx,"time (ns)","Charge (pC)",1.5);
	sprintf(buff,"MCP%d_all_charge:MCP%d_global_maximum_x>>hqx",i+1,i+1);
	t3->Draw(buff,"","colz");
	sprintf(buff,"%s_Qvsx.png",chname[i].data());
	gPad->SaveAs(buff);

	
	TH1F* hbl=new TH1F("hbl","",200,chR.at(i).bl.L,chR.at(i).bl.R);
    draw.Hist(hbl,"baseline level(V)","Counts",1.5);
	sprintf(buff,"MCP%d_baseline_level>>hbl",i+1);
	t3->Draw(buff,"","colz");
	sprintf(buff,"%s_bl.png",chname[i].data());
	gPad->SaveAs(buff);
	
	
	TH1F* hblrms=new TH1F("hblrms","",200,chR.at(i).blrms.L,chR.at(i).blrms.R);
    draw.Hist(hblrms,"baseline RMS(V)","Counts",1.5);
	sprintf(buff,"MCP%d_baseline_rms>>hblrms",i+1);
	t3->Draw(buff,"","colz");
	sprintf(buff,"%s_blrms.png",chname[i].data());
	gPad->SaveAs(buff);

	}
	
	/*TH2F* hqy = new TH2F("hqy","",2e3,chR.at(i).y.L,chR.at(i).y.R,200,chR.at(i).q.L,chR.at(i).q.R);
	TH2F* hqr = new TH2F("hqr","",2e3,,10);
	TH1F* hrise = new TH1F("hy","",2e3,-10,10);
	TH1F* hy = new TH1F("hy","",2e3,-10,10);
	TH1F* hy = new TH1F("hy","",2e3,-10,10);
	TH1F* hy = new TH1F("hy","",2e3,-10,10);*/
	
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