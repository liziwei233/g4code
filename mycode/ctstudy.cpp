#include "TH1F.h"
#include "Include/DrawMyClass.h"
/*
char path[1024] = "/mnt/f/XOPtest/4Anode/crosstalk/A3";
char name[1024] = "MPE2";
    int idindex[4]={2,0,1,3};
    const double step = 50;
    const int cutstart = 200;
    const int cutend = 700;


char path[1024] = "/mnt/f/XOPtest/4Anode/crosstalk/A4";
char name[1024] = "MPE";
    int idindex[4]={3,0,1,2};
    const double step = 50;
    const int cutstart = 250;
    const int cutend = 800;
*/
/*
char path[1024] = "/mnt/f/XOPtest/4Anode/crosstalk/A4";
char name[1024] = "MPE2";
    int idindex[4]={3,0,1,2};
    const double step = 50;
    const int cutstart = 200;
    const int cutend = 600;
    */
   /*
char path[1024] = "/mnt/f/XOPtest/4Anode/Rb";
char name[1024] = "RbA4";
    int idindex[4]={3,0,1,2};
    const double step = 50;
    const int cutstart = 300;
    const int cutend = 900;   
*/
char path[1024] = "/mnt/f/XOPtest/4Anode/Rb";
char name[1024] = "A4re";
    int idindex[4]={3,0,1,2};
    const double step = 50;
    const int cutstart = 250;
    const int cutend = 600;   
int draw_txt1();

void ctstudy()
{

    

    double AMP[4];
    double InvAMP[4];
    double trueSignal;
    sprintf(buff, "%s/%soutput.dat", path, name);
    ofstream op(buff, ios::trunc);
    sprintf(buff, "%s/%s.root", path, name);
    if (gSystem->AccessPathName(buff))
    {
        cout << "Error! The file isnt exist!" << endl;
        return;
    }
    TFile *f1 = new TFile(buff, "read");
    TTree *t1 = (TTree *)f1->Get("Pico");
    for (int i = 0; i < 4; i++)
    {
        sprintf(buff, "MCP%d_global_maximum_y", i);
        cout << buff << endl;
        t1->SetBranchAddress(buff, &(AMP[i]));
        sprintf(buff, "MCP%d_invert_maximum_y", i);
        cout << buff << endl;
        t1->SetBranchAddress(buff, &(InvAMP[i]));
    }
    //return;
    TH1F *hct[3];
    TH1F *hctInv[3];
    TH1F *hctNei[3];
    TF1 *fit[3];
    TCanvas *c;
    for (int i = 0; i < 3; i++)
    {
        sprintf(buff, "hct%d", i);
        hct[i] = new TH1F(buff, ";CtAmp/TrueAmp", 0.51e3, -0.01, 0.5);
        sprintf(buff, "hctInv%d", i);
        hctInv[i] = new TH1F(buff, ";CtInvAmp/TrueAmp", 0.21e3, -0.01, 0.2);
        sprintf(buff, "hctNei%d", i);
        hctNei[i] = new TH1F(buff, ";CtInv/CtAmp", 1.e3, -0.01, 1);
    }

    int Nstep = (cutend - cutstart) / step;
    double cut[Nstep + 3];
    cut[0] = 0;
    for (int i = 1; i < Nstep + 2; i++)
    {
        cut[i] = cutstart + step * (i - 1);
    }
    cut[Nstep + 2] = 1.5e3;
    double min = 0;
    double max = 0;
    int N = t1->GetEntries();
    for (int j = 0; j <  Nstep+2; j++)
    {
        for (int s = 0; s < 3; s++)
            {
                hct[s]->Reset();
                hctInv[s]->Reset();
                hctNei[s]->Reset();
            }
        for (int i = 0; i < N; i++)
        {
            
            t1->GetEntry(i);
            trueSignal = AMP[idindex[0]];
            
            if (trueSignal > cut[j] && trueSignal < cut[j + 1])
            {
                for (int i = 0; i < 3; i++)
                {
                    hct[i]->Fill(AMP[idindex[1+i]] / trueSignal);
                    hctInv[i]->Fill(-1 * InvAMP[idindex[ 1 + i]] / trueSignal);
                    hctNei[i]->Fill(-1 * InvAMP[idindex[1 + i]] / AMP[idindex[1 + i]]);
                }
            }
        }
        for (int i = 0; i < 3; i++)
        {
            fit[0]=NULL;
            fit[1]=NULL;
            fit[2]=NULL;
            c = cdC(i * 3);
            fit[0] = gausfit(hct[i], 0.1, 3, 2, 2, 0, 0.5);
            sprintf(buff, "%s/cut%g_%gct%d.png", path, cut[j], cut[j + 1], i);
            c->SaveAs(buff);
            if (fit[0])
                op << fit[0]->GetParameter(1) << "\t";
            else
                op << -999 << "\t";
            c = cdC(i * 3 + 1);
            fit[1] = gausfit(hctInv[i], 0.1, 3, 2, 1, 0, 0.2);
            sprintf(buff, "%s/cut%g_%gctInv%d.png", path, cut[j], cut[j + 1], i);
            c->SaveAs(buff);
            if (fit[1])
                op << fit[1]->GetParameter(1) << "\t";
            else
                op << -999 << "\t";

            c = cdC(i * 3 + 2);
            fit[2] = gausfit(hctNei[i], 0.2, 3, 2, 4, 0, 1);
            sprintf(buff, "%s/cut%g_%gctNei%d.png", path, cut[j], cut[j + 1], i);
            c->SaveAs(buff);
            if (fit[2])
                op << fit[2]->GetParameter(1) << "\t";
            else
                op << -999 << "\t";
        }
        op << endl;
    }
    draw_txt1();
}
int draw_txt1()
{
    TGraphErrors *g[10];
    
    //sprintf(path, "/mnt/f/XOPtest/4Anode/crosstalk/A2");
    setgStyle();
    gStyle->SetOptFit(000);
    const int N = 100;
    double ct[3][N];
    double ctInv[3][N];
    double ctNei[3][N];
    
    
    int Nstep = (cutend - cutstart) / step;
    double Amp[Nstep + 2];
    for (int i = 0; i < Nstep + 2; i++)
    {
        Amp[i] = cutstart + step * i;
    }
    double Amperr[N] = {0};

    //Size_t size[18] = {2, 1.7, 2.1, 2.1, 1.8, 2, 2.1, 2.1, 1.8, 2, 2.1, 2.1, 1.8, 2, 2.1, 2.1, 1.8, 2};
    //Color_t color[18] = {1, 2, 4, kOrange - 1, 6, 8, 4, kOrange - 1, 6, 8, 4, kOrange - 1, 6, 8, 4, kOrange - 1, 6, 8};
    //Style_t marker[18] = {20, 21, 22, 23, 25, 26, 22, 23, 25, 26, 22, 23, 25, 26, 22, 23, 25, 26};
    Size_t size[18] = {2,1.7,2.1,2,2, 1.7,1.7,1.7, 2.1,2.1,2.1};
    Color_t color[18] = {1,2,4,1,1, 2,2,2, 4,4,4};
    Style_t marker[18] = {20,21,22,20,21,22,20,21,22};

    ifstream input;
    char filename[1024];
    //char name[1024];

    sprintf(buff, "%s/%soutput.dat", path, name);
    input.open(buff);
    cout<<"Open your file:"<<buff<<endl;
    if (!input)
    {
        cout << "file can't be found: \n"
             << buff << endl;
        return 0;
    }

    int k = 0;
    //return 1;
    while (input >> ct[0][k] >> ctInv[0][k] >> ctNei[0][k]&& !input.eof())
    {
        for(int i=1;i<3;i++)
        {
            input >> ct[i][k] >> ctInv[i][k] >> ctNei[i][k];
        }
        cout << "enter the file: " << buff << endl;

        cout << ct[0][k] << "\t" << ct[1][k] << "\t" << ct[2][k]  << endl;
        k++;
        //if(input.eof()) break;
    }
    input.close();

    

    int n = k;
    for(int i=0;i<n;i++)
    {
        Amp[i]=Amp[i]*0.2586;
        cout<<Amp[i]<<endl;
    }

    TCanvas *c = cdC(0);
    DrawMyPad(gPad, "True Signal Amp (mV)", "Crosstalk Ratio ",  10, 300,0, 0.4);
    //DrawMyPad(gPad, xtitle, ytitle, x_low, x_high, y_low * 0.7, y_high * 1.1);
    //return 1;
    TLegend *leg = DrawMyLeg(0.6, 0.7, 0.95, 0.9, 62, 0.03);
    double fac=0;
    for (int i = 0; i < 3; i++)
    {
        
        g[i] = new TGraphErrors(n,  Amp,ct[i], 0, 0);
        DrawMyGraph(g[i], "", "", 1.2, marker[i], color[i]);
        g[i]->Draw("sameP");
        sprintf(buff,"AMP%d/AMP%d",idindex[i+1]+1,idindex[0]+1);
        leg->AddEntry(g[i], buff, "lp");
    }
    leg->Draw();
    sprintf(buff, "%s/ct.png", path);
    TGaxis::SetMaxDigits(3);
    gPad->SaveAs(buff);

    leg->Clear();
    c = cdC(1);
    DrawMyPad(gPad, "True Signal Amp (mV)", "Crosstalk Ratio ",  10, 300,0, 0.21);
     for (int i = 0; i < 3; i++)
    {
        g[i] = new TGraphErrors(n, Amp, ctInv[i], 0, 0);
        DrawMyGraph(g[i], "", "", 1.2, marker[i], color[i]);
        g[i]->Draw("sameP");
        sprintf(buff,"InvAMP%d/AMP%d",idindex[i+1]+1,idindex[0]+1);
        leg->AddEntry(g[i], buff, "lp");
    }
    leg->Draw();
sprintf(buff, "%s/ctInvert.png", path);
    gPad->SaveAs(buff);

    leg->Clear();
   c = cdC(2);
    DrawMyPad(gPad, "True Signal Amp (mV)", "Crosstalk Ratio ",  10, 300,0, 0.8);
    for (int i = 0; i < 3; i++)
    {
        g[i] = new TGraphErrors(n, Amp, ctNei[i], 0, 0);
        DrawMyGraph(g[i], "", "", 1.2, marker[i], color[i]);
        g[i]->Draw("sameP");
        sprintf(buff,"InvAMP%d/AMP%d",idindex[i+1]+1,idindex[i+1]+1);
        leg->AddEntry(g[i], buff, "lp");
    }
    leg->Draw();
    sprintf(buff, "%s/ctneighbour.png", path);
    gPad->SaveAs(buff);
    return 1;
}
void ctstudy2()
{

    
    sprintf(path,"%s","/mnt/f/R10754/MPE-crosstalk");
    sprintf(name,"%s","MPE-KT0864-2150V-A5A6A1A16-all");

    double AMP[4];
    double InvAMP[4];
    sprintf(buff, "%s/%soutput.dat", path, name);
    ofstream op(buff, ios::trunc);
    sprintf(buff, "%s/%s.root", path, name);
    if (gSystem->AccessPathName(buff))
    {
        cout << "Error! The file isnt exist:\t"<<buff << endl;
        return;
    }
    TFile *f1 = new TFile(buff, "read");
    TTree *t1 = (TTree *)f1->Get("Pico");
    for (int i = 1; i < 5; i++)
    {
        sprintf(buff, "MCP%d_global_maximum_y", i);
        cout << buff << endl;
        t1->SetBranchAddress(buff, &(AMP[i-1]));
        sprintf(buff, "MCP%d_invert_maximum_y", i);
        cout << buff << endl;
        t1->SetBranchAddress(buff, &(InvAMP[i-1]));
    }
    //return;
    TH1F *hct[3];
    TH1F *hctInv[3];
    TH1F *hctNei[3];
    TF1 *fit[3];
    TCanvas *c;
    for (int i = 0; i < 3; i++)
    {
        sprintf(buff, "hct%d", i);
        hct[i] = new TH1F(buff, ";CtAmp/TrueAmp", 1.e3, -0.01, 1);
        sprintf(buff, "hctInv%d", i);
        hctInv[i] = new TH1F(buff, ";CtInvAmp/TrueAmp", 1.e3, -0.01, 1);
        sprintf(buff, "hctNei%d", i);
        hctNei[i] = new TH1F(buff, ";CtInv/CtAmp", 1.e3, -0.01, 2);
    }

    double step = 0.015;
    double start = 0.035;
    double end = 0.23;
    int Nstep = (end - start) / step;
    double cut[Nstep + 3];
    cut[0] = 0;
    for (int i = 1; i < Nstep + 2; i++)
    {
        cut[i] = start + step * (i - 1);
    }
    cut[Nstep + 2] = 1.5e3;
    double min = 0;
    double max = 0;
    int N = t1->GetEntries();
    for (int j = 0; j < Nstep+ 2; j++)
    {
        for (int s = 0; s < 3; s++)
            {
                hct[s]->Reset();
                hctInv[s]->Reset();
                hctNei[s]->Reset();
            }
        for (int i = 0; i < N; i++)
        {
            
            t1->GetEntry(i);

            if (AMP[0] > cut[j] && AMP[0] < cut[j + 1])
            {
                for (int i = 0; i < 3; i++)
                {
                    hct[i]->Fill(AMP[1 + i] / AMP[0]);
                    hctInv[i]->Fill(-1 * InvAMP[1 + i] / AMP[0]);
                    hctNei[i]->Fill(-1 * InvAMP[1 + i] / AMP[1 + i]);
                }
            }
        }
        for (int i = 0; i < 3; i++)
        {
            fit[0]=NULL;
            fit[1]=NULL;
            fit[2]=NULL;
            c = cdC(i * 3);
            fit[0] = gausfit(hct[i], 0.01, 3, 2, 2, 0, 0.5);
            sprintf(buff, "%s/cut%g_%gct%d.png", path, cut[j], cut[j + 1], i);
            c->SaveAs(buff);
            if (fit[0])
                op << fit[0]->GetParameter(1) << "\t";
            else
                op << -999 << "\t";
            c = cdC(i * 3 + 1);
            fit[1] = gausfit(hctInv[i], 0.01, 3, 2, 1, 0, 0.2);
            sprintf(buff, "%s/cut%g_%gctInv%d.png", path, cut[j], cut[j + 1], i);
            c->SaveAs(buff);
            if (fit[1])
                op << fit[1]->GetParameter(1) << "\t";
            else
                op << -999 << "\t";

            c = cdC(i * 3 + 2);
            fit[2] = gausfit(hctNei[i], 0.01, 3, 2, 2, 0, 1);
            sprintf(buff, "%s/cut%g_%gctNei%d.png", path, cut[j], cut[j + 1], i);
            c->SaveAs(buff);
            if (fit[2])
                op << fit[2]->GetParameter(1) << "\t";
            else
                op << -999 << "\t";
        }
        op << endl;
    }
}



int draw_txt2()
{
TGraphErrors *g[10];
    
    sprintf(path, "/mnt/f/R10754/MPE-crosstalk");
    setgStyle();
    gStyle->SetOptFit(000);
    const int N = 100;
    double ct[3][N];
    double ctInv[3][N];
    double ctNei[3][N];
    int chID[4]={5,6,1,16};
    double step = 0.015;
    double start = 0.035;
    double end = 0.23;
    int Nstep = (end - start) / step;
    double Amp[Nstep + 2];
    for (int i = 0; i < Nstep + 2; i++)
    {
        Amp[i] = start + step * i;
    }
    double Amperr[N] = {0};

    //Size_t size[18] = {2, 1.7, 2.1, 2.1, 1.8, 2, 2.1, 2.1, 1.8, 2, 2.1, 2.1, 1.8, 2, 2.1, 2.1, 1.8, 2};
    //Color_t color[18] = {1, 2, 4, kOrange - 1, 6, 8, 4, kOrange - 1, 6, 8, 4, kOrange - 1, 6, 8, 4, kOrange - 1, 6, 8};
    //Style_t marker[18] = {20, 21, 22, 23, 25, 26, 22, 23, 25, 26, 22, 23, 25, 26, 22, 23, 25, 26};
    Size_t size[18] = {2,1.7,2.1,2,2, 1.7,1.7,1.7, 2.1,2.1,2.1};
    Color_t color[18] = {1,2,4,1,1, 2,2,2, 4,4,4};
    Style_t marker[18] = {20,21,22,20,21,22,20,21,22};

    ifstream input;
    char filename[1024];
    //char name[1024];

    sprintf(buff, "%s/%s.dat", path, "MPE-KT0864-2150V-A5A6A1A16-alloutput");
    input.open(buff);
    if (!input)
    {
        cout << "file can't be found: \n"
             << buff << endl;
        return 0;
    }

    int k = 0;
    //return 1;
    while (input >> ct[0][k] >> ctInv[0][k] >> ctNei[0][k]&& !input.eof())
    {
        for(int i=1;i<3;i++)
        {
            input >> ct[i][k] >> ctInv[i][k] >> ctNei[i][k];
        }
        cout << "enter the file: " << buff << endl;

        cout << ct[0][k] << "\t" << ct[1][k] << "\t" << ct[2][k]  << endl;
        k++;
        //if(input.eof()) break;
    }
    input.close();

    

    int n = k;
    for(int i=0;i<n;i++)
    {
        Amp[i]=Amp[i]*1e3;
        cout<<Amp[i]<<endl;
    }

    TCanvas *c = cdC(0);
    DrawMyPad(gPad, "True Signal Amp (mV)", "Crosstalk Ratio ",  0, 380,0, 0.4);
    //DrawMyPad(gPad, xtitle, ytitle, x_low, x_high, y_low * 0.7, y_high * 1.1);
    //return 1;
    TLegend *leg = DrawMyLeg(0.6, 0.7, 0.95, 0.9, 62, 0.03);
    double fac=0;
    for (int i = 0; i < 3; i++)
    {
        
        g[i] = new TGraphErrors(n,  Amp,ct[i], 0, 0);
        DrawMyGraph(g[i], "", "", 1.2, marker[i], color[i]);
        g[i]->Draw("sameP");
        sprintf(buff,"AMP%d/AMP%d",chID[i+1],chID[0]);
        leg->AddEntry(g[i], buff, "lp");
    }
    leg->Draw();
    sprintf(buff, "%s/ct.png", path);
    TGaxis::SetMaxDigits(3);
    gPad->SaveAs(buff);

    leg->Clear();
    c = cdC(1);
    DrawMyPad(gPad, "True Signal Amp (mV)", "Crosstalk Ratio ",   0, 380,0, 0.21);
     for (int i = 0; i < 3; i++)
    {
        g[i] = new TGraphErrors(n, Amp, ctInv[i], 0, 0);
        DrawMyGraph(g[i], "", "", 1.2, marker[i], color[i]);
        g[i]->Draw("sameP");
        sprintf(buff,"InvAMP%d/AMP%d",chID[i+1],chID[0]);
        leg->AddEntry(g[i], buff, "lp");
    }
    leg->Draw();
sprintf(buff, "%s/ctInvert.png", path);
    gPad->SaveAs(buff);

    leg->Clear();
   c = cdC(2);
    DrawMyPad(gPad, "True Signal Amp (mV)", "Crosstalk Ratio ",   0, 380,0.15, 2);
    for (int i = 0; i < 3; i++)
    {
        g[i] = new TGraphErrors(n, Amp, ctNei[i], 0, 0);
        DrawMyGraph(g[i], "", "", 1.2, marker[i], color[i]);
        g[i]->Draw("sameP");
        sprintf(buff,"InvAMP%d/AMP%d",chID[i+1],chID[i+1]);
        leg->AddEntry(g[i], buff, "lp");
    }
    leg->Draw();
    sprintf(buff, "%s/ctneighbour.png", path);
    gPad->SaveAs(buff);
    return 1;
}