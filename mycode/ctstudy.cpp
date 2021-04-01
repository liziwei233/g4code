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
char path[1024] = "/mnt/f/R10754/KT0881/newbase-crosstalkstudy/different_ground";
char name[1024] = "4chbase-HV2200-MPE-half-D";

const double step = 0.015;
const double cutstart = 0.03;
const double cutend = 0.14;
int idindex[2] = {2, 4};

//char path[1024] = "/mnt/f/XOPtest/4Anode/crosstalk";
//char name[1024] = "crosstalkA3MPE2";
//const double step = 50;
//const int cutstart = 200;
//const int cutend = 700;
//int idindex[4] = {2, 0, 1, 3};
int draw_txt1();

void ctstudy()
{

    gStyle->SetOptFit(111);
    double AMP[4];
    double ctAMP[4];
    double ringAMP[4];
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
        t1->SetBranchAddress(buff, &(ctAMP[i]));
        sprintf(buff, "MCP%d_secondinvertpeak_y", i);
        cout << buff << endl;
        t1->SetBranchAddress(buff, &(ringAMP[i]));
    }
    //return;
    TH1F *hct[4];
    TH1F *hring[4];
    TF1 *fit[4];
    TH1F *htemp;
    TCanvas *c;
    for (int i = 0; i < 4; i++)
    {
        sprintf(buff, "hct%d", i);
        hct[i] = new TH1F(buff, ";Crosstalk", 0.51e3, -0.01, 0.2);
        sprintf(buff, "hring%d", i);
        hring[i] = new TH1F(buff, ";Ringing", 0.21e3, -0.01, 0.5);
    }

    int Nstep = (cutend - cutstart) / step;
    double cut[Nstep + 3];
    cut[0] = 0.01;
    for (int i = 1; i < Nstep + 2; i++)
    {
        cut[i] = cutstart + step * (i - 1);
    }
    cut[Nstep + 2] = 1.5e3;
    double min = 0;
    double max = 0;
    int N = t1->GetEntries();
    for (int j = 0; j < Nstep + 2; j++)
    {
        for (int s = 0; s < 4; s++)
        {
            hct[s]->Reset();
            hring[s]->Reset();
        }
        for (int i = 0; i < N; i++)
        {

            t1->GetEntry(i);
            trueSignal = AMP[idindex[0]];

            if (trueSignal > cut[j] && trueSignal < cut[j + 1])
            {
                for (int i = 0; i < 4; i++)
                {
                    hring[i]->Fill(abs(ringAMP[idindex[i]]) / trueSignal);
                    hct[i]->Fill(abs(ctAMP[idindex[i]]) / trueSignal);
                }
            }
        }
        for (int i = 0; i < 4; i++)
        {
            fit[0] = NULL;
            fit[1] = NULL;
            c = cdC(i * 2);
            htemp = (TH1F *)gausfit(hct[i], 0.01, 3, 2, 2, 0, 0.2);
            fit[0] = (TF1 *)htemp->GetFunction("fitU");
            sprintf(buff, "%s/cut%g_%gct%d.png", path, cut[j], cut[j + 1], i);
            c->SaveAs(buff);
            if (fit[0])
                op << fit[0]->GetParameter(1) << "\t";
            else
                op << -999 << "\t";
            c = cdC(i * 2 + 1);
            htemp = (TH1F *)gausfit(hring[i], 0.01, 3, 2, 1, 0, 0.5);
            fit[1] = (TF1 *)htemp->GetFunction("fitU");
            sprintf(buff, "%s/cut%g_%gring%d.png", path, cut[j], cut[j + 1], i);
            c->SaveAs(buff);
            if (fit[1])
                op << fit[1]->GetParameter(1) << "\t";
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
    double ct[4][N];
    double ring[4][N];

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
    Size_t size[18] = {1.2, 1.1, 1.3, 1.3, 2, 1.7, 1.7, 1.7, 2.1, 2.1, 2.1};
    Color_t color[18] = {1, 2, 4, 8, 1, 2, 2, 2, 4, 4, 4};
    Style_t marker[18] = {20, 21, 22, 29, 21, 22, 20, 21, 22};

    ifstream input;
    char filename[1024];
    //char name[1024];

    sprintf(buff, "%s/%soutput.dat", path, name);
    input.open(buff);
    cout << "Open your file:" << buff << endl;
    if (!input)
    {
        cout << "file can't be found: \n"
             << buff << endl;
        return 0;
    }

    int k = 0;
    //return 1;
    while (input >> ct[0][k] >> ring[0][k] && !input.eof())
    {
        for (int i = 1; i < 4; i++)
        {
            input >> ct[i][k] >> ring[i][k];
        }
        cout << "enter the file: " << buff << endl;

        cout << ct[0][k] << "\t" << ct[1][k] << "\t" << ct[2][k] << "\t" << ct[3][k] << endl;
        k++;
        //if(input.eof()) break;
    }
    input.close();

    int n = k;
    for (int i = 0; i < n; i++)
    {
        Amp[i] = Amp[i] * 0.2586;
        cout << Amp[i] << endl;
    }

    TCanvas *c = cdC(0);
    DrawMyPad(gPad, "True Signal Amp (mV)", "Crosstalk Ratio ", 10, 300, 0, 0.21);
    //DrawMyPad(gPad, xtitle, ytitle, x_low, x_high, y_low * 0.7, y_high * 1.1);
    //return 1;
    TLegend *leg = DrawMyLeg(0.6, 0.6, 0.95, 0.9, 62, 0.04);
    double fac = 0;
    for (int i = 0; i < 4; i++)
    {

        g[i] = new TGraphErrors(n, Amp, ct[i], 0, 0);
        DrawMyGraph(g[i], "", "", size[i], marker[i], color[i]);
        g[i]->Draw("sameP");
        sprintf(buff, "crosstalk-A%d", idindex[i] + 1);
        leg->AddEntry(g[i], buff, "lp");
    }
    leg->Draw();
    sprintf(buff, "%s/%sct.png", path, name);
    TGaxis::SetMaxDigits(3);
    gPad->SaveAs(buff);

    leg->Clear();
    c = cdC(1);
    DrawMyPad(gPad, "True Signal Amp (mV)", "Ringing Ratio ", 10, 300, 0, 0.5);
    for (int i = 0; i < 4; i++)
    {
        g[i] = new TGraphErrors(n, Amp, ring[i], 0, 0);
        DrawMyGraph(g[i], "", "", size[i], marker[i], color[i]);
        g[i]->Draw("sameP");
        sprintf(buff, "Ringing-A%d", idindex[i] + 1);
        leg->AddEntry(g[i], buff, "lp");
    }
    leg->Draw();
    sprintf(buff, "%s/ring.png", path);
    gPad->SaveAs(buff);

    return 1;
}

void ctstudy2()
{

    double AMP[4];
    double InvAMP[4];
    double Inv2ndAMP[4];
    double charge[4];
    double xpos[4];

    sprintf(buff, "%s/%soutput.dat", path, name);
    ofstream op(buff, ios::trunc);
    sprintf(buff, "%s/%s.root", path, name);
    if (gSystem->AccessPathName(buff))
    {
        cout << "Error! The file isnt exist:\t" << buff << endl;
        return;
    }
    TFile *f1 = new TFile(buff, "read");
    TTree *t1 = (TTree *)f1->Get("Pico");
    //int CHNum[2] = {2, 4};
    int N = sizeof(idindex) / sizeof(idindex[0]);

    t1->SetBranchAddress("MCP2_all_charge", charge);
    for (int i = 0; i < N; i++)
    {
        sprintf(buff, "MCP%d_global_maximum_x", idindex[i]);
        cout << buff << endl;
        t1->SetBranchAddress(buff, &(xpos[i]));

        sprintf(buff, "MCP%d_global_maximum_y", idindex[i]);
        cout << buff << endl;
        t1->SetBranchAddress(buff, &(AMP[i]));
        sprintf(buff, "MCP%d_invert_maximum_y", idindex[i]);
        cout << buff << endl;
        t1->SetBranchAddress(buff, &(InvAMP[i]));
        sprintf(buff, "MCP%d_secondinvertpeak_y", idindex[i]);
        cout << buff << endl;
        t1->SetBranchAddress(buff, &(Inv2ndAMP[i]));
    }
    //return;
    TH1F *hct[4];
    TH1F *hring[4];
    TH1F *h2ndring[4];
    TH1F *htemp;
    TF1 *fit[4];
    TCanvas *c;
    for (int i = 0; i < N; i++)
    {
        sprintf(buff, "hct%d", i);
        hct[i] = new TH1F(buff, ";ctAmp/TrueAmp", 5.e3, -0.01, 1.01);

        sprintf(buff, "hring%d", i);
        hring[i] = new TH1F(buff, ";ringAmp/TrueAmp", 5e3, -0.01, 1.01);

        sprintf(buff, "h2ndring%d", i);
        h2ndring[i] = new TH1F(buff, ";2ndringAmp/TrueAmp", 5e3, -0.01, 1.01);
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
    int n = t1->GetEntries();
    for (int j = 0; j < Nstep + 2; j++)
    {
        cout << "cut: " << cut[j] << "\t" << cut[j + 1] << endl;
        for (int s = 0; s < N; s++)
        {
            hct[s]->Reset();
            hring[s]->Reset();
            h2ndring[s]->Reset();
        }
        for (int i = 0; i < n; i++)
        {

            t1->GetEntry(i);
            //cout<<charge[0]<<"\t"<<AMP[0]<<endl;
            if (charge[0] > 0.04 && AMP[0] > cut[j] && AMP[0] < cut[j + 1])
            {
                for (int i = 0; i < N; i++)
                {
                    //if (xpos[i] > 40 && xpos[i] < 43)
                    //if (xpos[i] > 23 && xpos[i] < 27)
                    if (xpos[i] > 75 && xpos[i] < 80)
                    {

                        hct[i]->Fill(-1 * InvAMP[i] / AMP[0]);
                        hring[i]->Fill(AMP[i] / AMP[0]);
                        h2ndring[i]->Fill(-1 * Inv2ndAMP[i] / AMP[0]);
                    }
                }
            }
        }
        //hct[0]->Draw();
        //return;
        for (int i = 0; i < N; i++)
        {
            fit[0] = NULL;
            fit[1] = NULL;
            fit[2] = NULL;
            c = cdC(i * 3);
            double UL = hct[i]->GetBinCenter(hct[i]->GetMaximumBin() - 2e2);
            double UR = hct[i]->GetBinCenter(hct[i]->GetMaximumBin() + 2e2);
            cout << "maximumbin: " << hct[i]->GetMaximumBin() << endl;
            cout << UL << "\t" << UR << endl;
            htemp = (TH1F *)gausfit(hct[i], 0.01, 3, 2, 2, UL, UR);
            if (!htemp)
            {
                cout << "the hct fit hist is NULL " << endl;
                op << hct[i]->GetBinCenter(hring[i]->GetMaximumBin()) << "\t";
                //return;
            }
            else
            {

                fit[0] = (TF1 *)htemp->GetFunction("fitU");
                sprintf(buff, "%s/cut%g_%gct%d.png", path, cut[j], cut[j + 1], i);
                c->SaveAs(buff);
                if (fit[0])
                    op << fit[0]->GetParameter(1) << "\t";
                else
                    op << -999 << "\t";
            }
            c = cdC(i * 3 + 1);
            UL = hring[i]->GetBinCenter(hring[i]->GetMaximumBin() - 2e2);
            UR = hring[i]->GetBinCenter(hring[i]->GetMaximumBin() + 2e2);
            cout << "maximumbin: " << hring[i]->GetMaximumBin() << endl;
            cout << UL << "\t" << UR << endl;
            htemp = (TH1F *)gausfit(hring[i], 0.01, 3, 2, 1, UL, UR);
            if (!htemp)
            {
                cout << "the hringfit hist is NULL " << endl;
                op << hring[i]->GetBinCenter(hring[i]->GetMaximumBin()) << "\t";
                //return;
            }
            else
            {

                fit[1] = (TF1 *)htemp->GetFunction("fitU");
                sprintf(buff, "%s/cut%g_%gring%d.png", path, cut[j], cut[j + 1], i);
                c->SaveAs(buff);
                if (fit[1])
                    op << fit[1]->GetParameter(1) << "\t";
                else
                    op << -999 << "\t";
            }

            c = cdC(i * 3 + 2);
            UL = h2ndring[i]->GetBinCenter(h2ndring[i]->GetMaximumBin() - 2e2);
            UR = h2ndring[i]->GetBinCenter(h2ndring[i]->GetMaximumBin() + 2e2);
            cout << "maximumbin: " << h2ndring[i]->GetMaximumBin() << endl;
            cout << UL << "\t" << UR << endl;
            htemp = (TH1F *)gausfit(h2ndring[i], 0.01, 3, 2, 2, UL, UR);
            if (!htemp)
            {
                cout << "the h2ndring hist is NULL " << endl;
                op << h2ndring[i]->GetBinCenter(hring[i]->GetMaximumBin()) << "\t";
                //return;
            }
            else
            {

                fit[2] = (TF1 *)htemp->GetFunction("fitU");
                sprintf(buff, "%s/cut%g_%g2ndring.png", path, cut[j], cut[j + 1], i);
                c->SaveAs(buff);
                if (fit[2])
                    op << fit[2]->GetParameter(1) << "\t";
                else
                    op << -999 << "\t";
            }
        }
        op << endl;
    }
}

int draw_txt2()
{
    TGraphErrors *g[10];
    //sprintf(path, "%s", "/mnt/f/R10754/KT0881/newbase-crosstalkstudy");
    //sprintf(name, "%s", "HV2100-A4-DN-close470pF");
    setgStyle();
    gStyle->SetOptFit(000);
    const int N = 100;
    double ct[4][N];
    double ring[4][N];
    double ring2nd[4][N];
    //int chID[2] = {2, 4};
    //double step = 0.015;
    //double start = 0.035;
    //double end = 0.23;
    int CHN = sizeof(idindex) / sizeof(idindex[0]);
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
    Size_t size[18] = {2, 1.7, 2.1, 2, 2, 1.7, 1.7, 1.7, 2.1, 2.1, 2.1};
    Color_t color[18] = {1, 2, 4, 1, 1, 2, 2, 2, 4, 4, 4};
    Style_t marker[18] = {20, 21, 22, 20, 21, 22, 20, 21, 22};

    ifstream input;
    char filename[1024];
    //char name[1024];
    sprintf(buff, "%s/%soutput.dat", path, name);
    input.open(buff);
    if (!input)
    {
        cout << "file can't be found: \n"
             << buff << endl;
        return 0;
    }

    int k = 0;
    //return 1;
    while (input >> ct[0][k] >> ring[0][k] >> ring2nd[0][k] && !input.eof())
    {
        for (int i = 1; i < CHN; i++)
        {
            input >> ct[i][k] >> ring[i][k] >> ring2nd[i][k];
        }
        cout << "enter the file: " << buff << endl;

        cout << ct[0][k] << "\t" << ct[1][k] << "\t" << ct[2][k] << endl;
        k++;
        //if(input.eof()) break;
    }
    input.close();

    int n = k;
    for (int i = 0; i < n; i++)
    {
        Amp[i] = Amp[i] * 1e3;
        cout << Amp[i] << endl;
    }

    TCanvas *c = cdC(0);
    DrawMyPad(gPad, "True Signal Amp (mV)", "Crosstalk Ratio ", (cutstart - step) * 1e3, (cutend + 2 * step) * 1e3, 0, 0.28);
    //DrawMyPad(gPad, xtitle, ytitle, x_low, x_high, y_low * 0.7, y_high * 1.1);
    //return 1;
    TLegend *leg = DrawMyLeg(0.6, 0.7, 0.95, 0.9, 62, 0.03);
    double fac = 0;
    for (int i = 0; i < CHN; i++)
    {

        g[i] = new TGraphErrors(n, Amp, ct[i], 0, 0);
        DrawMyGraph(g[i], "", "", 1.2, marker[i], color[i]);
        g[i]->Draw("sameP");
        sprintf(buff, "ct%d/AMP%d", idindex[i], idindex[0]);
        leg->AddEntry(g[i], buff, "lp");
    }
    leg->Draw();
    sprintf(buff, "%s/%sct.png", path, name);
    TGaxis::SetMaxDigits(3);
    gPad->SaveAs(buff);

    leg->Clear();
    c = cdC(1);
    DrawMyPad(gPad, "True Signal Amp (mV)", "Ringing Ratio ", (cutstart - step) * 1e3, (cutend + 2 * step) * 1e3, 0, 0.28);
    for (int i = 0; i < CHN; i++)
    {
        g[i] = new TGraphErrors(n, Amp, ring[i], 0, 0);
        DrawMyGraph(g[i], "", "", 1.2, marker[i], color[i]);
        g[i]->Draw("sameP");
        sprintf(buff, "ring%d/AMP%d", idindex[i], idindex[0]);
        leg->AddEntry(g[i], buff, "lp");
    }
    leg->Draw();
    sprintf(buff, "%s/%sring.png", path, name);
    gPad->SaveAs(buff);

    leg->Clear();
    c = cdC(2);
    DrawMyPad(gPad, "True Signal Amp (mV)", "2nd Ring Ratio ", (cutstart - step) * 1e3, (cutend + 2 * step) * 1e3, 0, 0.28);
    for (int i = 0; i < CHN; i++)
    {
        g[i] = new TGraphErrors(n, Amp, ring2nd[i], 0, 0);
        DrawMyGraph(g[i], "", "", 1.2, marker[i], color[i]);
        g[i]->Draw("sameP");
        sprintf(buff, "ring2ndAMP%d/AMP%d", idindex[i], idindex[0]);
        leg->AddEntry(g[i], buff, "lp");
    }
    leg->Draw();
    sprintf(buff, "%s/%sring2nd.png", path, name);
    gPad->SaveAs(buff);
    return 1;
}

int draw_txt3()
{
    TGraphErrors *g[10];
    sprintf(path, "%s", "/mnt/f/R10754/KT0881/newbase-crosstalkstudy/results");
    sprintf(name, "%s", "ImpactsOfPosOfCap");
    setgStyle();
    gStyle->SetOptFit(000);
    const int N = 100;
    double sigct, sigring, sigring2nd;
    double ct[10][N];
    double ring[10][N];
    double ring2nd[10][N];
    double Amp[10][N];
    int pointN[10];

    double yl=0;
    double yh=0.1;
    //int chID[2] = {2, 4};
    //double step = 0.015;
    //double start = 0.035;
    //double end = 0.23;
    int CHN = sizeof(idindex) / sizeof(idindex[0]);
    /*
    int Nstep = (cutend - cutstart) / step;
    double Amp[Nstep + 2];
    for (int i = 0; i < Nstep + 2; i++)
    {
        Amp[i] = cutstart + step * i;
    }
    double Amperr[N] = {0};
*/
    //Size_t size[18] = {2, 1.7, 2.1, 2.1, 1.8, 2, 2.1, 2.1, 1.8, 2, 2.1, 2.1, 1.8, 2, 2.1, 2.1, 1.8, 2};
    //Color_t color[18] = {1, 2, 4, kOrange - 1, 6, 8, 4, kOrange - 1, 6, 8, 4, kOrange - 1, 6, 8, 4, kOrange - 1, 6, 8};
    //Style_t marker[18] = {20, 21, 22, 23, 25, 26, 22, 23, 25, 26, 22, 23, 25, 26, 22, 23, 25, 26};
    Size_t size[18] = {2, 1.7, 2.1, 2, 2, 1.7, 1.7, 1.7, 2.1, 2.1, 2.1};
    Color_t color[18] = {1, 2, 4, 1, 1, 2, 2, 2, 4, 4, 4};
    Style_t marker[18] = {20, 21, 22, 20, 21, 22, 20, 21, 22};

    ifstream input;
    //char filename[1024];
    //char name[1024];
    string filename[]={
        "newbase-crosstalkstudyHV2100-A4-DN-D3.3output",
        "newbase-crosstalkstudyHV2100-A4-DN-D3.3-rm3Coutput",
        "HV2100-A4-DN-closeCoutput"
    };
    string legtitle[]={
        "Around","Farest","Nearest"
    };

/*
    string filename[]={
        "HV2100-A4-DN-close10nFoutput",
        "HV2100-A4-DN-closeCoutput",
        "HV2100-A4-DN-close470pFoutput"
    };
    string legtitle[]={
        "C=10nF","C=1nF","C=470pF"
    };
*/
    int fileN = sizeof(filename)/sizeof(filename[0]);
    for (int i = 0; i < fileN; i++)
    {
        sprintf(buff, "%s/%s.dat", path, filename[i].data());
        input.open(buff);
        if (!input)
        {
            cout << "file can't be found: \n"
                 << buff << endl;
            return 0;
        }
            cout << "enter the file: " << buff << endl;

        int k = 0;
        //return 1;
        while (input >> sigct >> sigring >> sigring2nd && !input.eof())
        {
            input >> ct[i][k] >> ring[i][k] >> ring2nd[i][k];
        Amp[i][k] = (cutstart + step * k) * 1e3;

            cout << sigct << "\t" << sigring << "\t" << sigring2nd << endl;
        cout << Amp[i][k] << endl;
            k++;
            //if(input.eof()) break;
        }
        pointN[i] = k;
        input.close();

    
    }

    TCanvas *c = cdC(0);
    DrawMyPad(gPad, "True Signal Amp (mV)", "Crosstalk Ratio ", (cutstart - step) * 1e3, (cutend + 4 * step) * 1e3, yl, yh);
    //DrawMyPad(gPad, xtitle, ytitle, x_low, x_high, y_low * 0.7, y_high * 1.1);
    //return 1;
    TLegend *leg = DrawMyLeg(0.6, 0.7, 0.95, 0.9, 62, 0.03);
    double fac = 0;
    for (int i = 0; i < fileN; i++)
    {

        g[i] = new TGraphErrors(pointN[i], Amp[i], ct[i], 0, 0);
        DrawMyGraph(g[i], "", "", 1.2, marker[i], color[i]);
        g[i]->Draw("sameP");
        //sprintf(buff, "ct%d/AMP%d", idindex[i], idindex[0]);
        sprintf(buff, "%s", legtitle[i].data());
        leg->AddEntry(g[i], buff, "lp");
    }
    leg->Draw();
    sprintf(buff, "%s/%sct.png", path, name);
    TGaxis::SetMaxDigits(3);
    gPad->SaveAs(buff);

    leg->Clear();
    c = cdC(1);
    DrawMyPad(gPad, "True Signal Amp (mV)", "Ringing Ratio ", (cutstart - step) * 1e3, (cutend + 4 * step) * 1e3, yl, yh);
    for (int i = 0; i < fileN; i++)
    {
        g[i] = new TGraphErrors(pointN[i], Amp[i], ring[i], 0, 0);
        DrawMyGraph(g[i], "", "", 1.2, marker[i], color[i]);
        g[i]->Draw("sameP");
        sprintf(buff, "%s", legtitle[i].data());
        leg->AddEntry(g[i], buff, "lp");
    }
    leg->Draw();
    sprintf(buff, "%s/%sring.png", path, name);
    gPad->SaveAs(buff);

    leg->Clear();
    c = cdC(2);
    DrawMyPad(gPad, "True Signal Amp (mV)", "2nd Ring Ratio ", (cutstart - step) * 1e3, (cutend + 4 * step) * 1e3, yl, yh);
    for (int i = 0; i < fileN; i++)
    {
        g[i] = new TGraphErrors(pointN[i], Amp[i], ring2nd[i], 0, 0);
        DrawMyGraph(g[i], "", "", 1.2, marker[i], color[i]);
        g[i]->Draw("sameP");
        sprintf(buff, "%s", legtitle[i].data());
        leg->AddEntry(g[i], buff, "lp");
    }
    leg->Draw();
    sprintf(buff, "%s/%sring2nd.png", path, name);
    gPad->SaveAs(buff);
    return 1;
}
