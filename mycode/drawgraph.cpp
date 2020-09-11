#include "Include/DrawMyClass.h"
TGraphErrors *g[10];
TGraphErrors *graph(double *x, double *y, double *xerr, double *yerr, const int n = 7, const char *xtitle = "Position (mm)", const char *ytitle = "Gain (#times10^{7})", const int id = 0)
{

    TCanvas *c = cdC(id);
    double xrange = TMath::MaxElement(n, x) - TMath::MinElement(n, x);
    double x_low = TMath::MinElement(n, x) - (xrange) / n;
    double x_high = TMath::MaxElement(n, x) + (xrange) / n;
    double yrange = TMath::MaxElement(n, y) - TMath::MinElement(n, y);
    double y_low = TMath::MinElement(n, y) - (yrange) / n;
    double y_high = TMath::MaxElement(n, y) + (yrange) / n;
    cout << x_low << "\t" << x_high << "\t" << y_low << "\t" << y_high << endl;

    DrawMyPad(gPad, xtitle, ytitle, x_low, x_high, y_low * 0.7, y_high * 1.1);

    g[id] = new TGraphErrors(n, x, y, xerr, yerr);
    DrawMyGraph(g[id], (char *)xtitle, (char *)ytitle, 1.2, 20, 1, 1, 2);
    g[id]->Draw("sameP");
    //g[id]->GetXaxis()->SetRangeUser(x[0]-1,x[n-1]);
    return g[id];
}

void drawgraph()
{
    setgStyle();
    double x[] = {0, 0.3, 0.8, 1.8, 3.8, 5.8, 6.8};
    double yGain[] = {2.33, 2.542, 2.624, 2.604, 2.513, 2.057, 1.334};
    double yGerr[] = {0.02, 0.017, 0.018, 0.02, 0.021, 0.021, 0.02};
    double yPurity[] = {0.9817, 0.9866, 0.9828, 0.9786, 0.9766, 0.8179, 0.5773};
    double yTTS[] = {64, 62, 63, 64, 68, 69, 74};
    double yTTSerr[] = {1.6, 1.3, 1.5, 1.7, 2.0, 2.1, 3.5};

    double y1[] = {1.43, 1.56, 1.61, 1.60, 1.55, 1.27, 0.81};
    double y2[] = {64, 62, 63, 64, 68, 69, 74};
    double y3[] = {344, 343, 338, 335, 338, 348, 360};
    const int n = sizeof(y1) / sizeof(y1[0]);
    double xerr[n] = {0}; // 0.5mm/sqrt(12)
    double yerr[n] = {0};

    double GainSyserr = 0.0257; //unit:pC
    double GainSysTTS = 1.89;   //ps;

    for (int i = 0; i < n; i++)
    {
        yGain[i] = yGain[i] / 1.6;
        yGerr[i] = yGerr[i] / yPurity[i];
        yGerr[i] = sqrt(yGerr[i] * yGerr[i] + 0.0257 * 0.0257) / 1.6;
        xerr[i] = 0.5 / sqrt(12);

        yTTSerr[i] = yTTSerr[i] / yPurity[i];
        yTTSerr[i] = sqrt(yTTSerr[i] * yTTSerr[i] + GainSysTTS * GainSysTTS);
    }

    g[0] = graph(x, yGain, xerr, yGerr, n, "Position (mm)", "Gain (#times10^{7})", 0);
    //g[0]=graph(x,y,xerr,yerr,n,"Position (mm)","Gain (#times10^{7})",0);
    gPad->SaveAs("Gain.png");
    gPad->SaveAs("Gain.pdf");
    g[1] = graph(x, yTTS, xerr, yTTSerr, n, "Position (mm)", "TTS (ps)", 1);
    //g[1]=graph(x,y2,xerr,yerr,n,"Position (mm)","TTS (ps)",1);
    gPad->SaveAs("TTS.png");
    gPad->SaveAs("TTS.pdf");
    //g[2]=graph(x,y3,xerr,yerr,n,"Position (mm)","Risetime (ps)",2);
    //gPad->SaveAs("risetime.png");
}
char path[1024] = "/mnt/f/R10754/MPE-crosstalk/twotubes";
int draw_txt()
{
    setgStyle();
    const int N = 100;
    double A1[6][N];
    double A2[6][N];
    double A3[6][N];
    double A4[6][N];

    double ct1[6][N]; // from the same tube
    double ct2[6][N]; // from the neighbour tube
    double xerr[6][N] = {0};
    double yerr[6][N] = {0};

    //string name[6] = {"KT0881_2100V_2mmgap", "KT0881_2100V_2mmgap_mental", "KT0881_2100V_10mmgap_plastic", "KT0881_2100V_5mmgap_plastic", "KT0881_2100V_5mmgap_mental", "KT0881_2100V_5mmgap_mental2"} ;
    string name[5] = {"KT0881_2100V_2mmgap",
    "KT0881_2100V_2mmgap_mental",  "KT0881_2100V_10mmgap_plastic", "KT0881_2100V_5mmgap_plastic", "KT0881_2100V_5mmgap_mental"} ;
    
    string info[5] = {"pitch=2mm", "pitch=2mm,with mental(No ground)", "pitch=10mm", "pitch=5mm", "pitch=5mm, with mental"} ;
    Size_t size[6] = {2,1.7,2.1,2.1,1.8,2};
    Color_t color[6] = {1,2,4,kOrange-1,6,8};
    Style_t marker[6] = {20,21,22,23,25,26};
    ifstream input;
    char filename[1024];
    char buff[1024];
    //char name[1024];
    for (int s = 0; s < 5; s++)
    {

        sprintf(buff, "%s/%s.dat", path, name[s].c_str());
        input.open(buff);
        if (!input)
        {
            cout << "file can't be found: \n"
                 << buff << endl;
            return 0;
        }
        long filelen;
        input.seekg(0, ios::end);
        filelen = input.tellg();
        cout << "filelength: " << filelen << endl;
        input.seekg(0, ios::beg);
        int k = 0;
        //return 1;
        while (input >> A1[s][k] >> A2[s][k] >> A3[s][k] >> A4[s][k] && !input.eof())
        {

            cout << "enter the file: " << buff << endl;

            cout << A1[s][k] << "\t" << A2[s][k] << "\t" << A3[s][k] << "\t" << A4[s][k] << endl;
            k++;
            //if(input.eof()) break;
        }
        input.close();

    int n = k;

    for (int i = 0; i < n; i++)
    {
        ct1[s][i] = A1[s][i] / A2[s][i] * 100;
        ct2[s][i] = A3[s][i] / A2[s][i] * 100;
    }
     g[s] = new TGraphErrors(n, A2[s],  ct2[s], xerr[s], yerr[s]);
    DrawMyGraph(g[s], "", "", 1.2, marker[s], color[s]);
    }
    TCanvas *c = cdC(0);
    DrawMyPad(gPad, "Signal Amp (mV)", "Crosstalk Ratio (%) ", -5, 295, -0.5, 10.1);
    //DrawMyPad(gPad, xtitle, ytitle, x_low, x_high, y_low * 0.7, y_high * 1.1);
    //return 1;
    TLegend * leg = DrawMyLeg(0.5,0.6,0.95,0.9,62,0.03);
    for(int i= 0; i<5; i++)
    {
    g[i]->Draw("sameP");
    leg->AddEntry(g[i],info[i].c_str(),"lp");
    }
    leg->Draw();
    sprintf(buff, "%s/ct2.png", path);
    gPad->SaveAs(buff);
    return 1;
}
void op(const char *name = "KT0881_2100V_2mmgap_2")
{
    char filename[1024];
    char buff[1024];
    sprintf(buff, "%s/%s.txt", path, name);
    ofstream output(buff);
    output << 1 << "\t" << 2;
    //output<<1<<"\t"<<2<<"\t"<<3<<"\t"<<4<<endl;
}