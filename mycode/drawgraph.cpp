#include "Include/DrawMyClass.h"
TGraphErrors *g[10];
TGraphErrors* graph(double* x,double* y,double* xerr,double* yerr,const int n=7,const char* xtitle="Position (mm)",const char* ytitle="Gain (#times10^{7})",const int id=0, bool logy=0)
{

    TCanvas *c = cdC(id);
    if(logy) c->SetLogy();
    double xrange=TMath::MaxElement(n,x)-TMath::MinElement(n,x);
    double xstep = xrange/n;
    double x_low = TMath::MinElement(n,x)-(xrange)/n; 
    double x_high = TMath::MaxElement(n,x)+(xrange)/n;
    if(TMath::RMS(n,x)/TMath::Mean(n,x)>10)  {
        x_low = TMath::MinElement(n,x)*0.2;
        x_high = TMath::MaxElement(n,x)*5;
    }
    double yrange=TMath::MaxElement(n,y)-TMath::MinElement(n,y);
    double y_low = TMath::MinElement(n,y)-(yrange)/n;
    double y_high = TMath::MaxElement(n,y)+(yrange)/n;
    if(TMath::RMS(n,y)/TMath::Mean(n,y)>0.7)  {
        y_low = TMath::MinElement(n,y)*0.2;
        y_high = TMath::MaxElement(n,y)*5;
    }
    cout<<"sigma/mean: "<<TMath::RMS(n,y)/TMath::Mean(n,y)<<endl;
    cout<<x_low<<"\t"<<x_high<<"\t"<<y_low<<"\t"<<y_high<<endl;
    
    //DrawMyPad(gPad,xtitle,ytitle,x_low,x_high,y_low*0.7,y_high*1.1);
    DrawMyPad(gPad,xtitle,ytitle,x_low,x_high,y_low*0.7,y_high*1.3);
    
    g[id] = new TGraphErrors(n,x,y,xerr,yerr);
    DrawMyGraph(g[id],(char*)xtitle,(char*)ytitle,1.2,20,kGreen+2,kGreen+2,2);
    g[id]->Draw("sameP");
    //g[id]->GetXaxis()->SetRangeUser(x[0]-1,x[n-1]);
    return g[id];
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
    double TTSSyserr = 1.89; //ps;

    for (int i = 0; i < n; i++)
    {
        yGain[i] = yGain[i]/1.6;
        yGerr[i] = yGerr[i]/yPurity[i];
        yGerr[i] = sqrt(yGerr[i]*yGerr[i]+0.0257*0.0257)/1.6;
        xerr[i]=0.5/sqrt(12);

        yTTSerr[i] = yTTSerr[i]/yPurity[i];
        yTTSerr[i] = sqrt(yTTSerr[i]*yTTSerr[i]+TTSSyserr*TTSSyserr);

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

void drawgraph1(const char* name="UV-2")
{
    setgStyle();

    //
    // ** Parameters of UV-6 
    // *
 /*   
    double x[] = {2700,2600,2500,2450,2400,2350,2300,2250,2200};
    double yTTS[]={44,  50, 58, 62, 69,75,80,87,94};
    double yTTSerr[] = {1.1,    1.1,    0.9,    1.0,    1.4, 1.4,1.6,1.8,2.4};
    double yRise[] = {344,  341,    314,    313,    309,    302,304,299,296};
    double yRiseerr[] = {0, 0.6,    0.5,    0.6,    0.8,    0.8,1.1,1.2,1.4};
    double yGain[] = {8.34e5,5.25e5,3.44e5,2.43e5,1.93e5,1.51e5,8.39e4,4.89e4,4.24e4};
    //double yGainerr[] = {0.0017, 0.0014, 0.0005, 0.0005,0.0005,0.0007,0.0002,0.0001,0.0001};
    // 
  */  
    /*
    //
    // ** Parameters of UV-2 
    // *
    double x[] = {2600,2500,2400,2350,2300,2250,2200,2150};
    double yGain[] = {1.17e6,8.2e5,5.25e5,4.5e5,3.58e5,2.35e5,1.79e5,1.37e4};
    double yTTS[]={47,46,52,58,64,79,84,93};
    double yTTSerr[] = {1.9,0.8,0.9,1.5,1.6,2.3,2.9,3.4};
    double yRise[] = {258,342,337,338,333,328,321,315};
    double yRiseerr[] = {1.3,0.3,0.4,0.9,1.2,1,1.8,1.9};
    //
*/
    //
    // ** Parameters of Bi-2 
    // *
    double x[] = {2600,2500,2400,2300,2200};
    double yGain[] = {3.5e6,2.1e6,1.26e6,7.3e5,3.74e5};
    double yTTS[]={52,61,70,71,105};
    double yTTSerr[] = {1.0,1.5,3.5,1.9,5.5};
    double yRise[] = {345,375,372,365,347};
    double yRiseerr[] = {0.8,1.9,0.9,1.3,2};
    //

    const int n=sizeof(x)/sizeof(x[0]);
    double xerr[n]={0}; 
    double yGainerr[n] = {0};

    //double GainSyserr = 0.0257; //unit:pC
    double TTSSyserr = 1.89; //ps;
    double GainSyserr = 0.0; //unit:pC
    for(int i =0; i<n; i++)
    {

        yTTSerr[i] = sqrt(yTTSerr[i]*yTTSerr[i]+TTSSyserr*TTSSyserr);
        yGainerr[i] = sqrt(yGainerr[i]*yGainerr[i]+GainSyserr*GainSyserr);
        yGainerr[i] = yGainerr[i]*1.e-12/1.6e-19;

    }

    g[0]=graph(x,yRise,xerr,yRiseerr,n,"HV (V)","Risetime (ps)",0);
    //g[0]=graph(x,y,xerr,yerr,n,"Position (mm)","Gain (#times10^{7})",0);
    sprintf(buff, "%s-Risetime.png",name);
    gPad->SaveAs(buff);
    g[1]=graph(x,yTTS,xerr,yTTSerr,n,"HV (V)","TTS (ps)",1);
    //g[1]=graph(x,y2,xerr,yerr,n,"Position (mm)","TTS (ps)",1);
    sprintf(buff, "%s-TTS.png",name);
    gPad->SaveAs(buff);
    //g[2]=graph(x,y3,xerr,yerr,n,"Position (mm)","Risetime (ps)",2);
    //gPad->SaveAs("risetime.png");

    g[2]=graph(x,yGain,xerr,yGainerr,n,"HV (V)","Gain",2,1);
    TF1 *fhv = new TF1("fhv", HVfun, 2100, 2700, 3);
    fhv->SetParNames("cons", "#delta", "#alpha");
    //fhv->SetParLimits(0,1,1e7);
    //fhv->SetParLimits(1,1,10);
    fhv->SetParameter(0, -10);
    fhv->FixParameter(2, 40);
    g[2]->Fit(fhv,"","",2200,2800);
    fhv->Draw("same");
    //gPad->SetLogy();
    //g1->GetXaxis()->SetRangeUser(1800, 2900);
    //g1->GetYaxis()->SetRangeUser(1e4, 1e7);
    //gPad->Update();
    //gPad->Modified();
    //g[1]=graph(x,y2,xerr,yerr,n,"Position (mm)","TTS (ps)",1);
    sprintf(buff, "%s-Gain.png",name);
    gPad->SaveAs(buff);
}

void drawgraph2(const char* name="MCP1-MCP2")
{
    setgStyle();

 /*   //
    // ** Parameters of Bi-2 MCP2-GND
    // *
    double x[] = {100,200,300};
    double yGain[] = {4.42e6,3.93e6,4.88e6};
    double yTTS[]={57,52,45};
    double yTTSerr[] = {2,3.6,1.7};
    double yRise[] = {537,538,537};
    double yRiseerr[] = {0,0,0};
    //
*/
/*
    //
    // ** Parameters of Bi-2 PK-MCP1
    // *
    double x[] = {100,200,300};
    double yGain[] = {2.86e6,3.8e6,3.77e6};
    double yTTS[]={50,52,53};
    double yTTSerr[] = {1.9,3.6,2.6};
    double yRise[] = {549,538,536};
    double yRiseerr[] = {0,0,0};
    //
*/
    //
    // ** Parameters of Bi-2 MCP1-MCP2
    // *
    double x[] = {1800,1900,2000,2150};
    double yGain[] = {2.94e5,9.38e5,1.36e6,3.77e6};
    double yTTS[]={137,64,53,53};
    double yTTSerr[] = {6,4.6,2.2,2.6};
    double yRise[] = {624,601,564,536};
    double yRiseerr[] = {0,0,0,0};
    //

    const int n=sizeof(x)/sizeof(x[0]);
    double xerr[n]={0}; 
    double yGainerr[n] = {0};

    //double GainSyserr = 0.0257; //unit:pC
    double TTSSyserr = 1.89; //ps;
    double GainSyserr = 0.0; //unit:pC
    for(int i =0; i<n; i++)
    {

        yTTSerr[i] = sqrt(yTTSerr[i]*yTTSerr[i]+TTSSyserr*TTSSyserr);
        yGainerr[i] = sqrt(yGainerr[i]*yGainerr[i]+GainSyserr*GainSyserr);
        yGainerr[i] = yGainerr[i]*1.e-12/1.6e-19;

    }

    g[0]=graph(x,yRise,xerr,yRiseerr,n,"Voltage of PK-MCP1 (V)","Risetime (ps)",0);
    //g[0]=graph(x,y,xerr,yerr,n,"Position (mm)","Gain (#times10^{7})",0);
    sprintf(buff, "%s-Risetime.png",name);
    gPad->SaveAs(buff);
    g[1]=graph(x,yTTS,xerr,yTTSerr,n,"Voltage of MCP1-MCP2 (V)","TTS (ps)",1);
    //g[1]=graph(x,y2,xerr,yerr,n,"Position (mm)","TTS (ps)",1);
    sprintf(buff, "%s-TTS.png",name);
    gPad->SaveAs(buff);
    //g[2]=graph(x,y3,xerr,yerr,n,"Position (mm)","Risetime (ps)",2);
    //gPad->SaveAs("risetime.png");

    g[2]=graph(x,yGain,xerr,yGainerr,n,"Voltage of PK-MCP1 (V)","Gain",2,1);
    //TGaxis::SetMaxDigits(3);
    /*
    TF1 *fhv = new TF1("fhv", HVfun, 2100, 2700, 3);
    fhv->SetParNames("cons", "#delta", "#alpha");
    //fhv->SetParLimits(0,1,1e7);
    //fhv->SetParLimits(1,1,10);
    fhv->SetParameter(0, -10);
    fhv->FixParameter(2, 40);
    g[2]->Fit(fhv,"","",2200,2800);
    fhv->Draw("same");
    */
    //gPad->SetLogy();
    //g1->GetXaxis()->SetRangeUser(1800, 2900);
    //g1->GetYaxis()->SetRangeUser(1e4, 1e7);
    //gPad->Update();
    //gPad->Modified();
    //g[1]=graph(x,y2,xerr,yerr,n,"Position (mm)","TTS (ps)",1);
    sprintf(buff, "%s-Gain.png",name);
    gPad->SaveAs(buff);
}
