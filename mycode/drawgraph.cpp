#include "Include/DrawMyClass.h"
TGraphErrors *g[10];
TGraphErrors* graph(double* x,double* y,double* xerr,double* yerr,const int n=7,const char* xtitle="Position (mm)",const char* ytitle="Gain (#times10^{7})",const int id=0)
{
    
    TCanvas *c = cdC(id);
    double xrange=TMath::MaxElement(n,x)-TMath::MinElement(n,x);
    double x_low = TMath::MinElement(n,x)-(xrange)/n; 
    double x_high = TMath::MaxElement(n,x)+(xrange)/n;
    double yrange=TMath::MaxElement(n,y)-TMath::MinElement(n,y);
    double y_low = TMath::MinElement(n,y)-(yrange)/n;
    double y_high = TMath::MaxElement(n,y)+(yrange)/n;
    cout<<x_low<<"\t"<<x_high<<"\t"<<y_low<<"\t"<<y_high<<endl;
    
    DrawMyPad(gPad,xtitle,ytitle,x_low,x_high,y_low*0.7,y_high*1.1);
    
    g[id] = new TGraphErrors(n,x,y,xerr,yerr);
    DrawMyGraph(g[id],(char*)xtitle,(char*)ytitle,1.2,20,1,1,2);
    g[id]->Draw("sameP");
    //g[id]->GetXaxis()->SetRangeUser(x[0]-1,x[n-1]);
    return g[id];
}

void drawgraph()
{
    setgStyle();
    double x[]={0,0.3,0.8,1.8,3.8,5.8,6.8};
    double yGain[]={2.33,2.542,2.624,2.604,2.513,2.057,1.334};
    double yGerr[]={0.02,0.017,0.018,0.02,0.021,0.021,0.02};
    double yPurity[]={0.9817,0.9866,0.9828,0.9786,0.9766,0.8179,0.5773};
    double yTTS[]={64,62,63,64,68,69,74};
    double yTTSerr[]={1.6,1.3,1.5,1.7,2.0,2.1,3.5};

    double y1[]={1.43,1.56,1.61,1.60,1.55,1.27,0.81};
    double y2[]={64,62,63,64,68,69,74};
    double y3[]={344,343,338,335,338,348,360};
    const int n=sizeof(y1)/sizeof(y1[0]);
    double xerr[n]={0}; // 0.5mm/sqrt(12)
    double yerr[n]={0};

    double GainSyserr = 0.0257; //unit:pC
    double GainSysTTS = 1.89; //ps;

    for(int i =0; i<n; i++)
    {
        yGain[i] = yGain[i]/1.6;
        yGerr[i] = yGerr[i]/yPurity[i];
        yGerr[i] = sqrt(yGerr[i]*yGerr[i]+0.0257*0.0257)/1.6;
        xerr[i]=0.5/sqrt(12);

        yTTSerr[i] = yTTSerr[i]/yPurity[i];
        yTTSerr[i] = sqrt(yTTSerr[i]*yTTSerr[i]+GainSysTTS*GainSysTTS);

    }

    g[0]=graph(x,yGain,xerr,yGerr,n,"Position (mm)","Gain (#times10^{7})",0);
    //g[0]=graph(x,y,xerr,yerr,n,"Position (mm)","Gain (#times10^{7})",0);
    gPad->SaveAs("Gain.png");
    gPad->SaveAs("Gain.pdf");
    g[1]=graph(x,yTTS,xerr,yTTSerr,n,"Position (mm)","TTS (ps)",1);
    //g[1]=graph(x,y2,xerr,yerr,n,"Position (mm)","TTS (ps)",1);
    gPad->SaveAs("TTS.png");
    gPad->SaveAs("TTS.pdf");
    //g[2]=graph(x,y3,xerr,yerr,n,"Position (mm)","Risetime (ps)",2);
    //gPad->SaveAs("risetime.png");
}