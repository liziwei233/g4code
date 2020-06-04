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
    
    DrawMyPad(gPad,xtitle,ytitle,x_low,x_high,y_low,y_high);
    
    g[id] = new TGraphErrors(n,x,y,xerr,yerr);
    DrawMyGraph(g[id],(char*)xtitle,(char*)ytitle,2,20,1);
    g[id]->Draw("sameLP");
    //g[id]->GetXaxis()->SetRangeUser(x[0]-1,x[n-1]);
    return g[id];
}

void drawgraph()
{
    setgStyle();
    double x[]={0,0.3,0.8,1.8,3.8,5.8,6.8};
    double y1[]={1.43,1.56,1.61,1.60,1.55,1.27,0.81};
    double y2[]={64,62,63,64,68,69,74};
    double y3[]={344,343,338,335,338,348,360};
    const int n=sizeof(y1)/sizeof(y1[0]);
    double xerr[n]={0};
    double yerr[n]={0};

    g[0]=graph(x,y1,xerr,yerr,n,"Position (mm)","Gain (#times10^{7})",0);
    gPad->SaveAs("Gain.png");
    g[1]=graph(x,y2,xerr,yerr,n,"Position (mm)","TTS (ps)",1);
    gPad->SaveAs("TTS.png");
    g[2]=graph(x,y3,xerr,yerr,n,"Position (mm)","Risetime (ps)",1);
    gPad->SaveAs("risetime.png");
}