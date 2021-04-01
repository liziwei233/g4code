//
//* Environment: 
//**     class LZWfunc
//**     class DrawMyfunc
//* Data: ASCII file.
//* Function: 
//**      draw figure simply.
//* Date: 2019.3.19
//*
//***********************************
//

#include <string>
#include <time.h>
#include <TString.h>

#include "Include/DrawMyClass.h"
using namespace std;

double HVfun(double *x, double *par){
    double val=0.;
    double A=par[1];
    double alpha=par[2];
    double C=par[0];
    //val = A*TMath::Power(x[0],beta);
    val = TMath::Exp(C+x[0]*A*alpha);
    return val;
}
double HVfun_complex(double *x, double *par)
{
    double val = 0.;
    double E0 = par[0]; //~1eV, initial energy of the secondary electron
    double alpha = par[1]; //=40, the ratio of length to diameter
    double K = par[2]; //~0.2, a constant from the ralation delta = K*Ec ,Ec is the collision energy
    //double C = par[3];
    //val = A*TMath::Power(x[0],beta);
    double a =4.*E0*alpha*alpha;
    val = TMath::Power(K*x[0]*x[0]/a,a/x[0]);
    return val;
}
double HVfun_complex2(double *x, double *par)
{
    double val = 0.;
    double E0 = par[0]; //~1eV, initial energy of the secondary electron
    double alpha = par[1]; //=40, the ratio of length to diameter
    double A = par[2]; //~0.2, is the propotionality constant in the relation delta = A*sqrt(Ec);
    //val = A*TMath::Power(x[0],beta);
    double a =2.*TMath::Sqrt(E0)*alpha;
    val = TMath::Power(A*x[0]/a,a*a/x[0]);
    return val;
}

//int main()
void drawHV(const char* name="")
{
    
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(1112);
    char str[1024];
    char buff[1024];
   
    

   
   


float x1[]={3100,3150,3200,3250,3300,3350}; //deep=width changed
float y1[]={2.65e5,3.01e5,5.43e5,6.7e5,8e5,1.23e6};
//float y2[]={2011,2011,2012,2011};
//float y3[]={148,152,108,92};


int n = sizeof(x1)/sizeof(x1[0]);

TGraph* g1 = new TGraph(n,x1,y1); 

/*
TH1F *h1 = new TH1F("h1","",100,-1,1);
h1->FillRandom("gaus",1e3);
TF1 *f1 = new TF1("f1","sin(x)",-1,1);
*/
TCanvas *c1;
c1 = cdC(1);
c1->SetLogy();
g1->Draw("AP");
//mydraw.Graph(g1,"NPE","TimeRes (ps)",1.5,20,4);

DrawMyGraph(g1,"Work voltage (kV)","Gain ",1.5,20,4);

/*
TF1* fhv= new TF1("fhv",HVfun_complex,3000,3400,3);
fhv->SetParNames("E0","#alpha","K");
//fhv->SetParLimits(0,1,1e7);
//fhv->SetParLimits(1,1,10);
fhv->SetParameter(0,1);
fhv->FixParameter(1,40);
fhv->SetParameter(2,0.2);
//fhv->SetParameter(3,0);
*/

/*
TF1* fhv= new TF1("fhv",HVfun_complex2,3000,3400,3);
fhv->SetParNames("E0","#alpha","A");
//fhv->SetParLimits(0,1,1e7);
//fhv->SetParLimits(1,1,10);
fhv->SetParameter(0,1);
fhv->FixParameter(1,40);
fhv->SetParameter(2,0.2);
*/

/*
TF1* fhv= new TF1("fhv",HVfun,3000,3400,3);
fhv->SetParNames("C","A","#alpha");
//fhv->SetParLimits(0,1,1e7);
//fhv->SetParLimits(1,1,10);
//fhv->SetParameter(0,-10);
fhv->FixParameter(2,40);
//fhv->SetParameter(1,0.2);
*/
g1->Fit(fhv);
g1->GetXaxis()->SetRangeUser(1800,2400);
g1->GetYaxis()->SetRangeUser(6e4,1e7);

/*
mydraw.Graph(g[1],"","",1.2,23,2);
mydraw.Graph(g[2],"","",1.2,21,1);
mydraw.Graph(g[3],"","",1.8,29,8);
//draw->Hist(h1,"x axis","y axis",3,2,4,2);

TLegend *leg=mydraw.Leg(0.6,0.7,0.85,0.9);
//DrawMyfunc::SetPad(c1,0.1,0.14,0.05,0.05);
for(i=0;i<4;i++){
if(i==0) g[i]->Draw("AP");
else g[i]->Draw("Psame");
sprintf(buff,"smear=%gps",sigma[3-i]*1e3);
leg->AddEntry(g[3-i],buff,"lp");
}
//h1->Draw();
leg->Draw();
//TLatex *l=draw->Latex("this is a text");
//l->Draw();
*/
sprintf(buff,"%s_HVscan.png",name);
c1->SaveAs(buff);

}
