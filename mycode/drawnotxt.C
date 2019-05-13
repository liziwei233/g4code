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

#include "Include/DrawMyfunc.h"
using namespace std;

//int main()
int drawnotxt(const char* name="")
{
    
    gStyle->SetOptTitle(0);
  
    char str[1024];
    char buff[1024];
   
    

   
   


float x1[]={2000,2100,2200,2300}; //deep=width changed
float y1[]={1.35e5,2.98e5,8.07e5,2.47e6};
//float y2[]={2011,2011,2012,2011};
//float y3[]={148,152,108,92};


int n = sizeof(x1)/sizeof(x1[0]);

TGraph* g1 = new TGraph(n,x1,y1); 

/*
TH1F *h1 = new TH1F("h1","",100,-1,1);
h1->FillRandom("gaus",1e3);
TF1 *f1 = new TF1("f1","sin(x)",-1,1);
*/
TCanvas *c1 = new TCanvas("c1","",800,600);
c1->cd();
g1->Draw("AP");
DrawMyfunc mydraw;
mydraw.SetPad(c1,0.1,0.14,0.05,0.05);
//mydraw.Graph(g1,"NPE","TimeRes (ps)",1.5,20,4);

mydraw.Graph(g1,"Work voltage (V)","Gain ",1.5,20,4);
//g1->GetYaxis()->SetRangeUser(0,50);

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

return 0;
}

int main(){
    return drawnotxt();
} 