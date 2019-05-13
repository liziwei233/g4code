//
//* Environment: 
//**     class LZWfunc
//**     class DrawMyfunc
//* Data: input manually.
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
   
//BJT_V2.1_boardS-ch1,M=50    
float x1={40,50,60,70,80,90,100,110,120,130,140,150,160,170}
float y1={4.95,4.7,4.17,3.66,3.29,2.98,2.67,2.07,1.34,0.77,0.47,0.33,0.27,0.22}   

//BJT_V2.1_boardL-ch4,M=50
float x1={40,50,60,70,80,90,100,110,120,130,140,150,160,170}
float y1={3.64,3.4,3.25,3.15,2.97,2.51,1.67,0.87,0.38,0.23,0.18,0.15,0.13,0.11}

//BJT_V3_boardL-ch1,M=50
float x1={10,13,15,18,20,23,25,28,30,33,35,38,40,45,50,55,60,65};
float y1={2.65,2.57,2.53,2.48,2.42,2.1,2.01,1.5,1.13,0.69,0.53,0.41,0.37,0.29,0.22,0.16,0.11,0.07};

//BJT_V3_BoardL-ch2,M=100
float x1={10,12,15,18,20,23,25,28,30,33,35,38,40,43,45};
float y1={3.46,3.26,2.95,2.63,2.43,2.00,1.61,0.99,0.68,0.43,0.36,0.27,0.23,0.18,0.16};

//BJT_V3_BoardS-ch2,M=100
float x1={15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100};
float y1={2.67,2.5,2.41,2.29,2.24,2.14,1.98,1.71,1.34,0.92,0.55,0.37,0.29,0.25,0.22,0.19,0.16,0.14} 


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