#include "DrawMyfunc.h"

int main()
//void draw()
{
    
float x1[]={0,0.1,6,12}; //deep=width changed
float y1[]={3.9,3.7,3.8,3.6};
float y2[]={2011,2011,2012,2011};
float y3[]={148,152,108,92};

int n=sizeof(x1)/sizeof((x1)[0]);

TGraph *g1 = new TGraph(n,x1,y1); 
TGraph *g2 = new TGraph(n,x1,y2); 
TGraph *g3 = new TGraph(n,x1,y3);

TH1F *h1 = new TH1F("h1","",100,-1,1);
TF1 *f1 = new TF1("f1","sin(x)",-1,1);
TCanvas *c1 = new TCanvas("c1","",800,600);
c1->cd();
DrawMyfunc* draw;
draw->Hist(h1,"x axis","y axis",3,2,4,2);
h1->Draw();
c1->SaveAs("c1.png");

return 0;
} 