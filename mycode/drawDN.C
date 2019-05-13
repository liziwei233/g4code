
//* Environment: 
//**     class LZWfunc
//**     class DrawMyfunc
//* Data: ASCII file.
//* Function: 
//**      draw figure of darknoise .
//* Date: 2019.3.29
//*
//***********************************
//

#include <string>
#include <time.h>
#include <TString.h>

#include "Include/DrawMyfunc.h"
using namespace std;

//int main()
int drawDN(const char* name="")
{
    
    gStyle->SetOptTitle(0);
    ifstream input;
    char buff[1024];
    sprintf(buff,"rawdata/darknoise.txt");
    input.open(buff);
    

    double x[21],y[21];
    double th=0,dn=0;
    int i=0;
    if(!input){
		cout<<"file can't be found!"<<endl;
		return 0;
    }
    while(!input.eof()){
        /*for(i=0;i<charN;i++)
        {
            input>>str;
            cout<<str<<endl;
        }
        */
       cout<<"enter the file: "<<endl;
        input>>th>>dn;
        cout<<th<<"\t"<<dn<<endl;
        x[i]=th;
        y[i]=dn/10.;
        i++;
        if(i==21) break;
    }
    input.close();


//float x1[]={0,0.1,6,12}; //deep=width changed
//float y1[]={3.9,3.7,3.8,3.6};
//float y2[]={2011,2011,2012,2011};
//float y3[]={148,152,108,92};

int n=sizeof(x)/sizeof(x[0]);

TGraph* g1 = new TGraph(n,x,y); 

/*
TH1F *h1 = new TH1F("h1","",100,-1,1);
h1->FillRandom("gaus",1e3);
TF1 *f1 = new TF1("f1","sin(x)",-1,1);
*/
TCanvas *c1 = new TCanvas("c1","",800,600);
c1->cd();
g1->Draw("AP");
DrawMyfunc mydraw;
mydraw.SetPad(c1,0.12,0.14,0.08,0.08);
//mydraw.Graph(g1,"NPE","TimeRes (ps)",1.5,20,4);

mydraw.Graph(g1,"Threshold (mV)","Dark Noise (s^{-1})",1.5,20,1);
TGaxis::SetMaxDigits(3);
g1->GetYaxis()->SetRangeUser(0,2500);
//draw->Hist(h1,"x axis","y axis",3,2,4,2);

//DrawMyfunc::SetPad(c1,0.1,0.14,0.05,0.05);
//h1->Draw();
//TLatex *l=draw->Latex("this is a text");
//l->Draw();
sprintf(buff,"%s_dk.png",name);
c1->SaveAs(buff);

return 0;
}

int main(){
    return drawDN();
} 