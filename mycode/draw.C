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
int draw(const char* name="")
{
    
    gStyle->SetOptTitle(0);
    double sample[10]={5,10,15,20,25,30,35,40,45,50};
    double sigma[4]={0,0.050,0.075,0.100};
    ifstream input;
    int charN=7;
    int i=0,j=0,k=0;
    char str[1024];
    char buff[1024];
    sprintf(buff,"%s_results.dat",name);
    input.open(buff);
    double smearsigma;
    double sNum;
    double hitmean;
    double hitsigma;
    

    double tr[4][10]={0};
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
        for( j=0;j<4;j++){
        for( k=0;k<10;k++){
        input>>smearsigma>>sNum>>hitmean>>hitsigma;
        cout<<smearsigma<<"\t"<<sNum<<"\t"<<hitmean<<"\t"<<hitsigma<<endl;
        tr[j][k]=hitsigma*1e3;
        }
        }
        if(j==4&&k==10) break;
    }
    input.close();


//float x1[]={0,0.1,6,12}; //deep=width changed
//float y1[]={3.9,3.7,3.8,3.6};
//float y2[]={2011,2011,2012,2011};
//float y3[]={148,152,108,92};

int n=sizeof(sample)/sizeof(sample[0]);

TGraph* g[4];
for(i=0;i<4;i++){
g[i] = new TGraph(n,sample,tr[i]); 
}
TGraph* g1 = new TGraph(n,sample,tr[1]); 

/*
TH1F *h1 = new TH1F("h1","",100,-1,1);
h1->FillRandom("gaus",1e3);
TF1 *f1 = new TF1("f1","sin(x)",-1,1);
*/
TCanvas *c1 = new TCanvas("c1","",800,600);
c1->cd();
g[0]->Draw("AP");
DrawMyfunc mydraw;
mydraw.SetPad(c1,0.1,0.14,0.05,0.05);
//mydraw.Graph(g1,"NPE","TimeRes (ps)",1.5,20,4);

mydraw.Graph(g[0],"NPE","TimeRes (ps)",1.5,20,4);
g[0]->GetYaxis()->SetRangeUser(0,50);
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
sprintf(buff,"%s_results.png",name);
c1->SaveAs(buff);

return 0;
}

int main(){
    return draw();
} 