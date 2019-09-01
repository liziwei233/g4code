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
/*
Gian: 102083
rise time: 0.289988
transit time: 16.5259
time res: 0.13314
 */

#include <string>
#include <time.h>
#include <TString.h>

#include "Include/DrawMyfunc.h"
using namespace std;

//int main()
int drawSiPMTR(const char* name="/mnt/f/SiPM/labtest/Scintillator/8-13/8-13source_bar2_v2.1_re")
{
    
    gStyle->SetOptTitle(0);
    
    //**
    //** set your parameters **//
    //
    
    double peArray[]={60,90,120,150,180,210,240,270,300};

    //****************************
    //

    ifstream input;
    int i=0,j=0,k=0;
    char str[1024];
    char buff[1024];
    double chN,frac,iterN;
    
    int N = sizeof(peArray)/sizeof(peArray[0]);
	double tts[N];
    double ttsErr[N];
    double zero[N];

    for (int k=0; k<N; k++){
        sprintf(buff,"%s_TR_thred%gpe",name,peArray[k]);
        
    sprintf(buff,"%s.dat",buff);
    //sprintf(buff,"%s%d.dat",name,k);
    input.open(buff);
    
    

    if(!input){
		cout<<"file can't be found: \n"<<buff<<endl;
		return 0;
    }
    while(!input.eof()){
        /*for(i=0;i<charN;i++)
        {
            input>>str;
            cout<<str<<endl;
        }
        */
       cout<<"enter the file: "<<buff<<endl;
        //for( j=0;j<4;j++){
        //for( k=0;k<10;k++){
        
        input>>str>>chN>>frac>>iterN>>tts[k]>>ttsErr[k];
        cout<<tts[k]<<"\t"<<ttsErr[k]<<endl;
        //tr[j][k]=hitsigma*1e3;
        //}
        //}
        //if(j==4&&k==10) 
        break;
    }
    input.close();
    zero[k]=0;
    }

//float x1[]={0,0.1,6,12}; //deep=width changed
//float y1[]={3.9,3.7,3.8,3.6};
//float y2[]={2011,2011,2012,2011};
//float y3[]={148,152,108,92};

//int N=sizeof(sample)/sizeof(sample[0]);

TGraphErrors *g1 = new TGraphErrors(N,peArray,tts,zero,ttsErr);



/*
TH1F *h1 = new TH1F("h1","",100,-1,1);
h1->FillRandom("gaus",1e3);
TF1 *f1 = new TF1("f1","sin(x)",-1,1);
*/
TCanvas *c1 = new TCanvas("c1","",800,600);
c1->SetGrid();
c1->cd();
g1->Draw("AP");
DrawMyfunc mydraw;
mydraw.SetPad(c1,0.12,0.14,0.05,0.05);
//mydraw.Graph(g1,"NPE","TimeRes (ps)",1.5,20,4);

mydraw.Graph(g1,"Threshold (pe)","Time Res (ns)",1.5,20,4);
g1->GetYaxis()->SetRangeUser(0.15,0.35);

sprintf(buff,"%sTR.png",name);
c1->SaveAs(buff);
return 0;
}