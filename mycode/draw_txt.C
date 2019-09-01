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
int draw_txt(int colNum, int Nevents)
{
    
    gStyle->SetOptTitle(0);
    
    //**
    //** set your parameters **//
    //
    int startpoint = 0;
    const int N=Nevents+startpoint;
    
    double gain[N];
    double risetime[N];
    double tt[N];
    double tts[N];
    double x[N];
    //****************************
    //

    ifstream input;
    int i=0,j=0,k=0;
    char str[1024];
    char buff[1024];
    //**
	//** G3 colum 1,1-7
	//const char* name="/mnt/f/MCP/7-26/KT0900_G3_Scanspot";
    //**
	//** G3 colum 2,1-7
	//const char* name="/mnt/f/MCP/7-27/KT0900_G3_scan/KT0900_G3_scanspot";
	
	//**
	//G4 colum 1 & colum 2,0-5
	//const char* name="/mnt/f/MCP/7-27/KT0900_G4_scan/KT0900_scanG4_";

    //**
	//** colum 3 , 1-7
	//**
	//const char* name="/mnt/f/MCP/7-28/col/colspot";

    //**
	//** anode , 1-16
	//**
	//const char* name="/mnt/f/MCP/7-28/scan/scanA";

    //**
	//** invert the tube, and scan column1
	//**
	//const char* name="/mnt/f/MCP/7-29/colum/columinvert_spot1_";

    //**
	//** change the base to version 5_a, and scan column1, 0-3 of group 1
	//**
	//const char* name="/mnt/f/MCP/7-30/base5_a/base5_aA";

    //**
	//** change the base to version 5_a, and scan column1, 0-3 of group 1
	//**
	const char* name="/mnt/f/MCP/7-30/KT0881_base5_b/KT0881_base5_bG4_";

    for (int k=startpoint; k<N; k++){
        x[k]=k*3;
    sprintf(buff,"%s%d_%d.dat",name,colNum,k*3);
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
        input>>gain[k]>>risetime[k]>>tt[k]>>tts[k];
        cout<<gain[k]<<"\t"<<risetime[k]<<"\t"<<tt[k]<<"\t"<<tts[k]<<endl;
        //tr[j][k]=hitsigma*1e3;
        //}
        //}
        //if(j==4&&k==10) 
        break;
    }
    input.close();

    }

//float x1[]={0,0.1,6,12}; //deep=width changed
//float y1[]={3.9,3.7,3.8,3.6};
//float y2[]={2011,2011,2012,2011};
//float y3[]={148,152,108,92};

//int N=sizeof(sample)/sizeof(sample[0]);


TGraph* g1 = new TGraph(N,x,gain); 
TGraph* g2 = new TGraph(N,x,risetime); 
TGraph* g3 = new TGraph(N,x,tt); 
TGraph* g4 = new TGraph(N,x,tts); 

/*
TH1F *h1 = new TH1F("h1","",100,-1,1);
h1->FillRandom("gaus",1e3);
TF1 *f1 = new TF1("f1","sin(x)",-1,1);
*/
TCanvas *c1 = new TCanvas("c1","",800,600);
TCanvas *c2 = new TCanvas("c2","",800,600);
TCanvas *c3 = new TCanvas("c3","",800,600);
TCanvas *c4 = new TCanvas("c4","",800,600);
c1->SetGrid();
c2->SetGrid();
c3->SetGrid();
c4->SetGrid();
c1->cd();
g1->Draw("AP");
DrawMyfunc mydraw;
mydraw.SetPad(c1,0.12,0.14,0.05,0.05);
//mydraw.Graph(g1,"NPE","TimeRes (ps)",1.5,20,4);

mydraw.Graph(g1,"position (1.5mm)","Gain",1.5,20,4);
g1->GetYaxis()->SetRangeUser(500,3.0e6);

sprintf(buff,"%scol%dgain.png",name,colNum);
c1->SaveAs(buff);

c2->cd();
g2->Draw("AP");
mydraw.SetPad(c2,0.12,0.14,0.05,0.05);
//mydraw.Graph(g1,"NPE","TimeRes (ps)",1.5,20,4);

mydraw.Graph(g2,"position (1.5mm)","Risetime  (ns)",1.5,20,4);
g2->GetYaxis()->SetRangeUser(0.05,0.7);

sprintf(buff,"%scol%drisetime.png",name,colNum);
c2->SaveAs(buff);

c3->cd();
g3->Draw("AP");
mydraw.SetPad(c3,0.12,0.14,0.05,0.05);
//mydraw.Graph(g1,"NPE","TimeRes (ps)",1.5,20,4);

mydraw.Graph(g3,"position (1.5mm)","Transit time (ns)",1.5,20,4);
g3->GetYaxis()->SetRangeUser(15,17);

sprintf(buff,"%scol%d_transittime.png",name,colNum);
c3->SaveAs(buff);

c4->cd();
g4->Draw("AP");
mydraw.SetPad(c4,0.12,0.14,0.05,0.05);
//mydraw.Graph(g1,"NPE","TimeRes (ps)",1.5,20,4);

mydraw.Graph(g4,"position (1.5mm)","TTS (ns)",1.5,20,4);
g4->GetYaxis()->SetRangeUser(0.01,0.2);

sprintf(buff,"%scol%d_tts.png",name,colNum);
c4->SaveAs(buff);



//draw->Hist(h1,"x axis","y axis",3,2,4,2);

//TLegend *leg=mydraw.Leg(0.6,0.7,0.85,0.9);
//DrawMyfunc::SetPad(c1,0.1,0.14,0.05,0.05);
/* 
for(i=0;i<4;i++){
if(i==0) g[i]->Draw("AP");
else g[i]->Draw("Psame");
sprintf(buff,"smear=%gps",sigma[3-i]*1e3);
leg->AddEntry(g[3-i],buff,"lp");
}
*/
//h1->Draw();
//leg->Draw();
//TLatex *l=draw->Latex("this is a text");
//l->Draw();
return 0;
}
/* 
int main(){
    return draw_txt(int colNum, int Nevents);
} 
*/