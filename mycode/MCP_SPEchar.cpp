//
//* Environment: 
//**     class LZWfunc
//**     class DrawMyfunc
//* Data: 
//**     .root file obtaied by .trc files of Oscilloscope.
//* Function:
//**     Draw Single photonelectron characteres of the test MCP.
//* Date: 2019.3.19
//*
//***********************************
//

#include <string>
#include <time.h>
#include <TString.h>

#include "Include/DrawMyfunc.h"
#include "Include/LZWfunc.h"

const char path[2014]="/mnt/e/sub-system/Analysis/Rootfile/calibration/Separate/KT890_A5";

using namespace std;
const int Nch = 1;
void drawline(float x1, float x2){

	TLine* line1 = new TLine(x1,gPad->VtoPixel(gPad->GetUymin()),x1,gPad->VtoPixel(gPad->GetUymax()));
    TLine* line2 = new TLine(x2,gPad->VtoPixel(gPad->GetUymin()),x2,gPad->VtoPixel(gPad->GetUymax()));
    
    line1->SetLineStyle(2);
    line1->SetLineColor(6);
    line1->Draw("same");
    line2->SetLineStyle(2);
    line2->SetLineColor(6);
    line2->Draw("same");
}

void MCP_SPEchar(int mych=0,const char* rootname="KT881_2100_A5A6A7"){

	gStyle->SetOptFit(1112);


	
//** Set your parameters**
//*	
	char name[1024];
    char buff[1024];
    
	string Chname[Nch];
	string str(rootname);
	for (int i = 0; i < Nch; i++ ){

    Chname[Nch-1-i]=str.substr(str.length()-2*(i+1),2);
cout<< " the ch name is: "<<Chname[i]<<endl;
	} 
//return;

    //**
	//**

	//double M=10;
	double M=1;
	int qN=0;
	int rbu=4;
	double fac2gaus=4;
	int rbq=4;
	double spefac1 = 15;
	double spefac2 = 12;
	double ampfac1 = 5;
	double ampfac2 = 15;
	RANGE def={-9999,9999};
	RANGE x={0,50};
	//RANGE x={0,100};
	//RANGE y={-5e-3,150e-3};
	RANGE y={-5e-3,150e-3};
    //RANGE q={-0.1,3};
    RANGE q={-0.05,2};
	RANGE r={0,1};
    //RANGE t={15,17};
    RANGE t={16,18};
	
	//**
	//** Hist range !! **
    charRANGE chR={
		x,
        y,
        q,
        r,
        def,
        def,
        t
	};//x,y,q,r,bl,blrms,t
	

	//**
	//** CUT RANGE !!**
    charRANGE chcut={
		def,
        def,
        {0.1,9999},
        //{0.4,9999},
        {0.3,9999},
        def,
        def,
        def
	};//blrms,y 
	chcut.x.L = 36;
	chcut.x.R = 37.5;

	int rbt=128,rbU=16;
	int rb=2;
	double res=2;
    double fac=1.8;
	double fac2=0.4;
    int iter=1;

	OPTION opt={rb,rbU,iter,fac2};
//*
//** end **
//

    EVENT A[3];
	EVENT B;
    double t0;
	double p[20];
	//sprintf(buff,"/mnt/f/experiment/FTOF/SIPM/labtest/0108/BFP650_v2_s2/123.dat");
	ofstream op;
	//cout<<"Build your data file "<<endl;
    DrawMyfunc draw;
    LZWfunc lzw;

	TCanvas *c[4];
	for (int i=0;i<4;i++){
    	sprintf(buff,"C%d",i);
        c[i]= new TCanvas(buff,buff,800,600);
    	draw.SetPad(gPad,0.12,0.12,0.1,0.1);
	}
//return;

    TH1F *hq = new TH1F("hq",";charge (pC);Counts",1000,chR.q.L,chR.q.R);
	draw.Hist(hq,"charge (pC)","Counts",2);

	TH1F *ha = new TH1F("ha",";Amp(V);Counts",500,chR.y.L,chR.y.R);
	draw.Hist(ha,"Amp(V)","Counts",2);
	
	TH1F *hr = new TH1F("hr",";risetime (ns);Counts",200,chR.r.L,chR.r.R);
	draw.Hist(hr,"Risetime (ns)","Counts",2);
	
	TH1F *ht = new TH1F("ht","",400,chR.t.L,chR.t.R);
	draw.Hist(ht,"STR (ps)","Counts",2);

	TH1F *hbl = new TH1F("hbl","",400,chR.bl.L,chR.bl.R);
	draw.Hist(hbl,"baseline (V)","Counts",2);
	TGaxis::SetMaxDigits(3);
	
	TH1F *hblrms = new TH1F("hblrms","",800,chR.blrms.L,chR.blrms.R);
	draw.Hist(hblrms,"baseline RMS (V)","Counts",2);

	TH2F* hqy=new TH2F("hqy","",200,chR.y.L,chR.y.R,200,chR.q.L,chR.q.R);
    draw.Hist(hqy,"Amp (V)","Charge (pC)");

	TH1F *hx_origin = new TH1F("hx_origin",";Time (ns);Counts",2e3,chR.x.L,chR.x.R);
	draw.Hist(hx_origin,"Time (ns)","Counts",2);
	
	sprintf(buff,"%s/%s.root",path,rootname);

	cout<<"=====>> Start to read rootfile: "<<buff<<endl;
	
    TFile *f1 = new TFile(buff,"READ");
    TTree *t1 = (TTree*)f1->Get("Pico");
	//string path = f1->GetPath();
	//path=path.substr(0,path.length()-7);
	
	
	
	//double Q[4];
	//t1->SetBranchAddress("MCP2_CFDtime",&t0);
	
	for(int i = 0; i<Nch; i++){
	sprintf(buff,"MCP%d_CFDtime",i+1);
    t1->SetBranchAddress(buff,A[i].CFD);
	sprintf(buff,"MCP%d_all_charge",i+1);
	t1->SetBranchAddress(buff,A[i].charge);
	sprintf(buff,"MCP%d_rise_time",i+1);
	t1->SetBranchAddress(buff,&A[i].rise);
	sprintf(buff,"MCP%d_global_maximum_y",i+1);
	t1->SetBranchAddress(buff,&A[i].Amp);
	sprintf(buff,"MCP%d_global_maximum_x",i+1);
	t1->SetBranchAddress(buff,&A[i].x);
	sprintf(buff,"MCP%d_baseline_rms",i+1);
	t1->SetBranchAddress(buff,&A[i].blrms);
	sprintf(buff,"MCP%d_baseline_level",i+1);
	t1->SetBranchAddress(buff,&A[i].bl);
	}
	
	t1->SetBranchAddress("MCP2_CFDfrac",p);
	t1->SetBranchAddress("MCP2_CFDtime",B.CFD);
	t1->SetBranchAddress("MCP2_all_charge",B.charge);
	t1->SetBranchAddress("MCP2_rise_time",&B.rise);
	t1->SetBranchAddress("MCP2_global_maximum_y",&B.Amp);
	t1->SetBranchAddress("MCP2_global_maximum_x",&B.x);
	t1->SetBranchAddress("MCP2_baseline_rms",&B.blrms);
	t1->SetBranchAddress("MCP2_baseline_level",&B.bl);

	//TCut c_x = "MCP2_global_maximum_x>410&&MCP2_global_maximum_x<412.5";
	//TCut c_x = "MCP2_global_maximum_x>530&&MCP2_global_maximum_x<532";
	//TCut c_blrms = "MCP2_baseline_rms<0.1e-3";
	TH1F* hx=(TH1F*)hx_origin->Clone("hx");
	hx->Draw();
	//return;
	//for(int j = mych ; j<mych+1; j++ ){
RANGE xrange;
int other1;
int other2;
	
	for(int j = mych ; j<mych+1; j++ ){
	other1 = j+1>=Nch? j+1-Nch:j+1;
	//other2 = j+2>=Nch? j+2-Nch:j+2;
	
	chcut.x = x;
	xrange = x;
	
	hq->Reset();
	ha->Reset();
	hr->Reset();
	hbl->Reset();
	hblrms->Reset();
	hqy->Reset();
	ht->Reset();
	
	

	cout<< "=====>> Analysis channel: "<<Chname[j]<<endl;

	sprintf(buff,"%s/%s_%s.dat",path,rootname,Chname[j].c_str());
	//op.open(buff,ios::app);
	op.open(buff,ios::trunc);
	cout<<"=====>> Start to open datafile: "<<buff<<endl;
	
	c[0]->cd();
	hx->Reset();
	
	sprintf(buff,"MCP%d_global_maximum_x>>hx",j+1);
	t1->Draw(buff);
    //find the cut range;
	
	//hx->GetXaxis()->SetRangeUser(2,30);
    xrange.L = hx->GetBinCenter(hx->GetMaximumBin())-1.5;
    xrange.R = hx->GetBinCenter(hx->GetMaximumBin())+1.5;
	cout<<xrange.L<<"\t"<<xrange.R<<endl;
//return;
    TF1* gx=lzw.gausfit(hx,1,1.5,xrange);
//return;
	
    xrange.L = gx->GetParameter(1)-1.5*gx->GetParameter(2);
	xrange.R = gx->GetParameter(1)+2.*gx->GetParameter(2);
	cout<<"global maximum x cut range = "<<xrange.L<<"\t"<<xrange.R<<endl;
	
	op<<chcut.x.L<<"\t"<<chcut.x.R<<endl;
	drawline(chcut.x.L,chcut.x.R);
	
	sprintf(buff,"%s/%s_%s_Xcut.png",path,rootname,Chname[j].c_str());
    c[0]->SaveAs(buff);
	//chcut.x.L = xrange.L;
	//chcut.x.R = xrange.R;
	//return;


	//gPad->Modified(); 
	//gPad->Update(); 
 
 //return;
 
 int N=t1->GetEntries();
	for(int i=0;i<N;i++){
		t1->GetEntry(i);
		//if(Q2>qlimit1&&Q2<qlimit2)


        //pick up the true signal from noise
			if (A[j].x>chcut.x.L&&
            A[j].x<chcut.x.R&&
            A[j].blrms>chcut.blrms.L&&
            A[j].blrms<chcut.blrms.R&&
            A[j].bl<chcut.bl.R&&
            A[j].bl>chcut.bl.L&&
			A[j].rise>chcut.r.L&&
			A[j].rise<chcut.r.R
			//&&A[other1].charge[0]<0.08
			//&&A[other2].charge[0]<0.05
			//A[other1].rise>chcut.r.L&&
			//A[other2].rise>chcut.r.L
			

			//&&A[j].Q>chcut.q.L
			)
			
			{
				hq->Fill(A[j].charge[qN]);
				ha->Fill(A[j].Amp);
				//cout<<Q[1]<<endl;
				hbl->Fill(A[j].bl);
				hblrms->Fill(A[j].blrms);
				hqy->Fill(A[j].Amp,A[j].charge[qN]);
				if(A[j].charge[qN]>chcut.q.L) {
					hr->Fill(A[j].rise);
					ht->Fill(A[j].CFD[3]-B.CFD[3]);
				}
			}
	}


	//t1->Draw("MCP2_global_maximum_y>>ha",c_x&&c_blrms);	
	//c[1]->cd();
	cout<<"program running normly"<<endl;	

//	Charge spectrum
	c[2]->cd();
	c[2]->SetLogy();
	//sprintf(buff,"%s_charge",name);
    //lzw.Set_name(buff);
	//hq->Draw();
	//return;
	//TF1* fq=lzw.SPSfit(hq,4,chR.q,facspe);
	TF1* fq=lzw.mcpSPfit(hq,rbq,chR.q,spefac1,spefac2);
	double Gain=(fq->GetParameter(4)-fq->GetParameter(1))*1e-12/1.6e-19;
	//double Gain=fq->GetParameter(4)*1e-12/1.6e-19;
	cout<<"Gain="<<Gain<<endl;
	sprintf(buff,"Gain=%0.2e",Gain);
	TLatex *l=draw.Latex(0.4,0.2,buff);
	l->Draw();
	sprintf(buff,"%s/%s_%s_charge.png",path,rootname,Chname[j].c_str());
    c[2]->SaveAs(buff);
//return;

//	amplitude spectrum
	sprintf(buff,"%s_amp",name);
    lzw.Set_name(buff);
	c[1]->cd();
	ha->Draw();
	c[1]->SetLogy();
	TF1* fa=lzw.mcpSPfit(ha,2,chR.y,ampfac1,ampfac2);
	double SPEAmp=(fa->GetParameter(4)-fa->GetParameter(1))*1e3;
	cout<<"SPEAmp="<<SPEAmp<<endl;
	sprintf(buff,"SPE Amp=%0.2f mV",SPEAmp);
	TLatex *l1=draw.Latex(0.4,0.2,buff);
	l1->Draw();
	sprintf(buff,"%s/%s_%s_amp.png",path,rootname,Chname[j].c_str());
    c[1]->SaveAs(buff);
//	risetime spectrum
	c[3]->cd();
	TF1* r1=lzw.gausfit(hr,2,2,chcut.r);
	//return;
	double risetime=r1->GetParameter(1);
	sprintf(buff,"%s/%s_%s_risetime.png",path,rootname,Chname[j].c_str());
    c[3]->SaveAs(buff);

//**	time resolution spectrum
	c[0]->cd();
	ht->Draw();
	//return;
	//TF1* gt=lzw.gausfit(ht,2,2,0.7,&chR.t);
	TF1* gt=lzw.twogausfit(ht,fac2gaus,2,chR.t);
	sprintf(buff,"%s/%s_%s_TR.png",path,rootname,Chname[j].c_str());
	c[0]->SaveAs(buff);
	double TT=gt->GetParameter(1);
	double TTS=gt->GetParameter(2);
	op<<SPEAmp<<endl;
	op<<Gain<<endl;
	op<<risetime<<endl;
	op<<TT<<endl;
	op<<TTS<<endl;
	op.close();
	
	}

	f1->Close();
	

}


//int main()
void drawresults(int mych, int tubeid = 881, const int pointNum = 4, double HVbase=2200  )
{
	char name[1024];
    char buff[1024];
	
	ifstream input;
	double SPEAmp[pointNum];
	double Gain[pointNum];
	double risetime[pointNum];
	double TT[pointNum];
	double TTS[pointNum];
	double XL[pointNum];
	double XR[pointNum]; 
	double HV[pointNum];

	;
	int HVstep = 50;
	for (int i = 0; i < pointNum; i++ ){
	HV[i] = HVbase +i*HVstep;
	sprintf(buff,"%s/KT%d%g_KT%d_A%d_A%d.dat",path,tubeid,HV[i],tubeid,mych,mych);
	cout<< "open datafile: "<<buff<<endl; 

	//op.open(buff,ios::app);
	input.open(buff);
	if(!input){
		cout<<"file can't be found: \n"<<buff<<endl;
		return 0;
    }
    while(input>>XL[i]>>XR[i]){
		cout<<"enter the file: "<<buff<<endl;
		

		input>>SPEAmp[i];
		input>>Gain[i];
		input>>risetime[i];
		input>>TT[i];
		input>>TTS[i];
		risetime[i]=risetime[i]*1e3;
		TTS[i]=abs(TTS[i])*1e3;
		 cout<<Gain[i]<<"\t"<<risetime[i]<<"\t"<<TT[i]<<"\t"<<TTS[i]<<endl;
	}
	input.close();
	}

sprintf(name,"%s/KT%d_A%d",path,tubeid,mych);

TGraph* g0 = new TGraph(pointNum,HV,SPEAmp); 
TGraph* g1 = new TGraph(pointNum,HV,Gain); 
TGraph* g2 = new TGraph(pointNum,HV,risetime); 
TGraph* g3 = new TGraph(pointNum,HV,TT); 
TGraph* g4 = new TGraph(pointNum,HV,TTS); 

/*
TH1F *h1 = new TH1F("h1","",100,-1,1);
h1->FillRandom("gaus",1e3);
TF1 *f1 = new TF1("f1","sin(x)",-1,1);
*/
TCanvas *c0 = new TCanvas("c0","",800,600);
TCanvas *c1 = new TCanvas("c1","",800,600);
TCanvas *c2 = new TCanvas("c2","",800,600);
TCanvas *c3 = new TCanvas("c3","",800,600);
TCanvas *c4 = new TCanvas("c4","",800,600);
c0->SetGrid();
c1->SetGrid();
c2->SetGrid();
c3->SetGrid();
c4->SetGrid();

c0->cd();
g0->Draw("AP");
DrawMyfunc mydraw;
mydraw.SetPad(c0,0.12,0.14,0.1,0.05);
//mydraw.Graph(g1,"NPE","TimeRes (ps)",1.5,20,4);

mydraw.Graph(g0,"WorkVol (V)","SPEAmp (mV)",1.5,20,4);
g0->GetYaxis()->SetRangeUser(0,120);
//g0->GetYaxis()->SetRangeUser(0,SPEAmp[pointNum-1]*1.2);

sprintf(buff,"%sSPEAmp.png",name);
c0->SaveAs(buff);

c1->cd();
g1->Draw("AP");
mydraw.SetPad(c1,0.12,0.14,0.1,0.05);
//mydraw.Graph(g1,"NPE","TimeRes (ps)",1.5,20,4);

mydraw.Graph(g1,"WorkVol (V)","Gain",1.5,20,4);
g1->GetYaxis()->SetRangeUser(1e5,20.0e6);

sprintf(buff,"%sgain.png",name);
c1->SaveAs(buff);

c2->cd();
g2->Draw("AP");
mydraw.SetPad(c2,0.12,0.14,0.1,0.05);
//mydraw.Graph(g1,"NPE","TimeRes (ps)",1.5,20,4);

mydraw.Graph(g2,"WorkVol (V)","Risetime  (ps)",1.5,20,4);
g2->GetYaxis()->SetRangeUser(0.05*1e3,0.6*1e3);

sprintf(buff,"%srisetime.png",name);
c2->SaveAs(buff);
/*
c3->cd();
g3->Draw("AP");
mydraw.SetPad(c3,0.12,0.14,0.1,0.05);
//mydraw.Graph(g1,"NPE","TimeRes (ps)",1.5,20,4);

mydraw.Graph(g3,"WorkVol (V)","Transit time (ns)",1.5,20,4);
g3->GetYaxis()->SetRangeUser(25,37);

sprintf(buff,"%s_transittime.png",name);
c3->SaveAs(buff);
 */

c4->cd();
g4->Draw("AP");
mydraw.SetPad(c4,0.12,0.14,0.1,0.05);
//mydraw.Graph(g1,"NPE","TimeRes (ps)",1.5,20,4);

mydraw.Graph(g4,"WorkVol (V)","TTS (ps)",1.5,20,4);
g4->GetYaxis()->SetRangeUser(0.01*1e3,0.100*1e3);

sprintf(buff,"%s_tts.png",name);
c4->SaveAs(buff);


}


void MCP_SPEchar(int colNum, int rowNum){

	gStyle->SetOptFit(1112);


	
//** Set your parameters**
//*	
	char name[1024];
    char buff[1024];
    
    
	//**
	//********* G3 colum 1,1-7*******
	//**
	//const char* rootname="/mnt/f/MCP/7-26/KT0900_G3_Scanspot";

    //**
	//******** G3 colum 2,1-7****
	//**
	//const char* rootname="/mnt/f/MCP/7-27/KT0900_G3_scan/KT0900_G3_scanspot";
	
	//**
	//******* G4 colum 1 & colum 2,0-5 ***
	//**
	//const char* rootname="/mnt/f/MCP/7-27/KT0900_scan/KT0900_scanG4_";

  	//**
	//******** G3 verify **********
	//
	//const char* rootname="/mnt/f/MCP/7-27/KT0900_G3_scan_verify/KT0900_G3_scan_verifyspot";

	//**
	//** colum 3 , 1-7
	//**
	//const char* rootname="/mnt/f/MCP/7-28/col/colspot";

	//**
	//** anode scan
	//**
	//const char* rootname="/mnt/f/MCP/7-28/scan/scanA";

	//**
	//** invert the tube, and scan column3
	//**
	//const char* rootname="/mnt/f/MCP/7-29/invert_spot";

	//**
	//** invert the tube, and scan column1
	//**
	//const char* rootname="/mnt/f/MCP/7-29/colum/columinvert_spot1_";

	//**
	//** invert the tube, and scan column1
	//**
	//const char* rootname="/mnt/f/MCP/7-30/base2/base2A";

	//**
	//** change the base to version 5_a, and scan column1, 0-3 of group 1
	//**
	//const char* rootname="/mnt/f/MCP/7-30/base5_a/base5_aA";


    //**
	//** change the base to version 5_a, and scan column1, 0-3 of group 1
	//**
	const char* rootname="/mnt/f/MCP/7-30/KT0881_base5_b/KT0881_base5_bG4_";

	//double M=10;
	double M=1;
	int qN=0;
	int rbu=4;
	double facspe=6.5;
	int rbq=4;
	double spefac1 =12;
	double spefac2 = 20;
	RANGE def={-9999,9999};
	RANGE x={0,100};
	//RANGE y={-5e-3,100e-3};
	RANGE y={-5e-3,50e-3};
    //RANGE q={-0.1,1};
    RANGE q={-0.01,1.0};
	RANGE r={0,1};
    RANGE t={15,17};
	
	//**
	//** Hist range !! **
    charRANGE chR={
		x,
        y,
        q,
        r,
        def,
        def,
        t
	};//x,y,q,r,bl,blrms,t
	

	//**
	//** CUT RANGE !!**
    charRANGE chcut={
		def,
        def,
        def,
        //{0.4,9999},
        {0.3,9999},
        def,
        def,
        def
	};//blrms,y 

	int rbt=128,rbU=16;
	int rb=2;
	double res=2;
    double fac=1.8;
	double fac2=0.4;
    int iter=1;

	OPTION opt={rb,rbU,iter,fac2};
//*
//** end **
//

    EVENT A,B;
    double t0;
	double p[20];
	//sprintf(buff,"/mnt/f/experiment/FTOF/SIPM/labtest/0108/BFP650_v2_s2/123.dat");
	ofstream op;
	//cout<<"Build your data file "<<endl;
    DrawMyfunc draw;
    LZWfunc lzw;

	TCanvas *c[4];
	for (int i=0;i<4;i++){
    	sprintf(buff,"C%d",i);
        c[i]= new TCanvas(buff,buff,800,600);
    	draw.SetPad(gPad,0.12,0.12,0.1,0.1);
	}
//return;

    TH1F *hq = new TH1F("hq",";charge (pC);Counts",1000,chR.q.L,chR.q.R);
	draw.Hist(hq,"charge (pC)","Counts",2);

	TH1F *ha = new TH1F("ha",";Amp(V);Counts",500,chR.y.L,chR.y.R);
	draw.Hist(ha,"Amp(V)","Counts",2);
	
	TH1F *hr = new TH1F("hr",";risetime (ns);Counts",200,chR.r.L,chR.r.R);
	draw.Hist(hr,"Risetime (ns)","Counts",2);
	
	TH1F *ht = new TH1F("ht","",400,chR.t.L,chR.t.R);
	draw.Hist(ht,"STR (ps)","Counts",2);

	TH1F *hbl = new TH1F("hbl","",400,chR.bl.L,chR.bl.R);
	draw.Hist(hbl,"baseline (V)","Counts",2);
	TGaxis::SetMaxDigits(3);
	
	TH1F *hblrms = new TH1F("hblrms","",800,chR.blrms.L,chR.blrms.R);
	draw.Hist(hblrms,"baseline RMS (V)","Counts",2);

	TH2F* hqy=new TH2F("hqy","",200,chR.y.L,chR.y.R,200,chR.q.L,chR.q.R);
    draw.Hist(hqy,"Amp (V)","Charge (pC)");

	TH1F *hx_origin = new TH1F("hx_origin",";Time (ns);Counts",2e3,chR.x.L,chR.x.R);
	draw.Hist(hx_origin,"Time (ns)","Counts",2);
	
	//for (int k=0; k<7; k++){
	for (int k=rowNum; k<rowNum+1; k++){
	sprintf(name,"%s%d_%d",rootname,colNum,k);
	//sprintf(name,"%s%d",rootname,k);

	
	sprintf(buff,"%s.root",name);
    TFile *f1 = new TFile(buff,"READ");
    TTree *t1 = (TTree*)f1->Get("Pico");
	string path = f1->GetPath();
	path=path.substr(0,path.length()-7);
	sprintf(buff,"%s.dat",path.c_str());
	//op.open(buff,ios::app);
	op.open(buff,ios::trunc);
	cout<<name<<endl;
	
	
	
	//double Q[4];
	//t1->SetBranchAddress("MCP2_CFDtime",&t0);
	
    t1->SetBranchAddress("MCP1_CFDtime",A.CFD);
	t1->SetBranchAddress("MCP1_CFDfrac",p);
	t1->SetBranchAddress("MCP1_all_charge",A.charge);
	t1->SetBranchAddress("MCP1_rise_time",&A.rise);
	t1->SetBranchAddress("MCP1_global_maximum_y",&A.Amp);
	t1->SetBranchAddress("MCP1_global_maximum_x",&A.x);
	t1->SetBranchAddress("MCP1_baseline_rms",&A.blrms);
	t1->SetBranchAddress("MCP1_baseline_level",&A.bl);
	t1->SetBranchAddress("MCP2_CFDtime",B.CFD);
	t1->SetBranchAddress("MCP2_all_charge",B.charge);
	t1->SetBranchAddress("MCP2_rise_time",&B.rise);
	t1->SetBranchAddress("MCP2_global_maximum_y",&B.Amp);
	t1->SetBranchAddress("MCP2_global_maximum_x",&B.x);
	t1->SetBranchAddress("MCP2_baseline_rms",&B.blrms);
	t1->SetBranchAddress("MCP2_baseline_level",&B.bl);

	//TCut c_x = "MCP2_global_maximum_x>410&&MCP2_global_maximum_x<412.5";
	//TCut c_x = "MCP2_global_maximum_x>530&&MCP2_global_maximum_x<532";
	//TCut c_blrms = "MCP2_baseline_rms<0.1e-3";

	c[0]->cd();
	TH1F* hx=(TH1F*)hx_origin->Clone("hx");
	t1->Draw("MCP1_global_maximum_x>>hx");

    //find the cut range;
	
	hx->GetXaxis()->SetRangeUser(28,30);
    chR.x.L=hx->GetBinCenter(hx->GetMaximumBin())-0.5;
    chR.x.R = hx->GetBinCenter(hx->GetMaximumBin())+0.5;
	cout<<chR.x.L<<"\t"<<chR.x.R<<endl;
//return;
    TF1* gx=lzw.gausfit(hx,1,1.5,&chR.x);
//return;

    chcut.x.L = gx->GetParameter(1)-3.*gx->GetParameter(2);
	chcut.x.R = gx->GetParameter(1)+3.*gx->GetParameter(2);
	cout<<"global maximum x cut range = "<<chcut.x.L<<"\t"<<chcut.x.R<<endl;



	//gPad->Modified(); 
	//gPad->Update(); 
 
 //return;
 
 int N=t1->GetEntries();
	for(int i=0;i<N;i++){
		t1->GetEntry(i);
		//if(Q2>qlimit1&&Q2<qlimit2)


        //pick up the true signal from noise
			if (A.x>chcut.x.L&&
            A.x<chcut.x.R&&
            A.blrms>chcut.blrms.L&&
            A.blrms<chcut.blrms.R&&
            A.bl<chcut.bl.R&&
            A.bl>chcut.bl.L
			//&&A.Q>chcut.q.L
			)
			
			{
				hq->Fill(A.charge[qN]);
				ha->Fill(A.Amp);
				//cout<<Q[1]<<endl;
				hr->Fill(A.rise);
				hbl->Fill(A.bl);
				hblrms->Fill(A.blrms);
				hqy->Fill(A.Amp,A.charge[qN]);
				if(A.charge[0]>0.03) ht->Fill(A.CFD[3]-B.CFD[3]);
				
			}
	}

/*
	double xl=ha->GetBinCenter(ha->FindFirstBinAbove(2));
	double xr=ha->GetBinCenter(ha->FindLastBinAbove(5));
	sigmaA=(xr-xl)/5/5;
	cout<<"sigma Amp="<<sigmaA<<endl;

	xl=hq->GetBinCenter(hq->FindFirstBinAbove(2));
	xr=hq->GetBinCenter(hq->FindLastBinAbove(5));
	sigmaq=(xr-xl)/5/5;
	cout<<"sigma Charge="<<sigmaq<<endl;

    lzw.Set_name(name);
    lzw.Set_charcut(chcut);
    RANGE t=chR.t;
    RANGE U={0,0.05};
*/

	//t1->Draw("MCP2_global_maximum_y>>ha",c_x&&c_blrms);	
	//c[1]->cd();
	cout<<"program running normly"<<endl;	
	/*
	TSpectrum *s = new TSpectrum(20,0.8);
	int nfound = s->Search(ha,2e-3,"",0.05);
	printf("Found %d candidate peaks to fit\n",nfound);
	ha->Draw();
	//return;
	double *xpeaks;
	xpeaks = s->GetPositionX();
	*/
//	Charge spectrum
	c[2]->cd();
	c[2]->SetLogy();
	sprintf(buff,"%s_charge",name);
    lzw.Set_name(buff);
	//hq->Draw();
	//return;
	//TF1* fq=lzw.SPSfit(hq,4,chR.q,facspe);
	TF1* fq=lzw.mcpSPfit(hq,rbq,chR.q,spefac1,spefac2);
	double Gain=(fq->GetParameter(4)-fq->GetParameter(1))*1e-12/1.6e-19;
	//double Gain=fq->GetParameter(4)*1e-12/1.6e-19;
	cout<<"Gain="<<Gain<<endl;
	sprintf(buff,"Gain=%0.2e",Gain);
	TLatex *l=draw.Latex(0.4,0.2,buff);
	l->Draw();
	sprintf(buff,"%s_charge.png",name);
    c[2]->SaveAs(buff);
//return;

//	amplitude spectrum
	sprintf(buff,"%s_amp",name);
    lzw.Set_name(buff);
	c[1]->cd();
	ha->Draw();
	c[1]->SetLogy();
	TF1* fa=lzw.mcpSPfit(ha,1,chR.y,8,5);
	double SPEAmp=(fa->GetParameter(4)-fa->GetParameter(1))*1e3;
	cout<<"SPEAmp="<<SPEAmp<<endl;
	sprintf(buff,"SPE Amp=%0.2f mV",SPEAmp);
	TLatex *l1=draw.Latex(0.4,0.2,buff);
	l1->Draw();
	sprintf(buff,"%s_amp.png",name);
    c[1]->SaveAs(buff);



/*	TF1* fa=lzw.fpeaksfit(ha,maxpk,res,sigmaA,thrd);
	gausPAR* gPa = lzw.Get_PeakPar();
	double speA=(gPa+1)->m-gPa->m;
	chcut.y.L=gPa->m-2*gPa->s;
	chcut.y.R=gPa->m+2*gPa->s;
	cout<<"chcut y left ="<<chcut.y.L<<", chcut y right ="<<chcut.y.R<<endl;
*/

	//gausPAR* gPq = lzw.Get_PeakPar();
	//double speq=(gPq+1)->m-gPq->m;
	//double G=speq*1e-12/1.6e-19;
	


//	risetime spectrum
	c[3]->cd();
	TF1* r1=lzw.twogausfit(hr,fac2,rb,chR.r);
	//return;
	double risetime=r1->GetParameter(1);
	sprintf(buff,"%s_risetime.png",name);
    c[3]->SaveAs(buff);

//**	time resolution spectrum
c[0]->cd();
ht->Draw();
TF1* gt=lzw.gausfit(ht,2,2,0.7,&chR.t);
sprintf(buff,"%s_TR.png",name);
c[0]->SaveAs(buff);
double TT=gt->GetParameter(1);
double TTS=gt->GetParameter(2);
op<<Gain<<endl;
op<<risetime<<endl;
op<<TT<<endl;
op<<TTS<<endl;
op.close();
}
/* chcut.q.L=0.03;
int N=t1->GetEntries();
	for(int i=0;i<N;i++){
		t1->GetEntry(i);
		//if(Q2>qlimit1&&Q2<qlimit2)


        //pick up the true signal from noise
			if (A.x>chcut.x.L&&
            A.x<chcut.x.R&&
            A.blrms>chcut.blrms.L&&
            A.blrms<chcut.blrms.R&&
            A.bl<chcut.bl.R&&
            A.bl>chcut.bl.L
			&&A.Q>chcut.q.L
			)
			{
				ht->Fill(A.CFD[3]-B.CFD[3]);
				
			}
	}
	//lzw.Set_charcut(chcut);
    //lzw.CH2Correction(t1,ch,p,chR,chcut,opt,name);
	//double STR=tt->GetParameter(2);
	//double STRerr=tt->GetParError(2);
	


//	pedestal spectrum
	c[5]->cd();
	TF1* b1=lzw.twoguasfit(hbl,fac2,rb,chR.bl);
	double pedmean=b1->GetParameter(1)*1e3;
	double pedmeansigma=b1->GetParameter(2)*1e3;
	sprintf(buff,"%s_baseline.png",name);
    c[5]->SaveAs(buff);

//	pedestal rms spectrum
	c[6]->cd();
	TF1* b2=lzw.twoguasfit(hblrms,fac2,rb*2,chR.blrms);
	double pedRMS=b2->GetParameter(1)*1e3;
	sprintf(buff,"%s_baselinerms.png",name);
    c[6]->SaveAs(buff);

	c[0]->cd();
	hqy->Draw("colz");
	sprintf(buff,"%s_Qvsy.png",name);
	c[0]->SaveAs(buff);
*/

/*	
	
	cout<<"The spe amplitude (mV) = "<<speA*1e3<<endl;
	cout<<"The spe charge (pC) = "<<speq<<endl;
	cout<<"The Gain = "<<G<<endl;
	cout<<"The STR (ps) = "<<STR*1e3<<"\t"<<STRerr*1e3<<endl;
	cout<<"The ped amp mean (mV) = "<<pedmean<<endl;
	cout<<"The ped amp sigma (mV) = "<<pedmeansigma<<endl;
	cout<<"The ped amp RMS (mV) = "<<pedRMS<<endl;
	cout<<"The risetime (ns)= "<<risetime<<endl;
	
	string time = lzw.getTime();
	op<<"**************"<<endl;
	op<<"*"<<endl;
	op<<"* $ SiPM SPE Fit Function $"<<endl;
	op<<"*"<<endl;
    op<<"* Date: "<< time << endl;
	op<<"*"<<endl;
	op<<"* ROOT file:  "<<name<<endl;
	op<<"*"<<endl;
	op<<"* The cut set: "<<endl;
	op<<"*"<<endl;
	op<<"** The global x cut range: "<<chcut.x.L<<"\t"<<chcut.x.R<<endl;
	op<<"** The blrms cut range: "<<chcut.blrms.L<<"\t"<<chcut.blrms.R<<endl;
	op<<"** The amp cut range: "<<chcut.y.L<<"\t"<<chcut.y.R<<endl;
	op<<"*"<<endl;
	op<<"*"<<endl;
	op<<"* The results:  "<<endl;
	op<<"*"<<endl;
	op<<"** The spe amplitude (mV) = "<<speA*1e3<<endl;
	op<<"** The spe charge (pC) = "<<speq<<endl;
	op<<"** The Gain = "<<G<<endl;
	op<<"** The STR (ps) = "<<STR*1e3<<"\t"<<STRerr*1e3<<endl;
	op<<"** The ped amp mean (mV) = "<<pedmean<<endl;
	op<<"** The ped amp sigma (mV) = "<<pedmeansigma<<endl;
	op<<"** The ped amp RMS (mV) = "<<pedRMS<<endl;
	op<<"** The risetime (ns)= "<<risetime<<endl;
	op<<"*"<<endl;
	op<<"********* end *********"<<endl;
	op<<"\n\n\n"<<endl;
	//double pedqsigma=q1->GetParameter(2);
	//TF1* q2=gausfit(hq,1,qlimit1,qlimit2);
*/

/*
	TF1* a1=gausfit(ha,1,aRL,alimit1);
	TF1* a2=gausfit(ha,1,alimit1,alimit2);
	double pe1A=a1->GetParameter(1);
	double pe2A=a2->GetParameter(1);
	double speA=pe2A-pe1A;
	a1->Draw("same");
	a2->Draw("same");
    sprintf(buff,"%s_amp.png",name);
    c[1]->SaveAs(buff);

	c[2]->cd();
    //t1->Draw("MCP2_all_charge>>hq",c_x&&c_blrms);
    TF1* q1=gausfit(hq,1,qRL,qlimit1);
	//return;
	TF1* q2=gausfit(hq,1,qlimit1,qlimit2);
	
	double pe1q=q1->GetParameter(1);
	double pe2q=q2->GetParameter(1);
	double speq=pe2q-pe1q;
	//double speq=pe1q-pedq;
	double G=speq*1e-12/1.6e-19;
	q1->Draw("same");
	q2->Draw("same");
	sprintf(buff,"%s_charge.png",name);
    
    c[2]->SaveAs(buff);
	
	c[3]->cd();
	//t1->Draw("MCP2_rise_time>>hr",c_x&&c_blrms);
    //TF1* r1=gausfit(hr,1,0,5);
	TF1* r1=twoguasfit(hr,&rRL,&rRR,0.3,1);
	//return;
	double risetime=r1->GetParameter(1);
	//double pedqsigma=q1->GetParameter(2);
	//TF1* q2=gausfit(hq,1,qlimit1,qlimit2);
	
	//double pe1q=q2->GetParameter(1);
	//double speq=pe1q-pedq;
	//double G=speq*1e-12/1.6e-19;
	//q1->Draw("same");
	//q2->Draw("same");
	sprintf(buff,"%s_risetime.png",name);
    
    c[3]->SaveAs(buff);
	//return;
	
	c[4]->cd();		
	
	TF1* t=twoguasfit(ht,&tRL,&tRR,0.2,2);
	double STR=t->GetParameter(2);
	double STRerr=t->GetParError(2);
	sprintf(buff,"%s_STR.png",name);
    
    c[4]->SaveAs(buff);

	c[5]->cd();
	TF1* b1=twoguasfit(hbl,&bRL,&bRR,0.4,1);
	double pedmean=b1->GetParameter(1)*1e3;
	double pedmeansigma=b1->GetParameter(2)*1e3;
	sprintf(buff,"%s_baseline.png",name);
    
    c[5]->SaveAs(buff);

	c[6]->cd();
	TF1* b2=twoguasfit(hblrms,&brmsRL,&brmsRR,0.4,1);
	double pedRMS=b2->GetParameter(1)*1e3;
	sprintf(buff,"%s_baselinerms.png",name);
    
    c[6]->SaveAs(buff);

	c[0]->cd();
	hqy->Draw("colz");
	sprintf(buff,"%s_Qvsy",name);
	c[0]->SaveAs(buff);
	
	
	cout<<"The spe amplitude (mV) = "<<pe1A*1e3<<endl;
	cout<<"The spe charge (pC) = "<<speq<<endl;
	cout<<"The Gain = "<<G<<endl;
	cout<<"The STR (ps) = "<<STR*1e3<<"\t"<<STRerr*1e3<<endl;
	cout<<"The ped amp mean (mV) = "<<pedmean<<endl;
	cout<<"The ped amp sigma (mV) = "<<pedmeansigma<<endl;
	cout<<"The ped amp RMS (mV) = "<<pedRMS<<endl;
	cout<<"The risetime (ns)= "<<risetime<<endl;


	
	*/
}
