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

using namespace std;
void MCP_SPEchar(){

	gStyle->SetOptFit(1112);


	
//** Set your parameters**
//*	
	char name[1024];
    char buff[1024];
    //const char* rootname="/mnt/f/experiment/FTOF/MCP/4-12/systemR3809-SPE";
    const char* rootname="/mnt/f/experiment/FTOF/MCP/4-10/test/HV-SPE2200";
	sprintf(name,"%s",rootname);
  
	//double M=10;
	double M=1;
	int qN=0;
	int rbu=4;
	double facspe=6.5;

	RANGE def={-9999,9999};
	RANGE x={0,100};
	//RANGE y={-5e-3,100e-3};
	RANGE y={-5e-3,160e-3};
    //RANGE q={-0.1,1};
    RANGE q={-0.01,2};
    RANGE t=x;
    charRANGE chR={
		x,
        y,
        q,
        def,
        def,
        def,
        t
	};//x,y,q,r,bl,blrms,t
	
    charRANGE chcut={
		def,
        def,
        def,
        def,
        def,
        def,
        def
	};//blrms,y 
	int rbt=128,rbU=16;
	int rb=2;
	int maxpk=10;
	double res=2;
	double sigmaA=2e-3;
	double sigmaq=2;
	double thrd=0.01;
    double fac=1.8;
	double fac2=0.4;
    double iter=4;
//*
//** end **
//

	//sprintf(buff,"/mnt/f/experiment/FTOF/SIPM/labtest/0108/BFP650_v2_s2/123.dat");
	ofstream op;
	//cout<<"Build your data file "<<endl;
	sprintf(buff,"%s.root",name);
    TFile *f1 = new TFile(buff,"READ");
    TTree *t1 = (TTree*)f1->Get("Pico");
	string path = f1->GetPath();
	path=path.substr(0,path.length()-7);
	sprintf(buff,"%s.dat",path.c_str());
	op.open(buff,ios::app);
	cout<<name<<endl;
	
	
	
    EVENT A;
    //double t0;
	double p[20];
	double Q[4];
	//t1->SetBranchAddress("MCP2_twentypercent_time",&t0);
	
    t1->SetBranchAddress("MCP2_CFDtime",A.CFD);
	t1->SetBranchAddress("MCP2_CFDfrac",p);
	t1->SetBranchAddress("MCP2_all_charge",Q);
	t1->SetBranchAddress("MCP2_rise_time",&A.rise);
	t1->SetBranchAddress("MCP2_global_maximum_y",&A.Amp);
	t1->SetBranchAddress("MCP2_global_maximum_x",&A.x);
	t1->SetBranchAddress("MCP2_baseline_rms",&A.blrms);
	t1->SetBranchAddress("MCP2_baseline_level",&A.bl);
	
	//TCut c_x = "MCP2_global_maximum_x>410&&MCP2_global_maximum_x<412.5";
	//TCut c_x = "MCP2_global_maximum_x>530&&MCP2_global_maximum_x<532";
	//TCut c_blrms = "MCP2_baseline_rms<0.1e-3";
    DrawMyfunc draw;
    LZWfunc lzw;

	TCanvas *c[7];
	for (int i=0;i<3;i++){
    	sprintf(buff,"C%d",i);
        c[i]= new TCanvas(buff,buff,800,600);
    	draw.SetPad(gPad,0.12,0.12,0.1,0.1);
	}

    TH1F *hq = new TH1F("hq",";charge (pC);Counts",1000,chR.q.L,chR.q.R);
	draw.Hist(hq,"charge (pC)","Counts",2);

	TH1F *ha = new TH1F("ha",";Amp(V);Counts",500,chR.y.L,chR.y.R);
	draw.Hist(ha,"Amp(V)","Counts",2);
	
	TH1F *hr = new TH1F("hr",";risetime (ns);Counts",200,chR.r.L,chR.r.R);
	draw.Hist(hr,"Risetime (ns)","Counts",2);
	
	TH1F *ht = new TH1F("ht","",200,chR.t.L,chR.t.R);
	draw.Hist(ht,"STR (ps)","Counts",2);

	TH1F *hbl = new TH1F("hbl","",400,chR.bl.L,chR.bl.R);
	draw.Hist(hbl,"baseline (V)","Counts",2);
	TGaxis::SetMaxDigits(3);
	
	TH1F *hblrms = new TH1F("hblrms","",800,chR.blrms.L,chR.blrms.R);
	draw.Hist(hblrms,"baseline RMS (V)","Counts",2);

	TH2F* hqy=new TH2F("hqy","",200,chR.y.L,chR.y.R,200,chR.q.L,chR.q.R);
    draw.Hist(hqy,"Amp (V)","Charge (pC)");

	TH1F *hx = new TH1F("hx",";Time (ns);Counts",2e3,chR.x.L,chR.x.R);
	draw.Hist(hx,"Time (ns)","Counts",2);

	c[0]->cd();

    //find the cut range;
	t1->Draw("MCP2_global_maximum_x>>hx");
	hx->GetXaxis()->SetRangeUser(5,80);
    chR.x.L=hx->GetBinCenter(hx->GetMaximumBin())-3;
    chR.x.R = hx->GetBinCenter(hx->GetMaximumBin())+3;
    TF1* gx=lzw.gausfit(hx,1,1.5,&chR.x);
//return;

    chcut.x.L = gx->GetParameter(1)-2.5*gx->GetParameter(2);
	chcut.x.R = gx->GetParameter(1)+2.5*gx->GetParameter(2);
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
				hq->Fill(Q[qN]);
				ha->Fill(A.Amp);
				//cout<<Q[1]<<endl;
				hr->Fill(A.rise);
				hbl->Fill(A.bl);
				hblrms->Fill(A.blrms);
				hqy->Fill(A.Amp,Q[qN]);
				
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

/*
//	amplitude spectrum
	sprintf(buff,"%s_amp",name);
    lzw.Set_name(buff);
	c[1]->cd();
	ha->Draw();
	c[1]->SetLogy();
	TF1* fa=lzw.mcpSPfit(ha,3,chR.y,20);
	double SPEAmp=(fa->GetParameter(4)-fa->GetParameter(1))*1e3;
	cout<<"SPEAmp="<<SPEAmp<<endl;
	sprintf(buff,"SPE Amp=%0.2f mV",SPEAmp);
	TLatex *l1=draw.Latex(0.4,0.2,buff);
	l1->Draw();
	sprintf(buff,"%s_amp.png",name);
    c[1]->SaveAs(buff);
	return;
*/

/*	TF1* fa=lzw.fpeaksfit(ha,maxpk,res,sigmaA,thrd);
	gausPAR* gPa = lzw.Get_PeakPar();
	double speA=(gPa+1)->m-gPa->m;
	chcut.y.L=gPa->m-2*gPa->s;
	chcut.y.R=gPa->m+2*gPa->s;
	cout<<"chcut y left ="<<chcut.y.L<<", chcut y right ="<<chcut.y.R<<endl;
*/
//	Charge spectrum
	c[2]->cd();
	//c[2]->SetLogy();
	sprintf(buff,"%s_charge",name);
    lzw.Set_name(buff);
	//hq->Draw();
	//return;
	//TF1* fq=lzw.SPSfit(hq,4,chR.q,facspe);
	TF1* fq=lzw.mcpSPfit(hq,rbu,chR.q,facspe);
	double Gain=(fq->GetParameter(4)-fq->GetParameter(1))*1e-12/1.6e-19;
	cout<<"Gain="<<Gain<<endl;
	sprintf(buff,"Gain=%0.2e",Gain);
	TLatex *l=draw.Latex(0.4,0.2,buff);
	l->Draw();
	sprintf(buff,"%s_charge.png",name);
    c[2]->SaveAs(buff);

	//gausPAR* gPq = lzw.Get_PeakPar();
	//double speq=(gPq+1)->m-gPq->m;
	//double G=speq*1e-12/1.6e-19;

	return;


//	risetime spectrum
	c[4]->cd();
	TF1* r1=lzw.twoguasfit(hr,fac2,rb,chR.r);
	//return;
	double risetime=r1->GetParameter(1);
	sprintf(buff,"%s_risetime.png",name);
    c[4]->SaveAs(buff);

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

/*	
//	time resolution spectrum
	lzw.Set_charcut(chcut);
    TF1* tt=lzw.CH1Correction(t1,&A,&t0,t,U,rbU,2*rbt,fac,iter);
	double STR=tt->GetParameter(2);
	double STRerr=tt->GetParError(2);
	
	
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