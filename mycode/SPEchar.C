#include <string>
#include <time.h>
#include <TString.h>

#include "Include/DrawMyfunc.h"
#include "Include/LZWfunc.h"

using namespace std;
void SPEchar(const char *rootname="FEboardA_ch2"){

	gStyle->SetOptFit(1111);


	
	
	char name[1024];
	sprintf(name,"%s",rootname);
    char buff[1024];
  
    charRANGE chR={{0,200},{-1e-3,182e-3},{-0.5,80},{0,3},{-22e-3,22e-3},{-1e-3,3e-3},{20,40}};//x,y,q,r,bl,blrms,t
    charRANGE chcut={{-1e4,1e4},{22e-3,38e-3},{-1e4,1e4},{-1e4,1e4},{-1e4,1e4},{0e-3,2e-3},{-1e4,1e4}};//blrms,y 
	int rbt=128,rbU=16;
    double fac=1.8;
    double iter=4;

	sprintf(buff,"%s.dat",name);
	ofstream output(buff,ios::app);
	
	sprintf(buff,"%s.root",name);
    TFile *f1 = new TFile(buff,"READ");
    TTree *t1 = (TTree*)f1->Get("Pico");
	
	
	
    EVENT A;
    double t0;
	t1->SetBranchAddress("MCP1_twentypercent_time",&t0);
	
    t1->SetBranchAddress("MCP2_twentypercent_time",&A.time);
	t1->SetBranchAddress("MCP2_all_charge",&A.Q);
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
	for (int i=1;i<2;i++){
    	sprintf(buff,"C%d",i);
        c[i]= new TCanvas(buff,buff,800,600);
    	draw.SetPad(gPad,0.12,0.12,0.1,0.1);
	}

    TH1F *hq = new TH1F("hq",";charge (pC);Counts",200,chR.q.L,chR.q.R);
	draw.Hist(hq,"charge (pC)","Counts",2);

	TH1F *ha = new TH1F("ha",";Amp(V);Counts",400,chR.y.L,chR.y.R);
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

	c[1]->cd();

    //find the cut range;
	t1->Draw("MCP2_global_maximum_x>>hx");
    chR.x.L=hx->GetBinCenter(hx->GetMaximumBin())-10;
    chR.x.R = hx->GetBinCenter(hx->GetMaximumBin())+10;
    TF1* gx=lzw.gausfit(hx,1,1.5,&chR.x);

    chcut.x.L = gx->GetParameter(1)-2.5*gx->GetParameter(2);
	chcut.x.R = gx->GetParameter(1)+2.5*gx->GetParameter(2);
	cout<<"global maximum x cut range = "<<chcut.x.L<<"\t"<<chcut.x.R<<endl;
//return;



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
            A.bl>chcut.bl.L&&
            A.Q>chcut.q.L)
			{
				ha->Fill(A.Amp);
				hq->Fill(A.Q);
				hr->Fill(A.rise);
				hbl->Fill(A.bl);
				hblrms->Fill(A.blrms);
				hqy->Fill(A.Amp,A.Q);
				
			}
	}

    lzw.Set_name(name);
    lzw.Set_charcut(chcut);
    RANGE t=chR.t;
    RANGE U={0,0.05};
    //lzw.CH1Correction(t1,&A,&t0,t,U,rbU,rbt,fac,iter);


	//t1->Draw("MCP2_global_maximum_y>>ha",c_x&&c_blrms);	
	c[1]->cd();
	cout<<"program running normly"<<endl;	
	TSpectrum *s = new TSpectrum(20,10);
	int nfound = s->Search(ha,2e-3,"",0.05);
	printf("Found %d candidate peaks to fit\n",nfound);
	ha->Draw();
	return;
	double *xpeaks;
	xpeaks = s->GetPositionX();

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


	
	string time = getTime();
    output<< time << endl;
	output<<"The spe amplitude (mV) = "<<pe1A*1e3<<endl;
	output<<"The spe charge (pC) = "<<speq<<endl;
	output<<"The Gain = "<<G<<endl;
	output<<"The STR (ps) = "<<STR*1e3<<"\t"<<STRerr*1e3<<endl;
	output<<"The ped amp mean (mV) = "<<pedmean<<endl;
	output<<"The ped amp sigma (mV) = "<<pedmeansigma<<endl;
	output<<"The ped amp RMS (mV) = "<<pedRMS<<endl;
	output<<"The risetime (ns)= "<<risetime<<endl;
	output<<"\n"<<endl;
	*/
}