#include <string>
#include <time.h>
#include <TString.h>
using namespace std;

void SiPM_SPE_Character(const char *rootname="")
{
	void DrawMyHist1(TH1F *datahist,const char *xtitle,const char *ytitle, Color_t LColor=1, Width_t LWidth=3 , Color_t TitleColor=1);
	TH2F* DrawMyHist2d(const char* name,double x1,double x2,double y1,double y2, char *xtitle, char *ytitle, Color_t LColor=1, Width_t LWidth=1.5);
	TH1F* DrawMyHist(const char* name,double x1,double x2, char *xtitle, char *ytitle, Color_t LColor=1, Width_t LWidth=1.5 );
	TF1* gausfit(TH1* hU,int rbU,float U_RL, float U_RR);
    void SetMyPad(TVirtualPad *pad,float left, float right, float top, float bottom);
	TF1* twoguasfit(TH1 *ht,double* tRL,double* tRR,double fac=0.4, int rbt=1);
	string getTime();
	
	
	char name[1024];
	sprintf(name,"%s",rootname);
	
    char buff[1024];

/*==========>*******************************<========	
******************** Cut Range **********************
*****************************************************
**==========>*******************************<========*/
	double xL = 100; //base cut
	double xR = 120;
	/*double blrmsL = 0.e-3;
	double blrmsR = 0.1e-3;
	double blL = -0.5e-3;
	double blR = 0.1e-3;
	*/
	double blrmsL = 0.e-3;
	double blrmsR = 2.e-3;
	double blL = -20.5e-3;
	//double blR = 0;
	double blR = 20.05e-3;
	
	
	double xRL = 0;
	double xRR = 200;
	
	
	double tRL=26; //hist time range
	double tRR=35;
	
	double rRL=0;//hist risetime range
	double rRR=3;
	
	double bRL = blL-1e-3;  //hist baseline range
	double bRR = blR+1e-3;
	
	double brmsRL = blrmsL-1e-3;  //hist baseline range
	double brmsRR = blrmsR+1e-3;

	double aRL=-1e-3,aRR=182e-3; //find Amp peak range
	double alimit1=22e-3;
	double alimit2=38e-3;

	double qRL=-0.5,qRR=80;	//find Charge peak range
	double qlimit1=8;
	double qlimit2=14;
		
	gStyle->SetOptFit(1111);
	sprintf(buff,"%s.dat",name);
	ofstream output(buff,ios::app);
	
	sprintf(buff,"%s.root",name);
    TFile *f1 = new TFile(buff,"READ");
    TTree *t1 = (TTree*)f1->Get("Pico");
	
	
	
	double T2,T1,Q1,Q2;
	double x,y,blrms,rise,bl;
	t1->SetBranchAddress("MCP2_twentypercent_time",&T2);
	t1->SetBranchAddress("MCP1_twentypercent_time",&T1);
	t1->SetBranchAddress("MCP2_all_charge",&Q2);
	t1->SetBranchAddress("MCP1_all_charge",&Q1);
	t1->SetBranchAddress("MCP2_rise_time",&rise);
	t1->SetBranchAddress("MCP2_global_maximum_y",&y);
	t1->SetBranchAddress("MCP2_global_maximum_x",&x);
	t1->SetBranchAddress("MCP2_baseline_rms",&blrms);
	t1->SetBranchAddress("MCP2_baseline_level",&bl);
	
	//TCut c_x = "MCP2_global_maximum_x>410&&MCP2_global_maximum_x<412.5";
	//TCut c_x = "MCP2_global_maximum_x>530&&MCP2_global_maximum_x<532";
	//TCut c_blrms = "MCP2_baseline_rms<0.1e-3";
	TCanvas *c[7];
	for (int i=0;i<7;i++){
	sprintf(buff,"C%d",i);
    c[i]= new TCanvas(buff,buff,800,600);
	SetMyPad(gPad,0.12,0.1,0.1,0.12);
	}

    TH1F *hq = new TH1F("hq",";charge (pC);Counts",200,qRL,qRR);
	DrawMyHist1(hq,"charge (pC)","Counts",kBlue,2);

	TH1F *ha = new TH1F("ha",";Amp(V);Counts",200,aRL,aRR);
	DrawMyHist1(ha,"Amp(V)","Counts",kBlue,2);
	
	TH1F *hr = new TH1F("hr",";risetime (ns);Counts",200,rRL*0.8,rRR*1.2);
	DrawMyHist1(hr,"Risetime (ns)","Counts",kBlue,2);
	
	TH1F *ht = new TH1F("ht","",200,tRL,tRR);
	DrawMyHist1(ht,"STR (ps)","Counts",kBlue,2);
	
	TH1F *hbl = new TH1F("hbl","",400,bRL,bRR);
	DrawMyHist1(hbl,"baseline (V)","Counts",kBlue,2);
	TGaxis::SetMaxDigits(3);
	
	TH1F *hblrms = new TH1F("hblrms","",800,brmsRL,brmsRR);
	DrawMyHist1(hblrms,"baseline RMS (V)","Counts",kBlue,2);

	TH2F* hqy=DrawMyHist2d("hqy",aRL,aRR,qRL,qRR,"Amp (V)","Charge (pC)");

	TH1F *hx = new TH1F("hx",";Time (ns);Counts",2e3,xRL,xRR);
	DrawMyHist1(hx,"Time (ns)","Counts",kBlue,2);
	c[1]->cd();
	t1->Draw("MCP2_global_maximum_x>>hx");
	//return;
    xRL = hx->GetBinCenter(hx->GetMaximumBin())-10;
    xRR = hx->GetBinCenter(hx->GetMaximumBin())+10;
	TF1* gx=twoguasfit(hx,&xRL,&xRR,0.4,1);
	xL = gx->GetParameter(1)-2.5*gx->GetParameter(2);
	xR = gx->GetParameter(1)+2.5*gx->GetParameter(2);
	cout<<"global maximum x range = "<<xL<<"\t"<<xR<<endl;
	
	//gPad->Modified(); 
	//gPad->Update(); 
 
 
 
 int N=t1->GetEntries();
	for(int i=0;i<N;i++){
		t1->GetEntry(i);
		//if(Q2>qlimit1&&Q2<qlimit2)
			if (x>xL&&x<xR&&blrms>blrmsL&&blrms<blrmsR&&bl<blR&&bl>blL&&Q2>0.5)
			{
				ha->Fill(y);
				hq->Fill(Q2);
				hr->Fill(rise);
				hbl->Fill(bl);
				hblrms->Fill(blrms);
				hqy->Fill(y,Q2);
				if(y<alimit1&&y>aRL)
				//if(y<alimit2&&y>alimit1)
				//if(y<1.8e-3&&y>1.3e-3)
					ht->Fill(T2-T1);
			}
	}
	//t1->Draw("MCP2_global_maximum_y>>ha",c_x&&c_blrms);	
	c[1]->cd();
	cout<<"program running normly"<<endl;	
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
	/*	t1->Draw("MCP2_twentypercent_time-MCP1_twentypercent_time>>ht","MCP2_all_charge>0.2&&MCP2_all_charge<8");
	*/
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
	
}
TH2F* DrawMyHist2d(const char* name,double x1,double x2,double y1,double y2, char *xtitle, char *ytitle, Color_t LColor=1, Width_t LWidth=1.5)
{	
	TH2F *datahist= new TH2F(name,"",200,x1,x2,200,y1,y2);
	TCanvas *c1 =new TCanvas("c1","c1",800,600);
	gPad->SetMargin(0.14,0.1,0.14,0.1);
	datahist->SetLineColor( LColor );
    datahist->SetLineWidth( LWidth );
	datahist->GetXaxis()->SetTitle( xtitle);
     datahist->GetYaxis()->SetTitle( ytitle);
     datahist->GetXaxis()->SetAxisColor(1);
     datahist->GetYaxis()->SetAxisColor(1);
     datahist->GetXaxis()->SetLabelColor(1);
     datahist->GetYaxis()->SetLabelColor(1);
     datahist->GetXaxis()->SetLabelFont( 42 );
     datahist->GetYaxis()->SetLabelFont( 42 );
     datahist->GetXaxis()->SetLabelSize( 0.06 );
     datahist->GetYaxis()->SetLabelSize( 0.06 );
     datahist->GetXaxis()->SetLabelOffset( 0.01 );
	 datahist->GetXaxis()->SetTitleFont( 42 );
     datahist->GetYaxis()->SetTitleFont( 42 );
	 
	 Color_t TitleColor=1;
     datahist->GetXaxis()->SetTitleColor( TitleColor);
     datahist->GetYaxis()->SetTitleColor( TitleColor );
     datahist->GetXaxis()->SetTitleSize(0.06);
     datahist->GetYaxis()->SetTitleSize(0.06);
     datahist->GetXaxis()->SetTitleOffset(1.0);
     datahist->GetYaxis()->SetTitleOffset(1.0);
	 //datahist->GetXaxis()->SetBorderSize(5);
     datahist->GetXaxis()->SetNdivisions(510);
     datahist->GetYaxis()->SetNdivisions(510);
     datahist->GetXaxis()->CenterTitle();
     datahist->GetYaxis()->CenterTitle();
	 return datahist;
	 }

TH1F* DrawMyHist(const char* name,double x1,double x2, char *xtitle, char *ytitle, Color_t LColor=1, Width_t LWidth=1.5 ){
	
	TH1F *datahist= new TH1F(name,"",200,x1,x2);
	TCanvas *c2 =new TCanvas("c2","c2",800,600);
	gPad->SetMargin(0.14,0.1,0.14,0.1);
     datahist->SetLineColor( LColor );
     datahist->SetLineWidth( LWidth );
	datahist->GetXaxis()->SetTitle( xtitle);
     datahist->GetYaxis()->SetTitle( ytitle);
     datahist->GetXaxis()->SetAxisColor(1);
     datahist->GetYaxis()->SetAxisColor(1);
     datahist->GetXaxis()->SetLabelColor(1);
     datahist->GetYaxis()->SetLabelColor(1);
     datahist->GetXaxis()->SetLabelFont( 42 );
     datahist->GetYaxis()->SetLabelFont( 42 );
     datahist->GetXaxis()->SetLabelSize( 0.06 );
     datahist->GetYaxis()->SetLabelSize( 0.06 );
     datahist->GetXaxis()->SetLabelOffset( 0.01 );
     datahist->GetYaxis()->SetLabelOffset( 0.01 );
     datahist->GetXaxis()->SetTitleFont( 42 );
     datahist->GetYaxis()->SetTitleFont( 42 );
	 Color_t TitleColor=1;
     datahist->GetXaxis()->SetTitleColor( TitleColor);
     datahist->GetYaxis()->SetTitleColor( TitleColor );
     datahist->GetXaxis()->SetTitleSize(0.06);
     datahist->GetYaxis()->SetTitleSize(0.06);
     datahist->GetXaxis()->SetTitleOffset(1.0);
     datahist->GetYaxis()->SetTitleOffset(1.1);
     //datahist->GetXaxis()->SetBorderSize(5);
     datahist->GetXaxis()->SetNdivisions(510);
     datahist->GetYaxis()->SetNdivisions(510);
     datahist->GetXaxis()->CenterTitle();
     datahist->GetYaxis()->CenterTitle();
	 return datahist;
}
void DrawMyHist1(TH1F *datahist,const char *xtitle,const char *ytitle, Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1){
     datahist->SetLineColor( LColor );
     datahist->SetLineWidth( LWidth );
     datahist->GetXaxis()->SetTitle( xtitle);
     datahist->GetYaxis()->SetTitle( ytitle);
     datahist->GetXaxis()->SetAxisColor(1);
     datahist->GetYaxis()->SetAxisColor(1);
     datahist->GetXaxis()->SetLabelColor(1);
     datahist->GetYaxis()->SetLabelColor(1);
     datahist->GetXaxis()->SetLabelFont( 42 );
     datahist->GetYaxis()->SetLabelFont( 42 );
     datahist->GetXaxis()->SetLabelSize( 0.06 );
     datahist->GetYaxis()->SetLabelSize( 0.06 );
     datahist->GetXaxis()->SetLabelOffset( 0.01 );
     datahist->GetYaxis()->SetLabelOffset( 0.01 );
     datahist->GetXaxis()->SetTitleFont( 42 );
     datahist->GetYaxis()->SetTitleFont( 42 );
     datahist->GetXaxis()->SetTitleColor( TitleColor);
     datahist->GetYaxis()->SetTitleColor( TitleColor );
     datahist->GetXaxis()->SetTitleSize(0.06);
     datahist->GetYaxis()->SetTitleSize(0.06);
     datahist->GetXaxis()->SetTitleOffset(1.0);
     datahist->GetYaxis()->SetTitleOffset(1.0);
     //datahist->GetXaxis()->SetBorderSize(5);
     datahist->GetXaxis()->SetNdivisions(510);
     datahist->GetYaxis()->SetNdivisions(510);
     datahist->GetXaxis()->CenterTitle();
     datahist->GetYaxis()->CenterTitle();
}

void SetMyPad(TVirtualPad *pad,float left, float right, float top, float bottom){
  pad->SetFillColor(10);
  pad->SetBorderMode(0);
  pad->SetBorderSize(0);
  pad->SetFrameFillColor(10);
  //pad->SetFrameFillStyle(3003);
  pad->SetFrameBorderMode(0);
  pad->SetFrameBorderSize(0);
  pad->SetLeftMargin(left);
  pad->SetRightMargin(right);
  pad->SetTopMargin(top);
  pad->SetBottomMargin(bottom);
}


TF1* gausfit(TH1* hU,int rbU,float U_RL, float U_RR){
	
	
	
	hU->Draw();
	hU->Rebin(rbU);
	TF1 *fitU = new TF1("fitU","gaus",U_RL,U_RR);
	fitU->SetParameter(1,hU->GetMean());
	hU->Fit(fitU,"R");
	/*
	U_RL = fitU->GetParameter(1)-5*fitU->GetParameter(2);
	U_RR = fitU->GetParameter(1)+5*fitU->GetParameter(2);

	hU->GetXaxis()->SetRangeUser(U_RL,U_RR);
	hU->Fit(fitU);
	*/
	return fitU;

}

TF1* twoguasfit(TH1 *ht,double* tRL,double* tRR,double fac=0.4, int rbt=1){
        //First fit for ensuring the rangement of histgram;
		TH1* h = (TH1*)ht->Clone();
		h->Rebin(rbt);
        double mean = h->GetBinCenter(h->GetMaximumBin());
        double sigma = h->GetRMS();
		TF1 *fit = new TF1("fit","gaus",*tRL,*tRR);        
		h->GetXaxis()->SetRangeUser(*tRL,*tRR);
		h->Fit(fit);
        mean = fit->GetParameter(1);
        sigma = TMath::Abs(fit->GetParameter(2));
		if(*tRL<mean-5*sigma||sigma>1)
		{*tRL = mean-5*sigma;
		*tRR = mean+5*sigma;}
		cout<<h->GetName()<<"\t"<<*tRL<<"\t"<<*tRR<<endl;
       
		h->GetXaxis()->SetRangeUser(*tRL,*tRR);
		
        TF1 *fit2 = new TF1("fit2","gaus(0)+gaus(3)",*tRL,*tRR);
		fit2->SetParNames("C_{TR}","#mu_{TR}","#sigma_{TR}","C_{nz}","#mu_{nz}","#sigma_{nz}");
		fit2->SetParameter(1,mean);
		fit2->SetParameter(2,sigma);
		fit2->SetParLimits(3,0,fit->GetParameter(0)*fac);
		fit2->SetParameter(4,mean);
		fit2->SetParameter(5,1.5*sigma);
		
		h->Fit(fit2,"","",mean-5*sigma,mean+5*sigma);
		TF1 *fit_tr = new TF1("fit_tr","gaus",*tRL,*tRR);
		fit_tr->SetParameter(0,fit2->GetParameter(0));
		fit_tr->SetParameter(1,fit2->GetParameter(1));
		fit_tr->SetParameter(2,fit2->GetParameter(2));
		TF1 *fit_nz = new TF1("fit_nz","gaus",*tRL,*tRR);
		fit_nz->SetParameter(0,fit2->GetParameter(3));
		fit_nz->SetParameter(1,fit2->GetParameter(4));
		fit_nz->SetParameter(2,fit2->GetParameter(5));
		fit_tr->SetLineColor(3);
		fit_tr->SetLineStyle(7);
		fit_nz->SetLineColor(3);
		fit_nz->SetLineStyle(7);
		fit_tr->Draw("same");
		fit_nz->Draw("same");
		
		if(fit2->GetParameter(3)>=fit2->GetParameter(0)&&TMath::Abs(fit2->GetParameter(5))>=0.05)
		{
			mean = fit2->GetParameter(4);
			sigma = TMath::Abs(fit2->GetParameter(5));
			}
		else if(TMath::Abs(fit2->GetParameter(2))>=0.05){
			mean = fit2->GetParameter(1);
			sigma = TMath::Abs(fit2->GetParameter(2));
		}
		
		if(*tRL<mean-10*sigma)
		{
			*tRL = mean-10*sigma;
			}
		if(*tRR>mean+10*sigma)
		{
			*tRR = mean+10*sigma;
			}
		
		cout<<h->GetName()<<"\t"<<*tRL<<"\t"<<*tRR<<endl;
        h->GetXaxis()->SetRangeUser(*tRL,*tRR);
		return fit2;
}


string getTime()
{
    time_t timep;
    time (&timep);
    char tmp[64];
    strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep) );
    return tmp;
}

