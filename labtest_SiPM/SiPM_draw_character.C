void SiPM_draw_character(const char* rootname=""){
	
	TH2F* DrawMyHist2d(const char* name,double x1,double x2,double y1,double y2, char *xtitle, char *ytitle, Color_t LColor=1, Width_t LWidth=1.5);
	TH1F* DrawMyHist(const char* name,double x1,double x2, char *xtitle, char *ytitle, Color_t LColor=1, Width_t LWidth=1.5 );
	
	char buff[1024];
	char name[1024];
	
	sprintf(name,"%s",rootname);
	sprintf(buff,"%s.root",name);
	TFile *f1 = new TFile(buff,"update");
	//TTree *t1 = (TTree*)f1->Get("SiPM_MCP");
	//TTree *t2 = (TTree*)f1->Get("SiPM");
	TTree *t3 = (TTree*)f1->Get("Pico");
	
	
	double qL=0;
	double qR=1e3;
	double yL=0;
	double yR=1.2;
	double riseL=0;
	double riseR=5;
	double xL=250;
	double xR=270;
	double blL=-20e-3;
	double blR=20e-3;
	double blrmsL=0;
	double blrmsR=20e-3;
	
	string chname[2]={"A2","B3"};
	
	for (int i=0;i<2;i++)
	{
	TGaxis::SetMaxDigits(3);	
	TH2F* hqy=DrawMyHist2d("hqy",yL,yR,qL,qR,"Amp (V)","Charge (pC)");
	sprintf(buff,"MCP%d_all_charge:MCP%d_global_maximum_y>>hqy",i+1);
	t3->Draw(buff,"","colz");
	sprintf(buff,"%s_Qvsy.png",chname[i].data());
	gPad->SaveAs(buff);
	
	TH2F* hqr=DrawMyHist2d("hqr",riseL,riseR,qL,qR,"risetime (ns)","Charge (pC)");
	sprintf(buff,"MCP%d_all_charge:MCP%d_rise_time>>hqr",i+1);
	t3->Draw(buff,"","colz");
	sprintf(buff,"%s_Qvsrise.png",chname[i].data());
	gPad->SaveAs(buff);
	
	
	TH2F* hqx=DrawMyHist2d("hqx",xL,xR,qL,qR,"time (ns)","Charge (pC)");
	sprintf(buff,"MCP%d_all_charge:MCP%d_global_maximum_x>>hqx",i+1);
	t3->Draw(buff,"","colz");
	sprintf(buff,"%s_Qvsx.png",chname[i].data());
	gPad->SaveAs(buff);

	
	TH1F* hbl=DrawMyHist("hbl",blL,blR,"baseline level(V)","Counts");
	sprintf(buff,"MCP%d_baseline_level>>hbl",i+1);
	t3->Draw(buff,"","colz");
	sprintf(buff,"%s_bl.png",chname[i].data());
	gPad->SaveAs(buff);
	
	
	TH1F* hblrms=DrawMyHist("hblrms",blrmsL,blrmsR,"baseline RMS(V)","Counts");
	sprintf(buff,"MCP%d_baseline_rms>>hblrms",i+1);
	t3->Draw(buff,"","colz");
	sprintf(buff,"%s_blrms.png",chname[i].data());
	gPad->SaveAs(buff);

	}
	
	/*TH2F* hqy = new TH2F("hqy","",2e3,yL,yR,200,qL,qR);
	TH2F* hqr = new TH2F("hqr","",2e3,,10);
	TH1F* hrise = new TH1F("hy","",2e3,-10,10);
	TH1F* hy = new TH1F("hy","",2e3,-10,10);
	TH1F* hy = new TH1F("hy","",2e3,-10,10);
	TH1F* hy = new TH1F("hy","",2e3,-10,10);*/
	
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

void SetMyPad(TVirtualPad *pad,double left, double right, double top, double bottom){
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

TF1* gausfit(TH1* hU,int rbU,double U_RL, double U_RR){
	
	
	
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
		TF1 *fit = new TF1("fit","gaus",*tRL,*tRR);        
		h->GetXaxis()->SetRangeUser(*tRL,*tRR);
		h->Fit(fit);
        double mean = fit->GetParameter(1);
        double sigma = TMath::Abs(fit->GetParameter(2));
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