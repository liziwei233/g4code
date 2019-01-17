struct EVENT{
	double p20,A,Q,rise,blrms;
};


void SiPM_MultiCH_Calculate_TR(const char* rootname=""){
	void twoguasfit(TH1 *ht,double* tRL,double* tRR,double fac=0.4,int rbt=1);
	TF1* profilefit(TH2* Rt,double rbU,double rbt,double tRL,double tRR,double URL,double URR,char* name);

	TF1* CH1ACorrection(TTree* t1,EVENT* A,EVENT* B,EVENT* MCP,int rbU,int rbt,char* name,double fac,int iter);
	TF1* CH1BCorrection(TTree* t1,EVENT* A,EVENT* B,EVENT* MCP,int rbU,int rbt,char* name,double fac,int iter);

	TF1* CH2Correction(TTree* t1,EVENT* A,EVENT* B,EVENT* MCP,int rbU,int rbt,char* name,double fac,int iter);
	TF1* CH3Correction(TTree* t1,EVENT* A,EVENT* B,EVENT* MCP,int rbU,int rbt,char* name,double fac,int iter);

	
	//void draw_TR(TTree *t1,double** A1,double** A2,double** B1,double** B2,double tRL,double tRR,double URL,double URR,double RRL,double RRR,char* name);
	//void draw_TR(TTree *t1,double** A1,double** A2,double** B1,double** B2,double tRL,double tRR,double URL,double URR,double RRL,double RRR,double mcptRL,double mcptRR,double mcpURL,double mcpURR,double mcpRRL,double mcpRRR,char* name);

	gStyle->SetOptFit(1111);

	char buff[1024];
	char name[1024];

	int rbU=2;
	int rbt=1;




	//ofstream output("timeres.dat");



	//TCanvas *c1 = new TCanvas("c1","c1",800,600);
	//TCanvas *c2 = new TCanvas("c2","c2",800,600);
	//TCanvas *c3 = new TCanvas("c3","c3",800,600);

	//TCanvas *c5 = new TCanvas("c4","c5",800,600);

	//sprintf(name,"%dpe_57v",npe[s]);
	//sprintf(name,"OSCch3_A2vs_OSCch4_B3_Vsipm_Sep");
	for (int i = 0;i<6;i++){
		sprintf(name,"scanvoltage%dv",55+i);
	
	//sprintf(name,"%s",rootname);
	sprintf(buff,"%s.root",name);
	TFile *f1 = new TFile(buff,"update");
	//TTree *t1 = (TTree*)f1->Get("SiPM_MCP");
	//TTree *t2 = (TTree*)f1->Get("SiPM");
	TTree *t3 = (TTree*)f1->Get("Pico");


	EVENT A,B,MCP; 
	//double A1rise,A2rise,B1rise;
	/*	
		t1->SetBranchAddress("A1",&A.p20);
		t1->SetBranchAddress("A2",&MCP.p20);
		t1->SetBranchAddress("B1",&B.p20);


		t1->SetBranchAddress("A1Q",&A.Q);
		t1->SetBranchAddress("A2Q",&MCP.Q);
		t1->SetBranchAddress("B1Q",&B.Q);

		t1->SetBranchAddress("A1y",&A.A);
		t1->SetBranchAddress("A2y",&MCP.A);
		t1->SetBranchAddress("B1y",&B.A);

		t2->SetBranchAddress("A1",&A.p20);
		t2->SetBranchAddress("A2",&MCP.p20);
		t2->SetBranchAddress("B1",&B.p20);


		t2->SetBranchAddress("A1Q",&A.Q);
		t2->SetBranchAddress("A2Q",&MCP.Q);
		t2->SetBranchAddress("B1Q",&B.Q);

		t2->SetBranchAddress("A1y",&A.A);
		t2->SetBranchAddress("A2y",&MCP.A);
		t2->SetBranchAddress("B1y",&B.A);
		*/



	t3->SetBranchAddress("MCP1_twentypercent_time",&A.p20);
	t3->SetBranchAddress("MCP2_twentypercent_time",&B.p20);
	//t3->SetBranchAddress("MCP3_twentypercent_time",&B2.p20);


	t3->SetBranchAddress("MCP1_all_charge",&A.Q);
	t3->SetBranchAddress("MCP2_all_charge",&B.Q);
	//t3->SetBranchAddress("MCP3_all_charge",&B2.Q);

	t3->SetBranchAddress("MCP1_global_maximum_y",&A.A);
	t3->SetBranchAddress("MCP2_global_maximum_y",&B.A);
	//t3->SetBranchAddress("MCP3_global_maximum_y",&B2.A);

	t3->SetBranchAddress("MCP1_rise_time",&A.rise);
	t3->SetBranchAddress("MCP2_rise_time",&B.rise);
	//t3->SetBranchAddress("MCP3_rise_time",&B2.rise);
	
	t3->SetBranchAddress("MCP1_baseline_rms",&A.blrms);
	t3->SetBranchAddress("MCP2_baseline_rms",&B.blrms);




	TH1F* hA = new TH1F("hA","Time Resolution;T (ns);Counts",2e3,-10,10);
	TH1F* hB = new TH1F("hB","Time Resolution;T (ns);Counts",2e3,-10,10);
	TH1F* hBB = new TH1F("hBB","Time Resolution;T (ns);Counts",2e3,-10,10);
	TH1F* hAB = new TH1F("hAB","Time Resolution;T (ns);Counts",1.6e3,-10,10);	


	/*		EVENT* p1;
			p1 = &A;
			int N = t1->GetEntries();
			cout<< "Total entries is :"<<N<<endl;
			for(int i=0; i<1; i++)
			{
			t1->GetEntry(i);


			cout<<(*p1).Q<<endl;}	
			return;		
			*/		
	//>>>>------------TA correction LOOP-------<<<<//
	double thrd[4]={700,800,900,1.6e3};
	double fac = 0.4;
	for (int i=0;i<1;i++){

	sprintf(buff,"B1B2_B1unsat%s",name);
	//CH3Correction(t3,&B1,&B2,&MCP,2*rbU,4*rbt,buff,fac,5);
	
	sprintf(buff,"A-B_%s",name);
	CH2Correction(t3,&A,&B,&A,2*rbU,4*rbt,buff,fac,5);
	
	sprintf(buff,"T0B2_B1B2unsat%s",name);
	//cout<<"<====process check====>"<<endl;
	//CH1ACorrection(t3,&B1,&B2,&MCP,2*rbU,4*rbt,buff,fac,5);

	sprintf(buff,"T0B2_B1unsat%s",name);
	//CH1ACorrection(t3,&B1,&B2,&MCP,2*rbU,4*rbt,buff,fac,5);

	}
	}
	return;
	/*	
		double thrd[4]={700,800,900,1.6e3};
		for (int i=0;i<1;i++){
		sprintf(buff,"T0B1_B1trigger_B1saturated_%s",name);
		TF1* A2A1_cor = CH1Correction(t1,&B1,&B1y,&B1Q,&A2,&A2y,&A2Q,&A1,1*rbU,4*rbt,buff,thrd[i],4);
		}
		return;

	//if(!t1->GetBranchStatus(buff)){


	double thrd[4]={700,800,900,1.6e3};
	for (int i=0;i<1;i++){
	sprintf(buff,"T0A1_A1trigger_%s",name);
	TF1* A2A1_cor = CH1Correction(t1,&A1,&A1y,&A1Q,&A2,&A2y,&A2Q,&B1Q,1*rbU,4*rbt,buff,thrd[i],4);
	}
	return;

	double thrd[4]={700,800,900,1.6e3};
	for (int i=0;i<1;i++){
	sprintf(buff,"withMCP_A1B1_%s",name);
	TF1* A2A1_cor = CH2Correction_noT0(t1,&A1,&A1y,&A1Q,&B1,&B1y,&B1Q,1*rbU,4*rbt,buff,thrd[i],5);
	}
	return;

	double thrd[4]={700,800,900,1.6e3};
	for (int i=0;i<1;i++){
	sprintf(buff,"T0A1B1_%s",name);
	TF1* A2A1_cor = CH2Correction(t1,&A1,&A1y,&A1Q,&B1,&B1y,&B1Q,&A2,&A2y,&A2Q,1*rbU,4*rbt,buff,thrd[i],4);
	}
	return;

	double thrd[4]={700,800,900,1.6e3};
	for (int i=0;i<4;i++){
	sprintf(buff,"A2A1_%s",name);
	TF1* A2A1_cor = CH1Correction(t1,&A1,&A1y,&A1Q,&A2,&A2y,&A2Q,1*rbU,4*rbt,buff,thrd[i],4);
	}
	return;
	double thrd[4]={700,800,900,1.6e3};
	for (int i=0;i<4;i++){
	sprintf(buff,"A2B1_%s",name);
	TF1* A2A1_cor = CH1Correction(t1,&B1,&B1y,&B1Q,&A2,&A2y,&A2Q,1*rbU,4*rbt,buff,thrd[i],4);
	}
	return;

*/	
	//


}

TF1* gausfit(TH1* hU,int rbU,float* U_RL, float* U_RR){



	hU->Draw();
	hU->Rebin(rbU);
	TF1 *fitU = new TF1("fitU","gaus",*U_RL,*U_RR);
	fitU->SetParameter(1,hU->GetMean());
	hU->Fit(fitU,"R");
	*U_RL = fitU->GetParameter(1)-5*fitU->GetParameter(2);
	*U_RR = fitU->GetParameter(1)+5*fitU->GetParameter(2);

	hU->GetXaxis()->SetRangeUser(*U_RL,*U_RR);
	hU->Fit(fitU);
	return fitU;

}

void twoguasfit(TH1 *ht,double* tRL,double* tRR,double fac=0.4, int rbt=1){
	//First fit for ensuring the rangement of histgram;
	TH1* h = (TH1*)ht->Clone();
	h->Rebin(rbt);
	TF1 *fit = new TF1("fit","gaus",*tRL,*tRR);        
	h->GetXaxis()->SetRangeUser(*tRL,*tRR);
	h->Fit(fit);
	float mean = fit->GetParameter(1);
	float sigma = TMath::Abs(fit->GetParameter(2));
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

	h->Fit(fit2,"","",mean-10*sigma,mean+10*sigma);
	
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

	//if(*tRL<mean-5*sigma)
	//{
	*tRL = mean-10*sigma;
	*tRR = mean+10*sigma;
	//}

	cout<<h->GetName()<<"\t"<<*tRL<<"\t"<<*tRR<<endl;
	h->GetXaxis()->SetRangeUser(*tRL,*tRR);
}

TF1* profilefit(TH2* Rt,double rbU,double rbt,double tRL,double tRR,double URL,double URR,char* name){
	char buff[1024];


	TCanvas *c5 = new TCanvas("c5","c5",1600,600);
	c5->Clear();
	c5->Divide(2,1);
	c5->cd(1);
	TH2* Qt = (TH2*)Rt->Clone();
	//TH2* Qt = (TH2*) Rt->Clone("tmp");

	Qt->Draw("colz");

	Qt->RebinX(rbU);
	Qt->RebinY(rbt);
	Qt->GetYaxis()->SetRangeUser(tRL,tRR);
	Qt->GetXaxis()->SetRangeUser(URL,URR);
	Qt->ProfileX();
	c5->cd(2);


	TH1* Qpfx = Qt->ProfileX();
	//Qpfx->Reset();


	//Qpfx=(TH1*)gDirectory->Get("Qt_pfx");
	Qpfx->Draw();
	Qpfx->GetYaxis()->SetRangeUser(tRL,tRR);
	Qpfx->GetXaxis()->SetRangeUser(URL,URR);

	TF1* fitQt = new TF1("fitQt","[0]+[1]/TMath::Sqrt(abs(x))+[2]/abs(x)+[3]/abs(x)/TMath::Sqrt(abs(x))+[4]/abs(x)/abs(x)",0,2e3);

	Qpfx->Fit(fitQt,"R");

	sprintf(buff,"%s.png",name);
	c5->SaveAs(buff);
	Qpfx->Reset();
	//delete Qt;
	//delete Qpfx;
	//delete c5;
	return fitQt;


}

TF1* CH1ACorrection(TTree* t1,EVENT* A,EVENT* B,EVENT* MCP,int rbU,int rbt,char* name,double fac,int iter){
	//Correction QT
	char buff[1024]; 


	TCanvas* c1 = new TCanvas("c1","c1",800,600); 	


	sprintf(buff,"%s_ht",name);
	TH1D *ht = new TH1D(buff,"Time Resolution;T (ps);Counts",2e3,-10, 10);




	sprintf(buff,"%s_hAT",name);
	TH2D *hAT = new TH2D(buff,"",200,0,1.6e3,2e3,-10,10);


	TF1 *fitAT=new TF1("fitAT","0",-1e4,1e4);





	double tL=2;
	double tR=6;
	double QL=0;
	double QR=1.6e3;

	double TAcor=0;
	double Q=0;
	double T[1000000]={0};
	for(int s =0;s<iter;s++){



		int N = t1->GetEntries();
		cout<< "Total entries is :"<<N<<endl;
		for(int i=0; i<N; i++)
		{
			t1->GetEntry(i);

			Q=(*B).Q;
			TAcor=fitAT->Eval(Q);
			//if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3){
			if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3&&(*A).rise<4&&(*B).rise<4){


				if((*A).Q<600){
				//if(1){
					if (s==0) 
						//T[i]=(*A).p20-(*MCP).p20-TAcor;
						T[i]=(*B).p20-(*MCP).p20-TAcor;
					else 
						T[i]=T[i]-TAcor;


					ht->Fill(T[i]);		
//				cout<<"Q="<<Q<<" T="<<T[i]<<endl;
					hAT->Fill(Q,T[i]);
				}
			}
		}

//cout<<"<====process check====>"<<endl;

		c1->cd();
		c1->Clear();
		ht->Draw();
		//return;
		//tRL=-5;
		//tRR=5;
		tL=ht->GetMean()-3*ht->GetRMS();
		tR=ht->GetMean()+3*ht->GetRMS();
		//cout<<"iter="<<s<<", MEAN="<<ht->GetMean()<<"\t RMS="<<ht->GetRMS()<<endl;
		//cout<<"iter="<<s<<", tL="<<tL<<"\t tR="<<tR<<endl;
		//if (s==1) return;
		twoguasfit(ht,&tL,&tR,fac,rbt);
		sprintf(buff,"%s_TR_cor%d.png",name,s);
		c1->SaveAs(buff);

		sprintf(buff,"%s_At_pfx_cor%d",name,s);

		fitAT=profilefit(hAT,rbU,rbt*2,tL,tR,QL,QR,buff);

		ht->Reset();
		hAT->Reset();
	}


	return fitAT;

}
TF1* CH1BCorrection(TTree* t1,EVENT* A,EVENT* B,EVENT* MCP,int rbU,int rbt,char* name,double fac,int iter){
	//Correction QT
	char buff[1024]; 


	TCanvas* c1 = new TCanvas("c1","c1",800,600); 	


	sprintf(buff,"%s_ht",name);
	TH1D *ht = new TH1D(buff,"Time Resolution;T (ps);Counts",2e3,-10, 10);




	sprintf(buff,"%s_hAT",name);
	TH2D *hAT = new TH2D(buff,"",200,0,1.6e3,2e3,-10,10);


	TF1 *fitAT=new TF1("fitAT","0",-1e4,1e4);





	double tL=2;
	double tR=6;
	double QL=0;
	double QR=1.6e3;

	double TAcor=0;
	double Q=0;
	double T[100000]={0};
	for(int s =0;s<iter;s++){



		int N = t1->GetEntries();
		cout<< "Total entries is :"<<N<<endl;
		for(int i=0; i<N; i++)
		{
			t1->GetEntry(i);

			Q=(*B).Q;
			//cout<<Q<<endl;

			TAcor=fitAT->Eval(Q);
			//if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3){
			if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3&&(*A).rise<4&&(*B).rise<4){


				if((*B).Q<600&&(*B).rise<4){
					if (s==0) 
						//T[i]=(*A).p20-(*MCP).p20-TAcor;
						T[i]=(*B).p20-(*MCP).p20-TAcor;
					else 
						T[i]=T[i]-TAcor;


					ht->Fill(T[i]);		
					hAT->Fill(Q,T[i]);
				}
			}
		}


		c1->cd();
		c1->Clear();
		ht->Draw();
		//return;
		//tRL=-5;
		//tRR=5;
		tL=ht->GetMean()-3*ht->GetRMS();
		tR=ht->GetMean()+3*ht->GetRMS();
		//cout<<"iter="<<s<<", MEAN="<<ht->GetMean()<<"\t RMS="<<ht->GetRMS()<<endl;
		//cout<<"iter="<<s<<", tL="<<tL<<"\t tR="<<tR<<endl;
		//if (s==1) return;
		twoguasfit(ht,&tL,&tR,fac,rbt);
		sprintf(buff,"%s_TR_cor%d.png",name,s);
		c1->SaveAs(buff);

		sprintf(buff,"%s_At_pfx_cor%d",name,s);

		fitAT=profilefit(hAT,rbU,rbt*2,tL,tR,QL,QR,buff);

		ht->Reset();
		hAT->Reset();
	}


	return fitAT;

}


TF1* CH3Correction(TTree* t1,EVENT* A,EVENT* B,EVENT* MCP,int rbU,int rbt,char* name,double fac,int iter){
	//Correction QT
	char buff[1024]; 




	TCanvas* c4 = new TCanvas("c4","c4",800,600); 	
	sprintf(buff,"%s_ht",name);
	TH1D *ht = new TH1D(buff,"Time Resolution;T (ps);Counts",2e3,-10, 10);




	sprintf(buff,"%s_hAT",name);
	TH2D *hAT = new TH2D(buff,"",200,0,1.6e3,2e3,-10,10);


	TF1 *fitAT=new TF1("fitAT","0",-1e4,1e4);





	double tL=2;
	double tR=6;
	double QL=0;
	double QR=1.6e3;

	double TAcor=0;
	double Q1=0,Q2=0;
	double T[1000000]={0};
	const int chN=2;
	EVENT** ch[chN];

	ch[0]=&A;
	ch[1]=&B;
	//ch[2]=&MCP;
	for(int s =0;s<iter;s++){



		int N = t1->GetEntries();
		cout<< "Total entries is :"<<N<<endl;

		//correction of A1
		for (int h= 0;h<chN;h++){
		tL=-2;
		tR=2;
		
		for(int i=0; i<N; i++)
		{
			t1->GetEntry(i);

			//Q1=(*A).Q;
			//Q2=(*B).Q;
			if (h-1<0)
				TAcor=fitAT->Eval((**ch[chN-1]).Q);
			else
				TAcor=fitAT->Eval((**ch[h-1]).Q);
			
			//TAcor=fitAT->Eval(Q2);
			if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3&&(*A).rise<4&&(*B).rise<4){
				//if(1){
				if((*A).Q<600&&(*B).Q<600){
				if (s==0) 
					T[i]=((*A).p20+(*B).p20)/2-(*MCP).p20-TAcor;
				else 
					T[i]=T[i]-TAcor;


				ht->Fill(T[i]);		
				hAT->Fill((**ch[h]).Q,T[i]);
			}
			}

		}


		c4->cd();
		c4->Clear();
		ht->Draw();
		//return;
		//tRL=-5;
		//tRR=5;

		twoguasfit(ht,&tL,&tR,fac,rbt);
		sprintf(buff,"%s_TR_ch%d_cor%d.png",name,h,s);
		c4->SaveAs(buff);

		sprintf(buff,"%s_At_ch%d_pfx_cor%d",name,h,s);

		fitAT=profilefit(hAT,rbU,rbt*2,tL,tR,QL,QR,buff);

		ht->Reset();
		hAT->Reset();
		}

	}


	return fitAT;

}

TF1* CH2Correction(TTree* t1,EVENT* A,EVENT* B,EVENT* MCP,int rbU,int rbt,char* name,double fac,int iter){
	//Correction QT
	char buff[1024]; 




	TCanvas* c4 = new TCanvas("c4","c4",800,600); 	
	sprintf(buff,"%s_ht",name);
	TH1D *ht = new TH1D(buff,"Time Resolution;T (ps);Counts",2e3,-10, 10);




	sprintf(buff,"%s_hAT",name);
	TH2D *hAT = new TH2D(buff,"",200,0,1.6e3,2e3,-10,10);


	TF1 *fitAT=new TF1("fitAT","0",-1e4,1e4);





	double tL=-2;
	double tR=2;
	double QL=0;
	double QR=1.6e3;

	double TAcor=0;
	double Q1=0,Q2=0;
	double T[1000000]={0};
	const int chN=2;
	EVENT** ch[chN];

	ch[0]=&A;
	ch[1]=&B;
	//ch[2]=&MCP;
	for(int s =0;s<iter;s++){



		int N = t1->GetEntries();
		cout<< "Total entries is :"<<N<<endl;

		//correction of A1
		for (int h= 0;h<chN;h++){
		tL=-2;
		tR=2;
		
		for(int i=0; i<N; i++)
		{
			t1->GetEntry(i);

			//Q1=(*A).Q;
			//Q2=(*B).Q;
			if (h-1<0)
				TAcor=fitAT->Eval((**ch[chN-1]).Q);
			else
				TAcor=fitAT->Eval((**ch[h-1]).Q);
			
			//TAcor=fitAT->Eval(Q2);
			//if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*A).rise<4&&(*B).rise<4){
			if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*A).rise<4&&(*B).rise<4&&(*MCP).Q>0.3){
				//if(1){
				//if(1){
				if((*A).A>0.3&&(*B).A>0.3){
				if (s==0&&h==0) 
					//T[i]=(*A).p20-(*B).p20-TAcor;
					T[i]=(*A).p20-(*B).p20;
				else 
					T[i]=T[i]-TAcor;


				ht->Fill(T[i]);		
				hAT->Fill((**ch[h]).Q,T[i]);
			}
			}

		}


		c4->cd();
		c4->Clear();
		ht->Draw();
		//return;
		//tRL=-5;
		//tRR=5;

		twoguasfit(ht,&tL,&tR,fac,rbt);
		sprintf(buff,"%s_TR_ch%d_cor%d.png",name,h,s);
		c4->SaveAs(buff);

		sprintf(buff,"%s_At_ch%d_pfx_cor%d",name,h,s);

		fitAT=profilefit(hAT,rbU,rbt*2,tL,tR,QL,QR,buff);

		ht->Reset();
		hAT->Reset();
		}

	}


	return fitAT;

}


