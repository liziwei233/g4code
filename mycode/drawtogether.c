//
//* Environment: 
//**     class LZWfunc
//**     class DrawMyfunc
//* Data: root file.
//* Function: 
//*		draw histgrams come from several root file together.
//* Date: 2019.3.19
//*
//***********************************
//

#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>
#include <TRandom.h>


using namespace std;

void drawtogether()
{
	/*void MygStyle();
    void SetMyPad(TPad *pad,float left, float right, float top, float bottom);
    void DrawMyHist1(TH1F *datahist, char *xtitle, char *ytitle, int LColor=1, int LWidth=3, int TitleColor=1,int LStyle=1);
    void DrawMyHist2(TH2F *datahist, char *xtitle, char *ytitle, int LColor=1, int LWidth=3, int TitleColor=1,int LStyle=1);
    TLegend* DrawMyLeg(Double_t xlow=0.2, Double_t ylow=0.2, Double_t xup=0.5, Double_t yup=0.5, Int_t textFont=62, Double_t textSize=0.05);
    TLatex* DrawMyLatex(char* text, Double_t x=0.65, Double_t y=0.5, Int_t textFont=62, Double_t textSize=0.05, Int_t colorIndex=2);
	*/
	TLegend *DrawMyLeg(Double_t xlow = 0.2, Double_t ylow = 0.2, Double_t xup = 0.5, Double_t yup = 0.5, Int_t textFont = 62, Double_t textSize = 0.05);




	gStyle->SetOptFit(1111);

	char buff[1024];
	char name[1024];


	//ofstream output("timeres.dat");

	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	c1->Clear();

	//TCanvas *c2 = new TCanvas("c2","c2",800,600);
	//TCanvas *c3 = new TCanvas("c3","c3",800,600);

	//TCanvas *c5 = new TCanvas("c4","c5",800,600);

	//sprintf(name,"%dpe_57v",npe[s]);

	int radius[]={0,3,5,7,9};

	float UL = 0;
	float UR = 1000;
	int n = sizeof(radius)/sizeof(radius[0]);
	TH2D *h[n];

	Color_t clr[8] = {kRed + 3, kOrange + 10, kOrange, kGreen + 2, kAzure + 1, kBlue - 3, kViolet + 1, kViolet};
	double tL=1.4;
	double tR=1.6;
	int bint=200;
	double zL=-26;
	double zR=-21;
	int binz=200;
	double entry=1;
	for (int i = 0; i < n; i++)
	{

		//c1->Clear();

		sprintf(name, "bare_r%dmm", radius[i]);
		sprintf(buff, "%s.root", name);

		TFile *f1 = new TFile(buff, "READ");
		TTree *t1 = (TTree *)f1->Get("Run");
		
		c1->cd();
		TH2D *hzt= new TH2D("hzt",";Arrival Time (ns);Z axis (mm)",bint,tL,tR,binz,zL,zR);

		t1->Draw("PmtL.bpz:PmtL.t>>hzt","","colz");
		sprintf(buff,"%s_zt.png",name);
		if(i==0) {
			entry=hzt->GetEntries();
			cout<< "the entry is "<<entry<<endl;
		}
		
		h[i] = ht->Scale(1/entry);
		
	}
	//return;
	TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
	c2->Clear();
	c2->cd();

	TLegend *leg;
	leg = DrawMyLeg(0.6, 0.2, 0.8, 0.7);

	for (int s = 0; s < n; s++)
	{

		h[s]->SetLineWidth(2);
		h[s]->SetLineColor(1);
		//h[s]->SetFillColor(1+s);
		//h[s]->SetFillStyle(3020);
		h[s]->SetFillStyle(3017);
		h[s]->SetFillColorAlpha(clr[t], 1);
		h[s]->SetLineStyle(2);

		//h[i]->Draw();
		h[s]->Draw("same");
		sprintf(buff, "r=%dmm", radius[s]);
		leg->AddEntry(h[s], buff, "lpf");
	}
	leg->Draw();



	
}

TLegend *DrawMyLeg(Double_t xlow = 0.2, Double_t ylow = 0.2, Double_t xup = 0.5, Double_t yup = 0.5, Int_t textFont = 62, Double_t textSize = 0.05)
{
	TLegend *leg = new TLegend(xlow, ylow, xup, yup);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetFillColor(10);
	leg->SetTextFont(textFont);
	leg->SetTextSize(textSize);
	//leg->Draw("same");
	return leg;
};

TF1 *gausfit(TH1 *hU, int rbU, float *U_RL, float *U_RR)
{

	hU->Draw();
	hU->Rebin(rbU);
	TF1 *fitU = new TF1("fitU", "gaus", *U_RL, *U_RR);
	fitU->SetParameter(1, hU->GetMean());
	hU->Fit(fitU, "R");
	*U_RL = fitU->GetParameter(1) - 5 * fitU->GetParameter(2);
	*U_RR = fitU->GetParameter(1) + 5 * fitU->GetParameter(2);

	hU->GetXaxis()->SetRangeUser(*U_RL, *U_RR);
	hU->Fit(fitU);
	return fitU;
}

void twoguasfit(TH1 *ht, double *tRL, double *tRR, double fac = 0.4, int rbt = 1)
{
	//First fit for ensuring the rangement of histgram;
	TH1 *h = (TH1 *)ht->Clone();
	h->Rebin(rbt);
	TF1 *fit = new TF1("fit", "gaus", *tRL, *tRR);
	h->GetXaxis()->SetRangeUser(*tRL, *tRR);
	h->Fit(fit);
	float mean = fit->GetParameter(1);
	float sigma = TMath::Abs(fit->GetParameter(2));
	if (*tRL < mean - 5 * sigma || sigma > 1)
	{
		*tRL = mean - 5 * sigma;
		*tRR = mean + 5 * sigma;
	}
	cout << h->GetName() << "\t" << *tRL << "\t" << *tRR << endl;

	h->GetXaxis()->SetRangeUser(*tRL, *tRR);

	TF1 *fit2 = new TF1("fit2", "gaus(0)+gaus(3)", *tRL, *tRR);
	fit2->SetParNames("C_{TR}", "#mu_{TR}", "#sigma_{TR}", "C_{nz}", "#mu_{nz}", "#sigma_{nz}");
	fit2->SetParameter(1, mean);
	fit2->SetParameter(2, sigma);
	fit2->SetParLimits(3, 0, fit->GetParameter(0) * fac);
	fit2->SetParameter(4, mean);
	fit2->SetParameter(5, 3 * sigma);

	h->Fit(fit2, "", "", mean - 5 * sigma, mean + 5 * sigma);
	TF1 *fit_tr = new TF1("fit_tr", "gaus", *tRL, *tRR);
	fit_tr->SetParameter(0, fit2->GetParameter(0));
	fit_tr->SetParameter(1, fit2->GetParameter(1));
	fit_tr->SetParameter(2, fit2->GetParameter(2));
	TF1 *fit_nz = new TF1("fit_nz", "gaus", *tRL, *tRR);
	fit_nz->SetParameter(0, fit2->GetParameter(3));
	fit_nz->SetParameter(1, fit2->GetParameter(4));
	fit_nz->SetParameter(2, fit2->GetParameter(5));
	fit_tr->SetLineColor(3);
	fit_tr->SetLineStyle(7);
	fit_nz->SetLineColor(3);
	fit_nz->SetLineStyle(7);
	fit_tr->Draw("same");
	fit_nz->Draw("same");

	if (fit2->GetParameter(3) >= fit2->GetParameter(0) && TMath::Abs(fit2->GetParameter(5)) >= 0.05)
	{
		mean = fit2->GetParameter(4);
		sigma = TMath::Abs(fit2->GetParameter(5));
	}
	else if (TMath::Abs(fit2->GetParameter(2)) >= 0.05)
	{
		mean = fit2->GetParameter(1);
		sigma = TMath::Abs(fit2->GetParameter(2));
	}

	//if(*tRL<mean-5*sigma)
	//{
	*tRL = mean - 5 * sigma;
	*tRR = mean + 5 * sigma;
	//}

	cout << h->GetName() << "\t" << *tRL << "\t" << *tRR << endl;
	h->GetXaxis()->SetRangeUser(*tRL, *tRR);
}

TF1 *profilefit(TH2 *Rt, double rbU, double rbt, double tRL, double tRR, double URL, double URR, char *name)
{
	char buff[1024];

	TCanvas *c5 = new TCanvas("c5", "c5", 1600, 600);
	c5->Clear();
	c5->Divide(2, 1);
	c5->cd(1);
	TH2 *Qt = (TH2 *)Rt->Clone();
	//TH2* Qt = (TH2*) Rt->Clone("tmp");

	Qt->Draw("colz");

	Qt->RebinX(rbU);
	Qt->RebinY(rbt);
	Qt->GetYaxis()->SetRangeUser(tRL, tRR);
	Qt->GetXaxis()->SetRangeUser(URL, URR);
	Qt->ProfileX();
	c5->cd(2);

	TH1 *Qpfx = Qt->ProfileX();
	//Qpfx->Reset();

	//Qpfx=(TH1*)gDirectory->Get("Qt_pfx");
	Qpfx->Draw();
	Qpfx->GetYaxis()->SetRangeUser(tRL, tRR);
	Qpfx->GetXaxis()->SetRangeUser(URL, URR);

	TF1 *fitQt = new TF1("fitQt", "[0]+[1]/TMath::Sqrt(abs(x))+[2]/abs(x)+[3]/abs(x)/TMath::Sqrt(abs(x))+[4]/abs(x)/abs(x)", 0, 2e3);

	Qpfx->Fit(fitQt, "R");

	sprintf(buff, "%s.png", name);
	c5->SaveAs(buff);
	Qpfx->Reset();
	//delete Qt;
	//delete Qpfx;
	//delete c5;
	return fitQt;
}
TF1 *Correction(TTree *t1, double *A1Q, double *A1, double *A1rise, double rbU, double rbt, double rbR, double tRL, double tRR, double URL, double URR, double RRL, double RRR, char *name)
{
	//Correction QT
	char buff[1024];
	TH1D *ht_before = new TH1D("htb", "Time Resolution;T (ps);Counts", 3e3, -20, 230);

	TH1D *ht_afterRt = new TH1D("htaRt", "Time Resolution;T (ps);Counts", 200, -20, 20);

	TH1D *hR = new TH1D("hR", ";Risetime (ns);Counts", 2e3, 0, 5);
	TH1D *hU = new TH1D("hU", ";Charge (pC);Counts", 200, URL, URR);

	TH2D *Rt = new TH2D("Rt", "", 2e3, 0, 5, 3e3, -20, 230);

	hR->GetXaxis()->SetRangeUser(RRL, RRR);
	Rt->GetYaxis()->SetRangeUser(RRL, RRR);
	//hU->GetXaxis->SetRangeUser(*URL,*URR);

	int N = 0;
	double Q = 0, T = 0, R = 0;

	double TRcor = 0;

	N = t1->GetEntries();
	cout << "Total entries is :" << N << endl;
	for (int i = 0; i < N; i++)
	{

		t1->GetEntry(i);
		Q = *A1Q;
		T = *A1;
		R = *A1rise;
		//cout<<i<<"\t"<<*A1Q<<"\t"<<*A1<<"\t"<<*A1rise<<endl;

		//cout<<Q<<"\t"<<T<<"\t"<<R<<endl;
		//if(A1Q>URL&&A2Q>URL&&A1rise>RRL&&A1rise<RRR&&A2rise>RRL&&A2rise<RRR)	{
		ht_before->Fill(T);
		hR->Fill(R);
		hU->Fill(Q);
		//Qt->Fill(Q,T);
		Rt->Fill(R, T);
		//cout<<"Risetime:"<<R<<";Charge:"<<Q<<";20%Rise:"<<T<<endl;

		//}
	}

	//fit risetime hist to get the rangement paraments;
	TCanvas *cR = new TCanvas("cR", "cR", 800, 600);
	cR->cd();
	//hR->Draw();

	twoguasfit(hR, &RRL, &RRR);
	//return;
	sprintf(buff, "%s_risetime.png", name);
	cR->SaveAs(buff);
	cR->Close();

	//fit charge hist to get the rangement paraments;
	TCanvas *cU = new TCanvas("cU", "cU", 800, 600);
	cU->cd();
	//hU->Draw();
	//return;
	twoguasfit(hU, &URL, &URR);
	//return;
	sprintf(buff, "%s_Charge.png", name);
	cU->SaveAs(buff);
	cU->Close();

	TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
	//c4->Clear();
	c4->cd();
	//ht_before->Draw();
	//return;
	twoguasfit(ht_before, &tRL, &tRR);
	sprintf(buff, "%s_tr_before.png", name);
	c4->SaveAs(buff);

	//make Rt profileX and fit
	//--------------------------------
	sprintf(buff, "%s_Rt_pfx_before", name);
	TF1 *fitRt = profilefit(Rt, rbR, rbt, tRL, tRR, RRL, RRR, buff);

	//Rt->ProfileX()->Draw();
	//return;
	Rt->Reset();

	//add a branch after correction
	char buff2[1024];
	sprintf(buff, "%s_Tcor", name);
	sprintf(buff2, "%s_Tcor/D", name);
	TBranch *bcor = t1->Branch(buff, &T, buff2);

	for (int i = 0; i < N; i++)
	{
		t1->GetEntry(i);

		R = *A1rise;
		TRcor = fitRt->Eval(R);

		T = *A1 - TRcor;
		bcor->Fill(); //cout<<"T0-T0_cor="<<T<<"-"<<T_cor<<"="<<T-T_cor<<endl;

		//if(A1Q>URL&&A2Q>URL&&A1rise>RRL&&A1rise<RRR&&A2rise>RRL&&A2rise<RRR)	{
		ht_afterRt->Fill(T);
		Rt->Fill(R, T);

		//}
		//cout<<"U:T0-T0_cor="<<UR<<"\t"<<T0-T0_cor<<endl;
	}
	t1->Write();

	c4->cd();
	c4->Clear();
	ht_afterRt->Draw();
	//return;
	tRL = -5;
	tRR = 5;
	twoguasfit(ht_afterRt, &tRL, &tRR);
	sprintf(buff, "%s_tr_afterRt.png", name);
	c4->SaveAs(buff);

	sprintf(buff, "%s_Rt_pfx_after", name);

	profilefit(Rt, rbR, rbt, tRL, tRR, RRL, RRR, buff);

	/*
		TF1* fitRt2 = profilefit(Rt,rbR,rbt,*tRL,*tRR,*RRL,*RRR,buff);
		
		//Rt->ProfileX()->Draw();
		//return;
		Rt->Reset();
		ht_afterRt->Reset();
		double TRcor2;
		for(int i=0; i<N; i++)
		{
		t1->GetEntry(i);

		R=*A1rise;
		TRcor2=fitRt2->Eval(R);
		T=*A1-TRcor-TRcor2;		//cout<<"T0-T0_cor="<<T<<"-"<<T_cor<<"="<<T-T_cor<<endl;
		
		//if(A1Q>URL&&A2Q>URL&&A1rise>RRL&&A1rise<RRR&&A2rise>RRL&&A2rise<RRR)	{
		ht_afterRt->Fill(T);		
		Rt->Fill(R,T);
		
			//}
//cout<<"U:T0-T0_cor="<<UR<<"\t"<<T0-T0_cor<<endl;
		}
		
		c4->cd();
		c4->Clear();
		ht_afterRt->Draw();
		//return;
		twoguasfit(ht_afterRt,tRL,tRR);
		sprintf(buff,"%s_tr_afterRt2.png",name);
		c4->SaveAs(buff);
		
		sprintf(buff,"%s_Rt_pfx_after2",name);
		
		profilefit(Rt,rbR,rbt,*tRL,*tRR,*RRL,*RRR,buff);
		*/

	return fitRt;
}

TF1 *CH1Correction(TTree *t1, double *A1, double *A1y, double *A1Q, double *T0, double *T0y, double *T0Q, int rbU, int rbt, char *name, int iter)
{
	//Correction QT
	char buff[1024];

	sprintf(buff, "%s_ht", name);
	TH1D *ht = new TH1D(buff, "Time Resolution;T (ps);Counts", 3e3, -10, 20);

	sprintf(buff, "%s_hQ", name);
	TH1D *hQ = new TH1D(buff, ";charge (pC);Counts", 3e3, 0, 1.6e3);

	sprintf(buff, "%s_hAT", name);
	TH2D *hAT = new TH2D(buff, "", 3e3, 0, 1.6e3, 3e3, -10, 20);

	TF1 *fitAT = new TF1("fitAT", "0", -1e4, 1e4);

	TCanvas *c4 = new TCanvas("c4", "", 800, 600);

	double tL = 10;
	double tR = 16;
	double QL = 0;
	double QR = 1.6e3;

	double TAcor = 0;
	double Q = 0;
	double T[100000] = {0};
	for (int s = 0; s < iter; s++)
	{

		int N = t1->GetEntries();
		cout << "Total entries is :" << N << endl;
		for (int i = 0; i < N; i++)
		{
			t1->GetEntry(i);

			Q = *A1Q;

			TAcor = fitAT->Eval(Q);
			if (*A1y > 1 &&
				*A1y < 3)
			{
				if (s == 0)
					T[i] = *A1 - *T0 - TAcor;
				else
					T[i] = T[i] - TAcor;

				hQ->Fill(Q);
				ht->Fill(T[i]);
				hAT->Fill(Q, T[i]);
			}
		}

		c4->cd();
		c4->Clear();
		ht->Draw();
		//return;
		//tRL=-5;
		//tRR=5;
		tL = ht->GetMean() - 3 * ht->GetRMS();
		tR = ht->GetMean() + 3 * ht->GetRMS();
		//cout<<"iter="<<s<<", MEAN="<<ht->GetMean()<<"\t RMS="<<ht->GetRMS()<<endl;
		//cout<<"iter="<<s<<", tL="<<tL<<"\t tR="<<tR<<endl;
		//if (s==1) return;
		twoguasfit(ht, &tL, &tR, 0.2, 8);
		sprintf(buff, "%s_TR_cor%d.png", name, s);
		c4->SaveAs(buff);

		c4->Clear();
		hQ->Draw();
		//return;
		//tRL=-5;
		//tRR=5;
		QL = hQ->GetMean() - 3 * hQ->GetRMS();
		QR = hQ->GetMean() + 3 * hQ->GetRMS();
		//cout<<"iter="<<s<<", MEAN="<<ht->GetMean()<<"\t RMS="<<ht->GetRMS()<<endl;
		//cout<<"iter="<<s<<", tL="<<tL<<"\t tR="<<tR<<endl;
		//if (s==1) return;
		twoguasfit(hQ, &QL, &QR, 0.2, 1);
		sprintf(buff, "%s_charge_cor%d.png", name, s);
		c4->SaveAs(buff);

		sprintf(buff, "%s_At_pfx_cor%d", name, s);

		fitAT = profilefit(hAT, rbU, rbt, tL, tR, QL, QR, buff);

		ht->Reset();
		hAT->Reset();
	}

	return fitAT;
}
TF1 *CH2Correction(TTree *t1, double *A1, double *A1y, double *A1Q, double *A2, double *A2y, double *A2Q, double *T0, double *T0y, double *T0Q, int rbU, int rbt, char *name, double thrd, int iter)
{
	//Correction QT
	char buff[1024];

	TCanvas *c4 = new TCanvas("c4", "", 800, 600);

	sprintf(buff, "%s_ht", name);
	TH1D *ht = new TH1D(buff, "Time Resolution;T (ps);Counts", 2e3, -10, 10);

	sprintf(buff, "%s_hAT", name);
	TH2D *hAT = new TH2D(buff, "", 200, 0, 1.6e3, 2e3, -10, 10);

	TF1 *fitAT = new TF1("fitAT", "0", -1e4, 1e4);

	double tL = 2;
	double tR = 6;
	double QL = 0;
	double QR = 1.6e3;

	double TAcor = 0;
	double Q1 = 0, Q2 = 0;
	double T[100000] = {0};
	for (int s = 0; s < iter; s++)
	{

		int N = t1->GetEntries();
		cout << "Total entries is :" << N << endl;

		//correction of A1
		//TCut cut="";
		for (int i = 0; i < N; i++)
		{
			t1->GetEntry(i);

			Q1 = *A1Q;
			Q2 = *A2Q;
			TAcor = fitAT->Eval(Q2);
			if (*A1y > 0.9 &&
				*A1y < 1.05 &&
				*A2Q > 400 &&
				*A1<*A2 && * T0Q> 0.3)
			{
				if (s == 0)
					T[i] = (*A1 + *A2) / 2 - *T0 - TAcor;
				else
					T[i] = T[i] - TAcor;

				ht->Fill(T[i]);
				hAT->Fill(Q1, T[i]);
			}
		}

		c4->cd();
		c4->Clear();
		ht->Draw();
		//return;
		//tRL=-5;
		//tRR=5;

		twoguasfit(ht, &tL, &tR, 0.2, 2);
		sprintf(buff, "%s_TR_ch1_cor%d_thrd%g.png", name, s, thrd);
		c4->SaveAs(buff);

		sprintf(buff, "%s_At_ch1_pfx_cor%d_thrd%g", name, s, thrd);

		fitAT = profilefit(hAT, rbU, rbt, tL, tR, QL, QR, buff);

		ht->Reset();
		hAT->Reset();

		//correction of A2
		tL = -2;
		tR = 2;
		for (int i = 0; i < N; i++)
		{
			t1->GetEntry(i);

			Q1 = *A1Q;
			Q2 = *A2Q;
			TAcor = fitAT->Eval(Q1);
			if (*A1y > 0.9 &&
				*A1y < 1.05 &&
				*A2Q > 400 &&
				*A1<*A2 && * T0Q> 0.3)
			{
				T[i] = T[i] - TAcor;

				ht->Fill(T[i]);
				hAT->Fill(Q2, T[i]);
			}
		}

		//c4->cd();
		c4->Clear();
		ht->Draw();
		//return;
		//tRL=-5;
		//tRR=5;

		twoguasfit(ht, &tL, &tR, 0.2, 2);
		sprintf(buff, "%s_TR_ch2_cor%d_thrd%g.png", name, s, thrd);
		c4->SaveAs(buff);

		sprintf(buff, "%s_At_ch2_pfx_cor%d_thrd%g", name, s, thrd);

		fitAT = profilefit(hAT, rbU, rbt, tL, tR, QL, QR, buff);

		ht->Reset();
		hAT->Reset();
	}

	return fitAT;
}

TF1 *CH2Correction_noT0(TTree *t1, double *A1, double *A1y, double *A1Q, double *A2, double *A2y, double *A2Q, int rbU, int rbt, char *name, double thrd, int iter)
{
	//Correction QT
	char buff[1024];

	sprintf(buff, "%s_ht", name);
	TH1D *ht = new TH1D(buff, "Time Resolution;T (ns);Counts", 2e3, -10, 10);

	sprintf(buff, "%s_hAT", name);
	TH2D *hAT = new TH2D(buff, "", 200, 0, 1.6e3, 2e3, -10, 10);

	TF1 *fitAT = new TF1("fitAT", "0", -1e4, 1e4);

	double tL = -2;
	double tR = 2;
	double QL = 0;
	double QR = 1.6e3;

	double TAcor = 0;
	double Q1 = 0, Q2 = 0;
	double T[100000] = {0};
	for (int s = 0; s < iter; s++)
	{

		int N = t1->GetEntries();
		cout << "Total entries is :" << N << endl;

		//correction of A1

		for (int i = 0; i < N; i++)
		{
			t1->GetEntry(i);

			Q1 = *A1Q;
			Q2 = *A2Q;
			//cout<<Q1<<"\t"<<Q2<<endl;
			TAcor = fitAT->Eval(Q2);
			if (Q2 < thrd)
			{
				if (s == 0)
					T[i] = *A1 - *A2;
				else
					T[i] = T[i] - TAcor;

				ht->Fill(T[i]);
				hAT->Fill(Q1, T[i]);
			}
		}

		c4->cd();
		c4->Clear();
		ht->Draw();
		//return;
		//tRL=-5;
		//tRR=5;

		twoguasfit(ht, &tL, &tR, 0.2, 2);
		sprintf(buff, "%s_TR_ch1_cor%d_thrd%g.png", name, s, thrd);
		c4->SaveAs(buff);

		sprintf(buff, "%s_At_ch1_pfx_cor%d_thrd%g", name, s, thrd);

		fitAT = profilefit(hAT, rbU, rbt, tL, tR, QL, QR, buff);

		ht->Reset();
		hAT->Reset();

		//correction of A2
		tL = -2;
		tR = 2;
		for (int i = 0; i < N; i++)
		{
			t1->GetEntry(i);

			Q1 = *A1Q;
			Q2 = *A2Q;
			TAcor = fitAT->Eval(Q1);
			if (Q2 < thrd)
			{
				T[i] = T[i] - TAcor;

				ht->Fill(T[i]);
				hAT->Fill(Q2, T[i]);
			}
		}

		//c4->cd();
		c4->Clear();
		ht->Draw();
		//return;
		//tRL=-5;
		//tRR=5;

		twoguasfit(ht, &tL, &tR, 0.2, 2);
		sprintf(buff, "%s_TR_ch2_cor%d_thrd%g.png", name, s, thrd);
		c4->SaveAs(buff);

		sprintf(buff, "%s_At_ch2_pfx_cor%d_thrd%g", name, s, thrd);

		fitAT = profilefit(hAT, rbU, rbt, tL, tR, QL, QR, buff);

		ht->Reset();
		hAT->Reset();
	}

	return fitAT;
}

//void draw_TR(TTree *t1,double** A1,double** A2,double** B1,double** B2,double tRL,double tRR,double URL,double URR,double RRL,double RRR,char* name){

void draw_TR(TTree *t1, double **A1, double **A2, double **B1, double **B2, double *A1_cTQR_par, double *A2_cTQR_par, char *name)
{
	char buff[1024];
	TH1F *hA = new TH1F("hA", "Time Resolution;T (ps);Counts", 300, -10, 5);

	int N = 0;
	double var[4][4] = {0};
	double T = 0;

	double A1_cTQR[6], A2_cTQR[6];

	for (int i = 0; i < 6; i++)
	{
		A1_cTQR[i] = A1_cTQR_par[i];
		A2_cTQR[i] = A2_cTQR_par[i];
	}

	N = t1->GetEntries();

	for (int i = 0; i < N; i++)
	{

		t1->GetEntry(i);
		// 0, 1, 2 ,T, Q, R
		var[0][0] = *A1[0];
		var[0][1] = *A1[1];
		var[0][2] = *A1[2];
		var[0][3] = *A1[3];

		var[1][0] = *A2[0];
		var[1][1] = *A2[1];
		var[1][2] = *A2[2];
		var[1][3] = *A2[3];

		var[2][0] = *B1[0];
		var[2][1] = *B1[1];
		var[2][2] = *B1[2];
		var[2][3] = *B1[3];

		var[3][0] = *B2[0];
		var[3][1] = *B2[1];
		var[3][2] = *B2[2];
		var[3][3] = *B2[3];

		//cout<<var[0][0]<<"\t"<<var[0][1]<<"\t"<<var[0][2]<<endl;
		//cout<<var[1][0]<<"\t"<<var[0][1]<<"\t"<<var[0][2]<<"\t"<<endl;
		//cout<<var[0][2]<<"\t"<<var[1][2]<<"\t"<<var[2][2]<<"\t"<<var[3][2]<<"\t"<<endl;

		//(A1+B1)/2-(A2+B2)/2
		//(0+2)/2-(1+3)/2
		T = (var[0][0] + var[2][0]) / 2. - (var[1][0] + var[3][0]) / 2.;
		//cout<<T<<endl;
		//cout<<T<<endl;
		//cout<<i<<"\t"<<*A1Q<<"\t"<<*A1<<"\t"<<*A1rise<<endl;

		//cout<<Q<<"\t"<<T<<"\t"<<R<<endl;
		//cout<<var[0][1]<<"\t"<<mcpURL<<endl;
		//cout<<var[0][2]<<"\t"<<mcpRRL<<"\t"<<mcpRRR<<endl;
		//cout<<var[0][3]<<"\t"<<mcptRL<<"\t"<<mcptRR<<endl;
		if (var[0][1] < 0.5 && var[0][2] < 0.3)
		{
			if (var[0][1] > mcpURL && var[0][2] > mcpRRL && var[0][2] < mcpRRR && var[0][3] > mcptRL && var[0][3] < mcptRR &&
				var[2][1] > mcpURL && var[2][2] > mcpRRL && var[2][2] < mcpRRR && var[2][3] > mcptRL && var[2][3] < mcptRR &&
				var[1][1] > URL && var[1][2] > RRL && var[1][2] < RRR && var[1][3] > tRL && var[1][3] < tRR &&
				var[3][1] > URL && var[3][2] > RRL && var[3][2] < RRR && var[3][3] > tRL && var[3][3] < tRR)
			{
				//if(var[0][1]>URL&&var[0][2]>RRL&&var[0][2]<RRR&&var[0][3]>227&&var[1][1]>URL&&var[1][2]>RRL&&var[1][2]<RRR&&var[1][3]>227&&var[2][1]>URL&&var[2][2]>RRL&&var[2][2]<RRR&&var[2][3]>227&&var[3][1]>URL&&var[3][2]>RRL&&var[3][2]<RRR&&var[3][3]>227)	{

				hA->Fill(T);
				//cout<<T<<endl;
				//Qt->Fill(Q,T);

				//cout<<"Risetime:"<<R<<";Charge:"<<Q<<";20%Rise:"<<T<<endl;
			}
		}

		else
		{
			if (var[0][1] > URL && var[0][2] > RRL && var[0][2] < RRR && var[0][3] > tRL && var[0][3] < tRR && var[2][1] > URL && var[2][2] > RRL && var[2][2] < RRR && var[2][3] > tRL && var[2][3] < tRR && var[1][1] > URL && var[1][2] > RRL && var[1][2] < RRR && var[1][3] > tRL && var[1][3] < tRR && var[3][1] > URL && var[3][2] > RRL && var[3][2] < RRR && var[3][3] > tRL && var[3][3] < tRR)
			{
				hA->Fill(T);
			}
		}
	}
	//fit risetime hist to get the rangement paraments;
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	c1->cd();
	hA->Draw();
	//return;
	tRL = -5;
	tRR = 5;
	twoguasfit(hA, &tRL, &tRR);
	sprintf(buff, "%s.png", name);
	c1->SaveAs(buff);
}