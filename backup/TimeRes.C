#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>

void TimeRes(){
	
	//void MygStyle();
	//MygStyle();
	gStyle->SetOptFit(111);
	
	char name[100];
	char buff[1024];
	float RL=1.25e-9;
	float RR=1.35e-9;
	float U_RL=-100;
	float U_RR=-60;
	float thrd;
	for (int s = 14;s<15;s++)
	{
		thrd = (s+1)*0.2;
		sprintf(name,"gun_23cm_Rate_0.1_thrd_%g",thrd);
	sprintf(buff,"%s.root",name);
	TFile *f = new TFile(buff,"READ");
	TTree *t = (TTree*)f->Get("data");
	double UL,UR,T0L,T0R,T0;
	double T0_cor;
	//t->SetBranchAddress("UL",&UL);
	t->SetBranchAddress("UR",&UR);
	//t->SetBranchAddress("T0L",&T0L);
	t->SetBranchAddress("T0R",&T0R);
	t->SetBranchAddress("T0",&T0);
	/*
	
	*/
	TCanvas *c1 = new TCanvas("c1","",800,600);
	TCanvas *c2 = new TCanvas("c2","",800,600);
	TCanvas *c3 = new TCanvas("c3","",800,600);
	//TCanvas *c4 = new TCanvas("c4","",800,600);
	c1->cd();
	
	
	TH1D *h = new TH1D("h","",3600,-1.e-9,2.6e-9);
	t->Draw("T0>>h");
	h->GetXaxis()->SetRangeUser(RL,RR);
	TF1 *fit1 = new TF1("fit1","gaus",-2e-9,2.6e-9);
	h->Fit(fit1);
	c1->SaveAs("T0.png");
	
	ofstream output("timeres.dat",ios::app);
	output<<fit1->GetParameter(2)<<"\t"<<fit1->GetParError(2)<<endl;
	
	TH1D *hpfx;
	TH2D *ht = new TH2D("ht","",100,U_RL,U_RR,3600,-1e-9,2.6e-9);
	c2->cd();
	t->Draw("T0:UR>>ht","","colz");
	
	ht->RebinX(4);
	ht->RebinY(2);
	c2->SaveAs("T0_U.png");
	ht->ProfileX();
	c3->cd();
	hpfx=(TH1D*)gDirectory->Get("ht_pfx");
	hpfx->GetYaxis()->SetRangeUser(RL,RR);
	hpfx->Draw();
	
	TF1 *fit2 = new TF1("fit2","[0]+[1]/TMath::Sqrt(abs(x))+[2]/abs(x)+[3]/abs(x)/TMath::Sqrt(abs(x))+[4]/abs(x)/abs(x)",U_RL,U_RR);
	hpfx->Fit(fit2,"R");
	c3->SaveAs("profileX.png");
	//return;
	//h->Reset();
	
	
	
	TH2D *ht_cor= new TH2D("ht_cor","",100,U_RL,U_RR,3600,-1e-9,2.6e-9);
	hpfx->GetYaxis()->SetRangeUser(-0.3e-9,0.3e-9);
	hpfx->Reset();
	
	TH1D *h_cor = new TH1D("h_cor","",3600,-1e-9,2.6e-9);
	int n=t->GetEntries();
	for(int i=0; i<n; i++)
	{
		t->GetEntry(i);
		T0_cor = fit2->Eval(UR);
		//cout<<"T0-T0_cor="<<T0<<"-"<<T0_cor<<"="<<T0-T0_cor<<endl;
		int temp=-3;
		h_cor->Fill(T0-T0_cor);
		ht_cor->Fill(UR,T0-T0_cor);
		//cout<<"U:T0-T0_cor="<<UR<<"\t"<<T0-T0_cor<<endl;
	}
	c1->cd();
	c1->Clear();
	ht_cor->Draw();
	//return;
	h_cor->Draw();
	h_cor->GetXaxis()->SetRangeUser(-0.1e-9,0.1e-9);
	fit1->SetRange(-0.3e-9,0.3e-9);
	h_cor->Fit(fit1,"R");
	c1->SaveAs("T0_AfterCorrect.png");
	
	//output("timeres.dat",ios::app);
	output<<fit1->GetParameter(2)<<"\t"<<fit1->GetParError(2)<<endl;
	
	c2->cd();
	c2->Clear();
	ht_cor->Draw("colz");
	
	c2->SaveAs("T0_U_AfterCorrect.png");
	ht_cor->ProfileX();
	c3->cd();
	c3->Clear();
	hpfx=(TH1D*)gDirectory->Get("ht_cor_pfx");
	
	hpfx->Draw();

	
	hpfx->Fit(fit2,"R");
	c3->SaveAs("profileX_AfterCorrect.png");
	}
	
}

/*
void MygStyle(){
		//gStyle->SetOptStat(10);
        gStyle->SetOptFit(1111);
		//gStyle->SetOptTitle(0);
        gStyle->SetStatX(0.96);
        gStyle->SetStatY(0.97);
        gStyle->SetStatH(0.17);
        gStyle->SetStatW(0.20);
        gStyle->SetStatStyle(0);
        gStyle->SetPadGridX(0);
        gStyle->SetPadGridY(0);
		*/
        /*gStyle->SetTitleFont(42,"XYZ");
        gStyle->SetTitleFont(42,"T");
        gStyle->SetLabelFont(42,"XYZ");
        gStyle->SetLabelSize(0.05,"XYZ");
        gStyle->SetTitleSize(0.05,"XYZ");
	gStyle->SetTitleSize(0.05,"T");
	gStyle->SetTitleOffset(1.2,"x");
        gStyle->SetTitleOffset(1.5,"y");
        gStyle->SetTitleOffset(1,"z");
	gStyle->SetLabelOffset(0.01,"xyz");
	gStyle->SetLabelColor(kRed,"XYZ");
        gStyle->SetTitleColor(kRed,"XYZ");
        gStyle->SetTitleColor(kRed,"T");
        gStyle->SetAxisColor(kRed,"XYZ");
        gStyle->SetFrameLineColor(kRed);
	//gStyle->SetFrameFillStyle(3003);
	//gStyle->SetFillColor(kGreen);
        //gStyle->SetFrameFillColor(kGreen);*/
        //gStyle->SetTitleX(0.55);

/*
        gStyle->SetFuncWidth(1);
        gStyle->SetHistLineWidth(1);
        gStyle->SetFuncColor(kRed);
		gStyle->SetEndErrorSize(0);
	//TGaxis::SetMaxDigits(3);

}
*/	
	
