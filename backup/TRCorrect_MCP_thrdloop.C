#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>

void TRCorrect(){
	
	//void MygStyle();
	//MygStyle();
	gStyle->SetOptFit(111);
	
	char filename[100];
	char name[100];
	char buff[1024];

	float TH[6]={-5,-7.5,-10,-12.5,-15,-20};
	float certain = 0;
	//float RL[2]={0.25e-9,1.05e-9};
	float R_RL=0.02e-9;
	float R_RR=0.25e-9;
	
	float L_RL=1.1e-9;
	float L_RR=1.32e-9;
	
	float RL_cor = -0.1e-9;
	float RR_cor = 0.1e-9;


	float U_RL=-120.,U_RR=-10.;
	//float U_RL=-230., U_RR=-210.;
	float thrd;
	//float TTS = 175;
	certain = 5;
	for (int s = certain; s<certain+1; s++)
	//for (int s = 0;s<6;s++)
	{
		thrd = TH[s];
		sprintf(filename,"twoside");
		sprintf(name,"%sthrd_%g",filename,TMath::Abs(thrd));
		
		
	sprintf(buff,"analysed_%sthrd_%g.root",filename,TMath::Abs(thrd));
	TFile *f = new TFile(buff,"READ");
	TTree *t = (TTree*)f->Get("data");
	double UL,UR,T0L,T0R,T0;
	double T0L_cor,T0R_cor;
	t->SetBranchAddress("UL",&UL);
	t->SetBranchAddress("UR",&UR);
	t->SetBranchAddress("T0L",&T0L);
	t->SetBranchAddress("T0R",&T0R);
	//t->SetBranchAddress("T0",&T0);
	/*
	
	*/
	TCanvas *c1 = new TCanvas("c1","",1600,600);
	TCanvas *c2 = new TCanvas("c2","",1600,600);
	TCanvas *c3 = new TCanvas("c3","",1600,600);
	//TCanvas *c4 = new TCanvas("c4","",800,600);
	
	TH1D *hL = new TH1D("hL","",12e3,-1.e-9,5.0e-9);
	TH1D *hR = new TH1D("hR","",12e3,-1.e-9,5.0e-9);
	TH1D *hpfx1 = new TH1D("hpfx1","",200,U_RL,U_RR);
	TH1D *hpfx2 = new TH1D("hpfx2","",200,U_RL,U_RR);
	TH2D *htL = new TH2D("htL","",200,U_RL,U_RR,12e3,-1e-9,5.0e-9);
	TH2D *htR = new TH2D("htR","",200,U_RL,U_RR,12e3,-1e-9,5.0e-9);

	
	
	int N = t->GetEntries();

	for(int i=0;i<N;i++){
		t->GetEntry(i);
		hL->Fill(T0L);
		hR->Fill(T0R);
		htL->Fill(UL,T0L);
		htR->Fill(UR,T0R);
	}

	c1->Divide(2,1);
	
	c1->cd(1);
	//t->Draw("T0L>>hL");
	hL->Rebin(4);
	hL->Draw();
	hL->GetXaxis()->SetRangeUser(L_RL,L_RR);
	TF1 *fit1 = new TF1("fit1","gaus",-1e-9,5.0e-9);
	hL->Fit(fit1);


	c1->cd(2);
	//t->Draw("T0R>>hR");
	hR->Rebin(4);
	hR->Draw();
	hR->GetXaxis()->SetRangeUser(R_RL,R_RR);
	TF1 *fit2 = new TF1("fit2","gaus",-1e-9,5.0e-9);
	hR->Fit(fit2);
	
	sprintf(buff,"T0_%s.png",name);
	c1->SaveAs(buff);

	//return;
	
	ofstream output("timeres.dat",ios::app);
	output<<fit1->GetParameter(2)<<"\t"<<fit1->GetParError(2)<<endl;
	output<<fit2->GetParameter(2)<<"\t"<<fit2->GetParError(2)<<endl;



	c2->Divide(2,1);
	c2->cd(1);
	//t->Draw("T0L:UL>>htL","","colz");
	
	htL->RebinX(8);
	htL->RebinY(12);
	htL->Draw("colz");
	htL->GetYaxis()->SetRangeUser(L_RL,L_RR);
	c2->cd(2);
	//t->Draw("T0R:UR>>htR","","colz");
	htR->RebinX(8);
	htR->RebinY(12);
	htR->Draw("colz");
	htR->GetYaxis()->SetRangeUser(R_RL,R_RR);

	sprintf(buff,"TU_%s.png",name);
	c2->SaveAs(buff);
	//return;

	c3->Divide(2,1);
	//htL->GetXaxis()->SetRangeUser(-68,-66);

	htL->ProfileX();
	c3->cd(1);
	hpfx1=(TH1D*)gDirectory->Get("htL_pfx");
	hpfx1->GetYaxis()->SetRangeUser(L_RL,L_RR);
	hpfx1->Draw();
	//return;
	TF1 *fit3 = new TF1("fit3","[0]+[1]/TMath::Sqrt(abs(x))+[2]/abs(x)+[3]/abs(x)/TMath::Sqrt(abs(x))+[4]/abs(x)/abs(x)",U_RL,U_RR);
	//TF1 *fit3 = new TF1("fit3","pol4",U_RL,U_RR);
	//TF1 *fit3 = new TF1("fit3","[0]+[1]/TMath::Sqrt(abs(x))+[2]/abs(x)+[3]/abs(x)/TMath::Sqrt(abs(x))",U_RL,U_RR);
	//fit3->FixParameter(3,1e-18);
	hpfx1->SetStats(0);
	
	hpfx1->Fit(fit3,"R");
//return;

	//htR->GetXaxis()->SetRangeUser(-68,-65);
	htR->ProfileX();
	c3->cd(2);
	hpfx2 = (TH1D*)gDirectory->Get("htR_pfx");
	hpfx2->GetYaxis()->SetRangeUser(R_RL,R_RR);
	hpfx2->GetXaxis()->SetRangeUser(-110,-10);
	hpfx2->Draw();
	
	

	TF1 *fit4 = new TF1("fit4","[0]+[1]/TMath::Sqrt(abs(x))+[2]/abs(x)+[3]/abs(x)/TMath::Sqrt(abs(x))+[4]/abs(x)/abs(x)",U_RL,U_RR);
	//TF1 *fit4 = new TF1("fit4","[0]+[1]/TMath::Sqrt(abs(x))+[2]/abs(x)+[3]/abs(x)/TMath::Sqrt(abs(x))",U_RL,U_RR);
	//TF1 *fit4 = new TF1("fit4","pol3",U_RL,U_RR);
	//fit4->FixParameter(3,1e-18);
	hpfx2->SetStats(0);
	hpfx2->Fit(fit4,"R");
	sprintf(buff,"profileX_%s.png",name);
	c3->SaveAs(buff);
	
	
	
	//return;
	//h->Reset();
	
	
	
	//TH2D *ht_cor= new TH2D("ht_cor","",100,U_RL,U_RR,3600,-1e-9,2.6e-9);
	htL->Reset();
	htR->Reset();

	hpfx1->GetYaxis()->SetRangeUser(RL_cor,RR_cor);
	hpfx1->Reset();
	hpfx1->GetYaxis()->SetRangeUser(RL_cor,RR_cor);
	hpfx2->Reset();
	
	//TH1D *h_cor = new TH1D("h_cor","",3600,-1e-9,2.6e-9);
	hL->Reset();
	hR->Reset();

	int n=t->GetEntries();
	for(int i=0; i<n; i++)
	{
		t->GetEntry(i);
		T0L_cor = fit3->Eval(UL);
		T0R_cor = fit4->Eval(UR);
		//cout<<"T0-T0_cor="<<T0<<"-"<<T0_cor<<"="<<T0-T0_cor<<endl;
		
		hL->Fill(T0L-T0L_cor);
		hR->Fill(T0R-T0R_cor);

		htL->Fill(UL,T0L-T0L_cor);
		htR->Fill(UR,T0R-T0R_cor);
		//cout<<"U:T0-T0_cor="<<UR<<"\t"<<T0-T0_cor<<endl;
	}

	
	c1->Clear();
	c1->Divide(2,1);
	//return;
	c1->cd(1);
	hL->Draw();
	hL->GetXaxis()->SetRangeUser(RL_cor,RR_cor);
	//fit1->SetRange(-1e-9,1e-9);
	hL->Fit(fit1);

	c1->cd(2);
	hR->Draw();
	hR->GetXaxis()->SetRangeUser(RL_cor,RR_cor);
	fit2->SetRange(-1e-9,1e-9);
	hR->Fit(fit2);
	
	sprintf(buff,"T0_AfterCorrect_%s.png",name);
	c1->SaveAs(buff);
	//return;
	//output("timeres.dat",ios::app);
	output<<fit1->GetParameter(2)<<"\t"<<fit1->GetParError(2)<<endl;
	output<<fit2->GetParameter(2)<<"\t"<<fit2->GetParError(2)<<endl;

	
	c2->Clear();
	c2->Divide(2,1);
	c2->cd(1);
	htL->GetYaxis()->SetRangeUser(RL_cor,RR_cor);
	htL->Draw("colz");
	
	c2->cd(2);
	htR->GetYaxis()->SetRangeUser(RL_cor,RR_cor);
	htR->Draw("colz");
	//return;
	sprintf(buff,"TU_Correct_%s.png",name);
	c2->SaveAs(buff);

	
	
	c3->Clear();
	c3->Divide(2,1);
	c3->cd(1);
	htL->ProfileX();
	hpfx1=(TH1D*)gDirectory->Get("htL_pfx");
	hpfx1->GetYaxis()->SetRangeUser(RL_cor,RR_cor);
	hpfx1->Draw();
	hpfx1->Fit(fit3,"R");

	c3->cd(2);
	htR->ProfileX();
	hpfx2=(TH1D*)gDirectory->Get("htR_pfx");
	hpfx2->GetYaxis()->SetRangeUser(RL_cor,RR_cor);
	hpfx2->Draw();
	hpfx2->Fit(fit4,"R");

	sprintf(buff,"profileX_AfterCorrect_%s.png",name);
	c3->SaveAs(buff);
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
	
