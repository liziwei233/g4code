#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>
#include <TRandom.h>

void OPinformation(const char *rootname="")
{
	char name[1024];
	char buff[1024];
	
    gStyle->SetOptStat(1111);

    void DrawMyGraph(TGraph *datagraph, const char *xtitle, const char *ytitle, Float_t MSize=1, Int_t MStyle =28, Int_t MColor=1, Int_t LColor=1, Float_t LWidth=1, Int_t LStyle=1, Int_t FColor=16);
    void DrawMyHist1(TH1 *datahist, const char *xtitle,const char *ytitle, Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);
    void DrawMyHist2(TH2 *datahist, const char *xtitle,const char *ytitle, Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);
    void SetMyPad(TPad *pad,float left, float right, float top, float bottom);
    TLegend* DrawMyLeg(Double_t xlow=0.2, Double_t ylow=0.2, Double_t xup=0.5, Double_t yup=0.5, Int_t textFont=62, Size_t textSize=0.05);
	
	TF1* gausfit(TH1* hU,int rbU,float* U_RL, float* U_RR);
	void twoguasfit(TH1 *ht,double* tRL,double* tRR,double fac=0.4, int rbt=1);
    
	for(int j=0; j<1; j++){
		
		sprintf(name,"%s",rootname);
		sprintf(buff,"%s.root",name);
		
	TFile *f1 = new TFile(buff,"READ");
    TTree *t1 = (TTree*)f1->Get("Run");


    int counts = 0;

    std::vector<int>* fphID;
    std::vector<int>* fphB;

    std::vector<int>* fpmtID[2];  //0==left,1==right;


	
	
	
	
    fphID = new std::vector<int>;
    fphB = new std::vector<int>;
    fpmtID[0] = new std::vector<int>;
    fpmtID[1] = new std::vector<int>;


    t1->SetBranchAddress("PmtL.trackID",&fpmtID[0]);
    t1->SetBranchAddress("PmtR.trackID",&fpmtID[1]);

    t1->SetBranchAddress("ph.ID",&fphID);
    t1->SetBranchAddress("op.Bounce",&fphB);

	

    int N = 0;

	double halfx = 5;
	double halfy = 5;
	
	TCanvas *c1 = new TCanvas("c1","",800,600);
	SetMyPad(c1,0.12,0.05,0.05,0.12);
	
	TCanvas *c2 = new TCanvas("c2","",800,600);
	SetMyPad(c2,0.12,0.05,0.05,0.12);
	
	TCanvas *c3 = new TCanvas("c3","",800,800);
	SetMyPad(c3,0.12,0.05,0.05,0.12);

	TCanvas *c4 = new TCanvas("c4","",800,600);
	SetMyPad(c4,0.12,0.05,0.05,0.12);
	
	
    int B = 0;
    bool flagph=0;
	bool flagdet=0;
    int temp=0;


    //binN = ID+1;
    //TH1D *HB = new TH1D;

    //TH1D *HL = new TH1D("HL","length",200,0,0.4);



    //t1->Draw("op.L>>HL","op.ID==binN");

    ofstream outputdata("PhotonsInform.dat",ios::trunc);

    t1->GetEntry(0);
    N = fphID->size();

    cout<<"the photon quantity is "<<N<<endl;
//return;
    TH1D *HB = new TH1D("HB",";ID;Bounce",N+10,1,N+10);
	
	TH1D *hphB = new TH1D("hphB",";Bounce;Counts",500,1,500);
	
	TH1D *hpmtB[2];
	hpmtB[0] = new TH1D("hpmtBL",";Bounce;Counts",500,1,500);
	hpmtB[1] = new TH1D("hpmtBR",";Bounce;Counts",500,1,500);
	TH1D *hloseB = new TH1D("hotherB",";Bounce;Counts",500,1,500);
	
	
	
	
	
	/****************************************************************
	*
	*		Draw the Bounce details about ph
	*
	*****************************************************************/
	
    t1->Draw("op.Bounce>>HB","Entry$==0");
    for(int i=2; i<N; i++){
        //ID = i+1;
        //binN = ID+1;
        B = HB->GetBinContent(i);
		flagdet=0;
		flagph=0;
		
		temp = fphID->size();
		for(int s=0;s<temp;s++)
		{
			if((*fphID)[s]==i)
			{
				hphB->Fill(B);
				flagph=1;
			}
		}
		
		if (flagph){
		for (int k=0;k<2;k++)
		{
		temp = fpmtID[k]->size();
			for(int s=0;s<temp;s++)
			{
				if((*fpmtID[k])[s]==i)
				{
					hpmtB[k]->Fill(B);
					flagdet=1;
				}
			}
		
		}
			if(!flagdet) hloseB->Fill(B);
		}
		

    }
	c1->cd();
	c1->SetLogx();
	c1->Clear();
	DrawMyHist1(hphB,"Bounce","Counts",1,2,1);
	hphB->Draw();
	sprintf(buff,"%s_BounceStat_allphotons.png",name);
	c1->SaveAs(buff);
	c1->Clear();
	DrawMyHist1(hpmtB[0],"Bounce","Counts",1,2,1);
	hpmtB[0]->Draw();
	sprintf(buff,"%s_BounceStat_PMTleft.png",name);
	c1->SaveAs(buff);
	//return;
	c1->Clear();
	DrawMyHist1(hpmtB[1],"Bounce","Counts",1,2,1);
	hpmtB[1]->Draw();
	sprintf(buff,"%s_BounceStat_PMTright.png",name);
	c1->SaveAs(buff);
	c1->Clear();
	DrawMyHist1(hloseB,"Bounce","Counts",1,2,1);
	hloseB->Draw();
	sprintf(buff,"%s_BounceStat_losingph.png",name);
	c1->SaveAs(buff);
	
	
	
	
	
	
	/****************************************************************
	*
	*		Draw the Hit time about ph
	*
	*****************************************************************/
	TH1D *ht[2];
	ht[0]=new TH1D("htL",";Time (ns); Counts",200,0,15);
	ht[1]=new TH1D("htR",";Time (ns); Counts",200,0,15);
	
	
	DrawMyHist1(ht[0],"Time (ns)","Counts",1,2,1); //Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);
	DrawMyHist1(ht[1],"Time (ns)","Counts",1,2,1); //Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);
	c2->cd();
	c2->Clear();
	t1->Draw("PmtL.t>>htL","Entry$==0");
	sprintf(buff,"%s_hittime_PMTleft.png",name);
	c2->SaveAs(buff);
	c2->Clear();
	t1->Draw("PmtR.t>>htR","Entry$==0");
	sprintf(buff,"%s_hittime_PMTright.png",name);
	c2->SaveAs(buff);
	
	
	
	
	/****************************************************************
	*
	*		Draw the Hit positon about ph
	*
	*****************************************************************/
	TH2D *hpos[2];
	hpos[0]=new TH2D("hposL",";x (mm); y (mm)",100,-1*halfx,halfx,100,-1*halfy,halfy);
	hpos[1]=new TH2D("hposR",";x (mm); y (mm)",100,-1*halfx,halfx,100,-1*halfy,halfy);
	
	c3->cd();
	
	DrawMyHist2(hpos[0],"x (mm)","y (mm)",1,2,1); //Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);
	DrawMyHist2(hpos[1],"x (mm)","y (mm)",1,2,1); //Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);
	c3->Clear();
	t1->Draw("PmtL.x:PmtL.y>>hposL","Entry$==0","colz");
	sprintf(buff,"%s_hitpostion_PMTleft.png",name);
	c3->SaveAs(buff);
	c3->Clear();
	t1->Draw("PmtR.x:PmtR.y>>hposR","Entry$==0","colz");
	sprintf(buff,"%s_hitposition_PMTright.png",name);
	c3->SaveAs(buff);
	
	
	
	
	/****************************************************************
	*
	*		Draw the wavelength about ph
	*
	*****************************************************************/
	
	TH1D *hphwl = new TH1D("hphwl",";Wavelength (nm); Counts",200,180,800);
	TH1D *hwl[2];
	hwl[0]=new TH1D("hwlL",";Wavelength (nm); Counts",200,180,800);
	hwl[1]=new TH1D("hwlR",";Wavelength (nm); Counts",200,180,800);
	
	c4->cd();
	DrawMyHist1(hphwl,"Wavelength (nm)","Counts",1,2,1); //Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);
	DrawMyHist1(hwl[0],"Wavelength (nm)","Counts",1,2,1); //Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);
	DrawMyHist1(hwl[1],"Wavelength (nm)","Counts",1,2,1); //Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);
	c4->Clear();
	t1->Draw("1264/ph.E/1e6>>hphwl","Entry$==0");
	sprintf(buff,"%s_wavelength_allphotons.png",name);
	c4->SaveAs(buff);
	
	c4->Clear();
	t1->Draw("1264/PmtL.E/1e6>>hwlL","Entry$==0");
	sprintf(buff,"%s_wavelength_PMTleft.png",name);
	c4->SaveAs(buff);
	c4->Clear();
	t1->Draw("1264/PmtR.E/1e6>>hwlR","Entry$==0");
	sprintf(buff,"%s_wavelength_PMTright.png",name);
	c4->SaveAs(buff);
	
	
	}
}


void DrawMyGraph(TGraph *datagraph, const char *xtitle, const char *ytitle, Float_t MSize=1, Int_t MStyle =28, Int_t MColor=1, Int_t LColor=1, Float_t LWidth=1, Int_t LStyle=1, Int_t FColor=16){
    datagraph->SetLineColor( LColor );
    datagraph->SetLineWidth( LWidth );
    datagraph->SetLineStyle( LStyle );
    datagraph->SetMarkerSize( MSize );
    datagraph->SetMarkerStyle( MStyle );
    datagraph->SetMarkerColor( MColor );
    datagraph->SetFillColor( FColor );
    //datagraph->SetFillStyle( FStyle );
    datagraph->GetXaxis()->SetTitle( xtitle);
    datagraph->GetYaxis()->SetTitle( ytitle);
    datagraph->GetXaxis()->SetAxisColor(1);
    datagraph->GetYaxis()->SetAxisColor(1);
    datagraph->GetXaxis()->SetLabelColor(1);
    datagraph->GetYaxis()->SetLabelColor(1);
    datagraph->GetXaxis()->SetLabelFont( 42 );
    datagraph->GetYaxis()->SetLabelFont( 42 );
    datagraph->GetXaxis()->SetLabelSize( 0.05 );
    datagraph->GetYaxis()->SetLabelSize( 0.05 );
    datagraph->GetXaxis()->SetLabelOffset( 0.01 );
    datagraph->GetYaxis()->SetLabelOffset( 0.01 );
    datagraph->GetXaxis()->SetTitleFont( 42 );
    datagraph->GetYaxis()->SetTitleFont( 42 );
    //datagraph->GetXaxis()->SetTitleColor( TitleColor);
    //datagraph->GetYaxis()->SetTitleColor( TitleColor );
    datagraph->GetXaxis()->SetTitleSize(0.06);
    datagraph->GetYaxis()->SetTitleSize(0.06);
    datagraph->GetXaxis()->SetTitleOffset(0.8);
    datagraph->GetYaxis()->SetTitleOffset(0.8);
}

void DrawMyHist1(TH1 *datahist, const char *xtitle,const char *ytitle, Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1){
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
void DrawMyHist2(TH2 *datahist, const char *xtitle,const char *ytitle, Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1){
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
void SetMyPad(TPad *pad,float left, float right, float top, float bottom){
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
TLegend* DrawMyLeg(Double_t xlow=0.2, Double_t ylow=0.2, Double_t xup=0.5, Double_t yup=0.5, Int_t textFont=62, Size_t textSize=0.05){
    TLegend *leg = new TLegend(xlow,ylow,xup,yup);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetFillColor(10);
    leg->SetTextFont(textFont);
    leg->SetTextSize(textSize);
    //leg->Draw("same");
    return leg;
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
		fit2->SetParameter(5,3*sigma);
		
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
		
		//if(*tRL<mean-5*sigma)
		//{
			*tRL = mean-5*sigma;
			*tRR = mean+5*sigma;
			//}
		
		cout<<h->GetName()<<"\t"<<*tRL<<"\t"<<*tRR<<endl;
        h->GetXaxis()->SetRangeUser(*tRL,*tRR);
}
