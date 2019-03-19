#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>
#include <TRandom.h>

Double_t response(Double_t* xx, Double_t* par){
	double x=xx[0];
	Double_t G = par[0];
	Double_t Q = par[1];
	Double_t Ca = par[2];
	Double_t C = par[3];
	Double_t R = par[4];
	Double_t Z = par[5];
	Double_t rise = par[6];

	Double_t val = 0.;
	Double_t a,b,beta,gamma;
	
	// how to accumulate these single photon signal?
	// which distribution does t0 sample? 
	//!! get t0 by result of simulation!!
	
	beta = ((R+Z)*C+2.0*Ca*R)/(2.0*C*Ca*R*Z);
	gamma = 1.0/(2.0*C*Ca*R*Z);
	a = -(beta+TMath::Sqrt(beta*beta-4.0*gamma))/2.0;
	b = -(beta-TMath::Sqrt(beta*beta-4.0*gamma))/2.0;
	
	val = -0.5*G*Q*(b*TMath::Exp(b*x+0.5*b*b*rise*rise)*TMath::Erfc((-b*rise-x/rise)/TMath::Sqrt(2.0))-a*TMath::Exp(a*x+0.5*a*a*rise*rise)*TMath::Erfc((-a*rise-x/rise)/TMath::Sqrt(2.0)))/(Ca*(b-a));
	
	
	
	return val;


	
}
Double_t rsps(Double_t* xx, Double_t* par){
	double x=xx[0];
	
	Double_t G = par[0];
	Double_t Ca = par[1];
	Double_t C = par[2];
	Double_t R = par[3];
	Double_t rise = par[4];
	
	Double_t Q = 1.6e-19;
	Double_t Z = 50;

	double t=par[5];
	
	/*
	SPEpar[0]=1e6;  //Gain
	SPEpar[2]=2.65e-12;  //Ca  ??
	SPEpar[3]=12e-12; //C
	SPEpar[4]=10e3;    //R   ??
	SPEpar[6]=90e-12;  //rise time
*/

	Double_t val = 0.;
	Double_t a,b,beta,gamma;
	
	// how to accumulate these single photon signal?
	// which distribution does t0 sample? 
	//!! get t0 by result of simulation!!
	
	beta = ((R+Z)*C+2.0*Ca*R)/(2.0*C*Ca*R*Z);
	gamma = 1.0/(2.0*C*Ca*R*Z);
	a = -(beta+TMath::Sqrt(beta*beta-4.0*gamma))/2.0;
	b = -(beta-TMath::Sqrt(beta*beta-4.0*gamma))/2.0;
	
	val = -0.5*G*Q*(b*TMath::Exp(b*(x-t)+0.5*b*b*rise*rise)*TMath::Erfc((-b*rise-(x-t)/rise)/TMath::Sqrt(2.0))-a*TMath::Exp(a*(x-t)+0.5*a*a*rise*rise)*TMath::Erfc((-a*rise-(x-t)/rise)/TMath::Sqrt(2.0)))/(Ca*(b-a));
	
	
	
	return val;


	
}
void SPE(){
	Double_t rsps(Double_t* xx, Double_t* par);
    Double_t response(Double_t* xx, Double_t* par);
    

    void DrawMyGraph(TGraph *datagraph, const char *xtitle, const char *ytitle, Float_t MSize=1, Int_t MStyle =28, Color_t MColor=1, Color_t LColor=1, Float_t LWidth=1, Int_t LStyle=1, Color_t FColor=16);
    void SetMyPad(TVirtualPad *pad,float left, float right, float top, float bottom);
    TLegend* DrawMyLeg(Double_t xlow=0.2, Double_t ylow=0.2, Double_t xup=0.5, Double_t yup=0.5, Int_t textFont=62, Size_t textSize=0.05);
    Double_t RL = -5e-9;
	Double_t RR = 55e-9;
	Int_t range = 6e3;
    Double_t SPEpar[7];
    Double_t x[6000]={};
	Double_t y[6000]={};
	Double_t U = 0;

	double rise1=0,rise2=0;
	double fall1=0,fall2=0;
	double tRise=0,tFall=0;
	double sum=0,Q=0,Q_act=0;
	bool flag1 = 0,flag2 = 0;

 /*   
 	//SiPM
    SPEpar[0]=2.4e6;  //Gain
	SPEpar[1]=1.6e-19; //e
	//SPEpar[2]=150e-12;  //Ca  ??
	//SPEpar[3]=3e-12; //C
	SPEpar[2]=26e-11;  //Ca  ??
	SPEpar[3]=55e-11; //C
	SPEpar[4]=10e3;    //R   ??
	SPEpar[5]=50;      //Z
	SPEpar[6]=0.5e-9;  //rise time
*/

	//MCP R10754
	//rise time 180ps
	//fall time 420ps
	//regulation: Ca up!   >>  Q down! Tfall down!
	//	      C  down! >>  Q down! Tfall down!
	SPEpar[0]=1e6;  //Gain
	SPEpar[1]=1.6e-19; //e
	SPEpar[2]=2.65e-12;  //Ca  ??
	SPEpar[3]=12e-12; //C
	SPEpar[4]=10e3;    //R   ??
	SPEpar[5]=50;      //Z
	SPEpar[6]=90e-12;  //rise time


    for(int j=0;j<range;j++){
			x[j]=(RR-RL)/range*j+RL;
			y[j]=response(&x[j],SPEpar);
                        x[j]=x[j]*1e9;
                        y[j]=y[j]*1e3;
		}   

	U = TMath::MinElement(range,y);
	//cout<<"the amplitude of signal = "<<U*1000<<"mV"<<endl;
	cout<<"the amplitude of signal = "<<U<<"mV"<<endl;

	flag1 = 1;
	flag2 = 1;

	for( int q = 0 ;q < range; q++){
                    if(q>1){
			if(y[q]<U*0.1&& y[q-1]>U*0.1&&y[q-2]>U*0.1)
			//if(yR[q]<Rate*UR && flagR && yR[q]<thrd) 
			{
				rise1 = x[q];
				//flag1 = 0;
				//cout<<"		[+] selected xR = "<<xT0_R<<"\t"<<yR[q]<<endl;
				}
			if(y[q]<U*0.9&&y[q-1]>U*0.9&&y[q-2]>U*0.9) 
			//if(yL[q]<Rate*UL && flagL && yL[q]<thrd) 
			{
				rise2 = x[q];
				//flag2 = 0;
				//cout<<"		[+] selected xL = "<<xT0_L<<"\t"<<yL[q]<<endl;
				}
			if(y[q]>U*0.9&&y[q-1]<U*0.9&&y[q-2]<U*0.9)
			{
				fall2 = x[q];
				//flag2 = 1;
			}
			if(y[q]>U*0.1&&y[q-1]<U*0.1&&y[q-2]<U*0.1)
			{
				fall1 = x[q];
				//flag1 = 1;
			}
                        }
			sum+=y[q];
			

				//cout<<" q value"<<q<<endl;
		}
	tRise = (rise2-rise1)*1e3;
	tFall = (fall1-fall2)*1e3;
	//Q = -1*TMath::Mean(range,y)*range/50*(RR-RL)*1e12; //pC
	Q = (x[0]-x[1])*sum/50; //pC
	Q_act=SPEpar[0]*SPEpar[1]*1e12; //e

	cout<<"rise time = "<<tRise<<"(ps), fall time = "<<tFall<<"(ps)"<<endl;
	cout<<"the charge = "<<Q<<"pC"<<endl;
	cout<<"the charge (actual) = "<<Q_act<<"pC"<<endl;

    TGraph *g = new TGraph(range,x,y);   
    TCanvas *c = new TCanvas("c","",800,600);
    DrawMyGraph(g,"Time (ns)","Amplitude (mV)",1,28,1,kRed,3,1,16);
    SetMyPad(c,0.18,0.1,0.1,0.15);
	c->cd();
	g->Draw("AL");
       g->GetXaxis()->SetRangeUser(-1,3); 
        /*
	TAxis *xA = g->GetXaxis();
	TAxis *yA = g->GetYaxis();
	g->SetLineWidth(2);
	g->SetLineColor(kRed);
	xA->SetTitle("Time (s)");
	xA->SetRangeUser(RL,RR);
	//xA->SetRangeUser(5,10);
	yA->SetTitle("Amplitude (V)");
*/
    
}

void DrawMyGraph(TGraph *datagraph, const char *xtitle, const char *ytitle, Float_t MSize=1, Int_t MStyle =28, Color_t MColor=1, Color_t LColor=1, Float_t LWidth=1, Int_t LStyle=1, Color_t FColor=16){
     datagraph->SetTitle("");
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
     datagraph->GetXaxis()->SetTitleOffset(1.0);
     datagraph->GetYaxis()->SetTitleOffset(1.2);
     datagraph->GetXaxis()->SetNdivisions(510);     
     datagraph->GetYaxis()->SetNdivisions(505);
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
