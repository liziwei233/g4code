#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>
#include <TRandom.h>

Double_t response(Double_t x, Double_t par[7]){
	Double_t G = par[0];
	Double_t Q = par[1];
	Double_t Ca = par[2];
	Double_t C = par[3];
	Double_t R = par[4];
	Double_t Z = par[5];
	Double_t rise = par[6];

	Double_t val = 0;
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
void SPE(){
    Double_t response(Double_t x, Double_t par[7]);

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

/*
	//MCP
	SPEpar[0]=3e5;  //Gain
	SPEpar[1]=1.6e-19; //e
	SPEpar[2]=3e-12;  //Ca  ??
	SPEpar[3]=6.5e-12; //C
	SPEpar[4]=10e3;    //R   ??
	SPEpar[5]=50;      //Z
	SPEpar[6]=80e-12;  //rise time
*/

    for(int j=0;j<range;j++){
			x[j]=(RR-RL)/range*j+RL;
			y[j]=response(x[j],SPEpar);
		}   

	U = TMath::MinElement(range,y);
	cout<<"the amplitude of signal = "<<U<<endl;

	flag1 = 1;
	flag2 = 1;

	for( int q = 0 ;q < range; q++){
			if(y[q]<U*0.1&& flag1)
			//if(yR[q]<Rate*UR && flagR && yR[q]<thrd) 
			{
				rise1 = x[q];
				flag1 = 0;
				//cout<<"		[+] selected xR = "<<xT0_R<<"\t"<<yR[q]<<endl;
				}
			if(y[q]<U*0.9&&flag2) 
			//if(yL[q]<Rate*UL && flagL && yL[q]<thrd) 
			{
				rise2 = x[q];
				flag2 = 0;
				//cout<<"		[+] selected xL = "<<xT0_L<<"\t"<<yL[q]<<endl;
				}
			if(y[q]>U*0.9&&flag2==0)
			{
				fall2 = x[q];
				flag2 = 1;
			}
			if(y[q]>U*0.1&&flag1==0)
			{
				fall1 = x[q];
				flag1 = 1;
			}
			sum+=y[q];
			

				//cout<<" q value"<<q<<endl;
		}
	tRise = (rise2-rise1)*1e12;
	tFall = (fall1-fall2)*1e12;
	//Q = -1*TMath::Mean(range,y)*range/50*(RR-RL)*1e12; //pC
	Q = -1*sum/50*10; //pC
	Q_act=SPEpar[0]*SPEpar[1]*1e12; //e

	cout<<"rise time = "<<tRise<<"(ps), fall time = "<<tFall<<"(ps)"<<endl;
	cout<<"the charge = "<<Q<<endl;
	cout<<"the charge (actrul) = "<<Q_act<<endl;

    TGraph *g = new TGraph(range,x,y);   
    //TCanvas *c = new TCanvas("c","",1600,600);

	//c->cd();
	g->Draw();
	TAxis *xA = g->GetXaxis();
	TAxis *yA = g->GetYaxis();
	g->SetLineWidth(2);
	g->SetLineColor(kRed);
	xA->SetTitle("Time (s)");
	xA->SetRangeUser(RL,RR);
	//xA->SetRangeUser(5,10);
	yA->SetTitle("Amplitude (V)");

    
}
