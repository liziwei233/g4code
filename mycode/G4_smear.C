//
//* Environment: 
//**     class LZWfunc
//**     class DrawMyfunc
//* Data: 
//**     Geant4 simulation like TORCH for HIEPA 
//* Function:
//**     Draw results and use different sigmas to smear
//* Date: 2019.3.19
//*
//***********************************
//

#include <string>
#include <time.h>
#include <TString.h>

#include "Include/DrawMyfunc.h"
#include "Include/LZWfunc.h"

using namespace std;

void G4_smear(const char* name=""){
	TF1* gausfit(TH1* ht,int rbU,float* U_RL, float* U_RR);
	
	gStyle->SetOptFit(1111);
	
    char buff[1024];
    //char name[1024];

    int L[5]={10,20,30,40,50};
    double sigma[4]={0,0.050,0.075,0.100};
    int Nsi=sizeof(sigma)/sizeof(sigma[0]);
    int sample[]={5,10,15,20,25,30,35,40,45,50};
    int Nsa=sizeof(sample)/sizeof(sample[0]);
    //float t.L=-0.2,t.R=1.0;
    RANGE t={-0.2,1};
	int bin= (t.R-t.L)/1e-4;
	int pos=1; //responds the specify sigma value;
    int rb=8;
	
	vector<double>* hitTL;
    hitTL = new vector<double>;
    vector<double>* flyTL;
    flyTL = new vector<double>;
    vector<double>* hitTR;
    hitTR = new vector<double>;
    vector<double>* flyTR;
    flyTR = new vector<double>;


    vector<double> deltaT;
    vector<double> pickT;
    deltaT.reserve(2e3);

    double smear;
    double temp;
    double thetime;
    int pickup=0;
    TRandom r,r1;

    ofstream output;
    sprintf(buff,"%s_results.dat",name);
    output.open(buff,ios::trunc);
    
    TCanvas *c1 = new TCanvas("c1","",800,600);
    TCanvas *c2 = new TCanvas("c2","",800,600);
    //TCanvas *c3 = new TCanvas("c3","",1800,600);
    LZWfunc lzw;
	
    TH1D* htime = new TH1D("htime",";time distribution (ns);Counts",bin,t.L,t.R);
    TH1D* hhit = new TH1D("hhit",";time distribution (ns);Counts",bin,t.L,t.R);

    //TH1D* hrdom[4];
    
    //TH1D* hrdom2 = new TH1D("hrdom2","",200,0,1);

    for (int i = 4;i < 5;i++)
    {
        /*
        string time = lzw.getTime();
        output<<"**************"<<endl;
	    output<<"*"<<endl;
    	output<<"* $ Smearing Function $"<<endl;
    	output<<"*"<<endl;
        output<<"* Date: "<< time << endl;
    	output<<"*"<<endl;
    	output<<"* sigma    sampleNums  hittimeMean hittimeSigma"<<endl;
        */

        //sprintf(name,"%s%dmm",path,L[i]);
        //cout<<"file name: "<<name<<endl;
        sprintf(buff,"%s.root",name);
        if(gSystem->AccessPathName(buff))
    {
        cout<<"The root file isn't exist!"<<endl;
      //gROOT->ProcessLine(".x CreateData.C");
        return;
    }
    //output.open(buff,ios::app);
        TFile *f1 = new TFile(buff,"READ");
        TTree *t1 = (TTree*)f1->Get("Run");

        t1->SetBranchAddress("PmtL.t",&hitTL);
        t1->SetBranchAddress("PmtL.flyt",&flyTL);
        t1->SetBranchAddress("PmtR.t",&hitTR);
        t1->SetBranchAddress("PmtR.flyt",&flyTR);

        int N=t1->GetEntries();
        //* smearing
        for(int p=0;p<Nsi;p++){
        //for(int p=0;p<2;p++){
        
        //* Sampling
        for(int s=0; s<Nsa;s++){
        //for(int s=0; s<2;s++){
        hhit->Reset();
        //sprintf(buff,"hr%d",p*2+s);
        //hrdom[p*2+s] = new TH1D(buff,"",2000,-1,2);
        //* Events
        for(int j = 0; j < 8000;j++)
        {

            hitTL->clear();
            flyTL->clear();
            hitTR->clear();
            flyTR->clear();

            vector<double>().swap(deltaT);
            vector<double>().swap(pickT);

            htime->Reset();
            smear=0;
            temp=0;
            thetime=0;

            t1->GetEntry(j);
            temp = hitTL->size();
            //cout<<"* left pmt hits= "<<temp<<endl;
            for( int k = 0;k < temp;k++)
            {
                thetime = (*hitTL)[k]-(*flyTL)[k];                
                //r.SetSeed(k);
                deltaT.push_back(thetime);
                //smear = r.Gaus(0,sigma[pos]);
                //thetime = thetime+smear;
                //htime->Fill(thetime);
            }


            temp = hitTR->size();
            //cout<<"* right pmt hits= "<<temp<<endl;
            //cout<<"* Curren event num = "<<j<<endl;
            for( int k = 0;k < temp;k++)
            {
                thetime = (*hitTR)[k]-(*flyTR)[k];                
                deltaT.push_back(thetime);
                //r.SetSeed(-1-k);
                //smear = r.Gaus(0,sigma[pos]);
                //thetime = thetime+smear;
                //htime->Fill(thetime);
            }
            thetime=0;
            temp=deltaT.size();
            for(int k=0;k<sample[s];k++){
                        //r.SetSeed(k);
                        pickup = r1.Rndm()*temp;
                        
                        //hrdom1->Fill(pickup/temp);
                        //hrdom[p*2+s]->Fill(pickup/temp);
                        //thetime = (*hitTL)[pickup]-(*flyTL)[pickup];                
                        smear = r.Gaus(0,sigma[p]);
                        //hrdom2->Fill(smear);
                        
                        thetime = deltaT.at(pickup)+smear;
                        pickT.push_back(thetime);
                        htime->Fill(thetime);
			}

        
        c1->cd();
        c1->Clear();
        htime->Draw();
		//TF1 *fit1 = gausfit(htime,rb,&t.L,&t.R);
        hhit->Fill(TMath::Mean(pickT.begin(),pickT.end()));
        //output<<htime->GetRMS()<<endl;
        }
        c2->cd();
        c2->Clear();
        hhit->Draw();
		TF1 *fit2 = lzw.gausfit(hhit,rb,5,t);
        
	    
        output<<sigma[p]<<"\t"<<htime->GetEntries()<<"\t"<<fit2->GetParameter(1)<<"\t"<<fit2->GetParameter(2)<<endl;
        //output<<htime->GetEntries()<<"\t"<<fit1->GetParameter(2)<<endl;
		sprintf(buff,"%s_sigma%g_sample%d.png",name,sigma[p],sample[s]);
        c2->SaveAs(buff);
        }
        }

    }
    //output<<"*"<<endl;
	//output<<"********* end *********"<<endl;
	//output<<"\n\n\n"<<endl;

    output.close();
    
    /*c3->Divide(4,1);
    for(int rd=0;rd<4;rd++){
    c3->cd(rd+1);
    hrdom[rd]->Draw();
    }
	sprintf(buff,"%s_random.png",name);
    c3->SaveAs(buff);
    */
}
/*
TF1* gausfit(TH1* ht,int rbU,float* U_RL, float* U_RR){
	
	
	TH1* hU = (TH1*)ht->Clone();
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
*/