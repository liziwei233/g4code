//************************************************
//* Date: 2019.3.19
//* Environment: 
//**     class LZWfunc
//**     class DrawMyfunc
//* Data: 
//*      .root file of Geant4 simulation for Trigger of 
//*      Cosmicray test system.
//* Function:
//**     Caculate the time resolution
//**     Draw augular distribution
//*
//***********************************
//*

#include "Include/DrawMyfunc.h"
#include "Include/LZWfunc.h"
#include "TVector3.h"
using namespace std;

#define DRAWTR



void Trigger(){
    gStyle->SetOptFit(111);
    gStyle->SetOptStat("emr");
    TGaxis::SetMaxDigits(3);
    //gStyle->SetO
    char name[100];
    char buff[1024];

    //***************************************************//
    //--------------Configuration-----------------------//
    //***************************************************//
    
    RANGE t = {0e-9, 1.5e-9};
    RANGE u = {-1e4, 0};
    int bint = (t.R - t.L) / 1e-13;
    int binu = (u.R-u.L)/30;
    int rbt=1;
    double fac=3;
    /*
    RANGE posx = {-20, 20};
    RANGE posy = posx;
    int rbt = 2;
    int rbu = 16;

    // the range set After Correct
    //double L2 = -100e-12;
    //double R2 = 100e-12;
    int bint = (t.R - t.L) / 1e-12;
    int binu = (u.R - u.L) / 1;
    int binpx = (posx.R - posx.L) / 1;
    int binpy = (posy.R - posy.L) / 1;
    const int iter = 5; //the number of iteration
    double fac = 2;
    */
/* 
   EVENT PMTL,PMTR;
   vector<EVENT*> ch;
   ch.push_back(&PMTL);
   ch.push_back(&PMTR);

   charRANGE rL={{-5e3,0},
   {-1e-9,2e-9},
   {0,200}
   };
   charRANGE rR=rL;
   vector<charRANGE> range;
   range.push_back(rL);
   range.push_back(rR);

   CUT cL=rL;
   CUT cR=rR;
   vector<CUT> cut;
   cut.push_back(cL);
   cut.push_back(cR);

   OPTION opt={1,1,4,3};
*/

   double theta,phi,ag_g,ag_det;
   double PEtL,PEtR;
   double PET0;
   double UL,UR;
    //---------------------------------------------------//
    //***************************************************//

    cout << "Start excuate TA correction procedure =====>>>>" << endl;
    //return;

    //cout<<"<<---- Succeed excuating ---->>"<<endl;
    //return;
    //thrd = (s+1)*0.2;
    char rootname[1024]="CRYextenddata";
    sprintf(name, "%s", rootname);
    sprintf(buff, "%s.root", name);
    TFile *f = new TFile(buff, "READ");
    TTree *tr = (TTree *)f->Get("data");
   
    POSITION PCR; //momentum direction of CR
    DIRECTION LCR; //local position of CR

    TVector3 Pdet(0,-1,0); // towards the PmtL detector
    TVector3 Pg(0,-1,0);    // towards ground
    TVector3 Pcr;
    TVector3 Lcr;

    tr->SetBranchAddress("UL", &UL);
    tr->SetBranchAddress("UR", &UR);
    tr->SetBranchAddress("T01stpeL", &PEtL);
    tr->SetBranchAddress("T01stpeR", &PEtR);
    /*
    tr->SetBranchAddress("UL", &PMTL.Amp);
    tr->SetBranchAddress("UR", &PMTR.Amp);
    //tr->SetBranchAddress("T0L", &event[0].time);
    //tr->SetBranchAddress("T0R", &event[1].time);
    tr->SetBranchAddress("T01stpeL", &PMTL.time);
    tr->SetBranchAddress("T01stpeR", &PMTR.time);
    tr->SetBranchAddress("npeL", &PMTL.npe);
    tr->SetBranchAddress("npeR", &PMTR.npe);
    */
    tr->SetBranchAddress("InX", &LCR.x);
    tr->SetBranchAddress("InY", &LCR.y);
    tr->SetBranchAddress("InZ", &LCR.z);
    
    tr->SetBranchAddress("InpX", &PCR.x);
    tr->SetBranchAddress("InpY", &PCR.y);
    tr->SetBranchAddress("InpZ", &PCR.z);

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);

/*
    c1->cd();
    sprintf(buff,"%s_LR",name);
    lzw.SetADD(1);
    lzw.CH2Correction(tr,ch,range,cut,opt,name);
*/


    TH1D *hL = new TH1D("hL", ";Time (s);Counts", bint, t.L, t.R);
    TH1D *hR = new TH1D("hR", ";Time (s);Counts", bint, t.L, t.R);
    TH1D *ht = new TH1D("ht", "", bint, t.L, t.R);
    TH2D *hAT = new TH2D("hAT", "",binu, u.L,u.R, bint, t.L, t.R);
    
    TH1D *hAg1 = new TH1D("hAg1",";#theta (#circ);Counts",200,0,180);
    TH1D *hAg2 = new TH1D("hAg2",";#phi (#circ);Counts",200,-180,0);
    TH1D *hAg3 = new TH1D("hAg3",";#Ag_{det} (#circ);Counts",200,0,180);
    TH1D *hAg4 = new TH1D("hAg4",";#Ag_{g} (#circ);Counts",200,0,180);
    //TH1D *hLPR = new TH1D("hLPR", ";Time (s);Counts", bint, t.L, t.R);
    //TH1D *hLMR = new TH1D("hLMR", ";Time (s);Counts", bint, t.L, t.R);
    //TH2D *hpos = new TH2D("hpos", ";x (mm);y (mm)", binpx, posx.L, posx.R, binpy, posy.L, posy.R)
    LZWfunc lzw;
    DrawMyfunc draw;
    draw.Hist(hAg1,"#theta (#circ)","Counts",2,1,1);
    draw.Hist(hAg2,"#phi (#circ)","Counts",2,1,1);
    draw.Hist(hAg3,"Ag_{det} (#circ)","Counts",2,1,1);
    draw.Hist(hAg4,"Ag_{g} (#circ)","Counts",2,1,1);
    
    draw.Hist(hR,"Time (s)","Counts",2,1,1);
    draw.Hist(hL,"Time (s)","Counts",2,1,1);
    draw.Hist(ht,"Time (s)","Counts",2,1,1);

    int N = 0;
    N = tr->GetEntries();
    cout << "the Entries is :" << N << endl;
    for (int i = 0; i < N; i++)
    {
        tr->GetEntry(i);
        
        Pcr.SetXYZ(PCR.x,PCR.y,PCR.z);
        Lcr.SetXYZ(LCR.x,LCR.y,LCR.z);
        
        theta=Pcr.Theta()*TMath::RadToDeg();
        phi=Pcr.Phi()*TMath::RadToDeg();
        ag_det=Pcr.Angle(Pdet)*TMath::RadToDeg();
        ag_g=Pcr.Angle(Pg)*TMath::RadToDeg();

        hAg1->Fill(theta);
        hAg2->Fill(phi);
        hAg3->Fill(ag_det);
        hAg4->Fill(ag_g);
        //if(ag_det<28&&PEtR>0&&PEtL>0) 
        if(ag_g<2&&PEtR>0)
        {
            hR->Fill(PEtR);
            hL->Fill(PEtL);
            PET0=(PEtR+PEtL)/2;
            ht->Fill(PET0);
            hAT->Fill(UR,PET0);
        } 

    }
    //.x hR->Rebin(10);
#ifdef DRAWTR
    double max=0;
    double length=6e-12; 
    double std=0;
    c1->cd();
    draw.SetPad(gPad,0.12,0.14,0.1,0.1);
    ht->Draw();
    max = ht->GetBinCenter(ht->GetMaximumBin());
    ht->GetXaxis()->SetRangeUser(max-fac*length,max+fac*length);

/*    TF1 *fit = lzw.gausfit(ht,1,fac,t);
    //return;
    double TR=0;
    if(fit) TR = fit->GetParameter(2);*/
    
    double E=0;
    E= ht->GetEntries()/N;
    std = ht->GetStdDev();
    
    sprintf(buff,"Std= %.2fps",std*1e12);
    TLatex* l=draw.Latex(0.6,0.3,buff); 
    l->Draw();
    sprintf(buff,"E= %.2f%%",E*100);
    l->DrawLatex(0.6,0.4,buff); 
    sprintf(buff,"two.png");
    c1->SaveAs(buff);

    c2->cd();
    hAT->Draw("colz");
    sprintf(buff,"twoTArelationship.png");
    c2->SaveAs();

    c3->cd();
    draw.SetPad(gPad,0.12,0.14,0.1,0.1);
    hR->Draw();
    max = hR->GetBinCenter(hR->GetMaximumBin());
    hR->GetXaxis()->SetRangeUser(max-fac*length,max+fac*length);
    //TF1 *fit2 = lzw.gausfit(hR,rbt,fac,t);
    //if(fit2) TR= fit2->GetParameter(2);
    
     E= hR->GetEntries()/N;
    std = hR->GetStdDev();
    
    sprintf(buff,"Std= %.2fps",std*1e12);
    TLatex* l2=draw.Latex(0.6,0.3,buff); 
    l2->Draw();
    sprintf(buff,"E= %.2f%%",E*100);
    l2->DrawLatex(0.6,0.4,buff); 
    sprintf(buff,"singleRight.png");
    c3->SaveAs(buff);
    
    c4->cd();
    draw.SetPad(gPad,0.12,0.14,0.1,0.1);
    hL->Draw();
    max = hL->GetBinCenter(hL->GetMaximumBin());
    hL->GetXaxis()->SetRangeUser(max-fac*length,max+fac*length);
    //TF1 *fit3 = lzw.gausfit(hL,rbt,fac,t);
    //if(fit3) TR= fit3->GetParameter(2);
    
     E= hL->GetEntries()/N;
    std = hL->GetStdDev();
    
    sprintf(buff,"Std= %.2fps",std*1e12);
    TLatex* l3=draw.Latex(0.6,0.3,buff); 
    l3->Draw();
    sprintf(buff,"E= %.2f%%",E*100);
    l3->DrawLatex(0.6,0.4,buff); 
    sprintf(buff,"singleLeft.png");
    c4->SaveAs(buff);
#else    

    c1->cd();
    draw.SetPad(gPad,0.12,0.14,0.08,0.08);
    //draw.SetYTitleOffset(1.1);
    hAg1->Draw();
    sprintf(buff,"%s_theta.png",name);
    c1->SaveAs(buff);
    c2->cd();
    draw.SetPad(gPad,0.12,0.14,0.08,0.08);
    hAg2->Draw();
    sprintf(buff,"%s_phi.png",name);
    c2->SaveAs(buff);
    c3->cd();
    draw.SetPad(gPad,0.12,0.14,0.08,0.08);
    hAg3->Draw();
    sprintf(buff,"%s_angledet.png",name);
    c3->SaveAs(buff);
    c4->cd();
    draw.SetPad(gPad,0.12,0.14,0.08,0.08);
    hAg4->Draw();
    sprintf(buff,"%s_anglegrounf.png",name);
    c4->SaveAs(buff);
#endif

}