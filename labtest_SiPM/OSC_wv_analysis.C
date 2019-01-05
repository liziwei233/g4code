Double_t TRfun(Double_t *x, Double_t *par) {
        
    double c=par[0]*par[0];
    double sigma_sipm = par[1];
    double sigma_ele = par[2];
    double val=0;
    val = TMath::Sqrt(sigma_sipm*sigma_sipm/x[0]+sigma_ele*sigma_ele+c);
    return val;
}

void OSC_wv_analysis(){
    gStyle->SetOptFit(1111);
    int v0=55;
    int v=0;
    char buff[1024];
    char name[1024];
    const int n=9;
    int npe[6]={5,10,12,20,30,40};
    /*
    ofstream output("timeres.dat");
    for (int s=0;s<n;s++)
    {
        TCanvas *c1 = new TCanvas("c1","",800,600);
        c1->cd();
        v=v0+s;
        //sprintf(name,"%dpe_57v",npe[s]);
        sprintf(name,"hamamatsu_%d",s+1);
        sprintf(buff,"%s.root",name);
        TFile *f1 = new TFile(buff,"READ");
        TTree *t1 = (TTree*)f1->Get("Pico");

        TH1F *h = new TH1F("h","",4000,-10,10);
        h->SetTitle("time resolution");
        TCut cMCPy = "MCP1_global_maximum_y>0.0012";
        TCut cMCPQ = "MCP1_all_charge>0.02";

        TCut cSiPMy = "MCP2_global_maximum_y>0.001";
        TCut cSiPMQ = "MCP2_all_charge>0.59";
        TCut cSiPMQ_3 = "MCP2_all_charge>0.59";


        TCut cfail = "MCP1_twentypercent_failed!=1&&MCP2_twentypercent_failed!=1";
        if (s!=3) t1->Draw("MCP2_twentypercent_time-MCP1_twentypercent_time>>h",cMCPy&&cMCPQ&&cSiPMy&&cSiPMQ&&cfail);
        else t1->Draw("MCP2_twentypercent_time-MCP1_twentypercent_time>>h",cMCPy&&cMCPQ&&cSiPMy&&cSiPMQ_3&&cfail);

        TF1 *fit = new TF1("fit","gaus",-10,10);
        h->Fit(fit);
        float mean = fit->GetParameter(1);
        float sigma = fit->GetParameter(2);
        if(sigma>0.4) h->Rebin(8);
        else if(sigma>0.2) h->Rebin(4);
        
        h->Fit(fit);
        mean = fit->GetParameter(1);
        sigma = fit->GetParameter(2);
        h->GetXaxis()->SetRangeUser(mean-5*sigma,mean+5*sigma);
        h->Fit(fit,"","",mean-sigma,mean+sigma);
        output<<v<<"\t"<<fit->GetParameter(2)<<endl;
        sprintf(buff,"%s_tr.png",name);
        c1->SaveAs(buff);

    }
    output.close();
   */ 
    TCanvas *c2 = new TCanvas("c2","",800,600);
    c2->cd();
    double vop[n];
    double tr[n];
    double yFit[n];
    double yCal[n];
    double par[3];
    double x[n]={2.46,7,14,25,31,96.8,221,492,695};

    //TH2F *gtr = new TH2F("gtr",";#V_{op} (V); TimeRes (ns)",7,54,61,200,0.3,1);
    ifstream input;
    input.open("timeres.dat");
    for(int i=0;i<n;i++){
        input>>vop[i]>>tr[i];
        cout<<vop<<"\t"<<tr<<endl;
    }
    TGraph *gtr = new TGraph(n,x,tr);
    gtr->SetMarkerStyle(20);
    gtr->SetMarkerSize(1.5);
    gtr->SetMarkerColor(2);
    gtr->Draw("AP");
    //return;
    gtr->SetTitle("");
    gtr->GetXaxis()->SetTitle("npe");
    gtr->GetYaxis()->SetTitle("Time Res (ns)");
    gtr->GetYaxis()->SetRangeUser(-0.1,1);
    gtr->GetXaxis()->SetRangeUser(-10,1e3);
    TF1 *myfun = new TF1("myfun",TRfun,0,1000,3);
    myfun->SetParLimits(0,0,0.1);
    myfun->SetParLimits(2,0.005,0.3);
    myfun->SetParameter(1,0.400);
    //myfun->SetParameter(2,0.020);
    myfun->SetParNames("c","#sigma_{SiPM}","#sigma_{ele}");
    gtr->Fit(myfun,"","",20,1e3);
    myfun->SetLineColor(4);
    myfun->Draw("same");
    /*
    par[0]=myfun->GetParameter(0);
    par[1]=myfun->GetParameter(1);
    par[2]=myfun->GetParameter(2);
    for(int i=0;i<n;i++){

    yFit[i]=myfun->Eval(x[i]);
    yCal[i]=sqrt(par[1]*par[1]/x[i]+par[2]*par[2]+par[0]);
    
    cout<<par[0]<<"\t"<<par[1]<<"\t"<<par[2]<<"\t"<<tr[i]<<"\t"<<yFit[i]<<"\t"<<yCal[i]<<endl;
    }
    TGraph *gtF = new TGraph(n,x,yFit);
    gtF->Draw("sameLP");
*/
    //gtr->Draw("AP");
    sprintf(buff,"TR_Distribution.png");
    c2->SaveAs(buff);
}

