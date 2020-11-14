void DrawWaveform(int id = 467)
{
    int CHid1 = 2;
    int CHid2 = 4;
    char buff[1024];
    const char *name = "testbeam_dumpfile.root";
    TFile *f1 = new TFile(name,"read");
    sprintf(buff,"CH%dgraph_%d",CHid1,id);
TGraph *g1 = (TGraph*)f1->Get(buff);
    sprintf(buff,"CH%dgraph_%d",CHid2,id);
TGraph *g2 = (TGraph*)f1->Get(buff);

    g1->SetLineWidth(3);
    g2->SetLineWidth(3);
    g2->SetLineColor(2);

    g1->GetXaxis()->SetRangeUser(72,82);
    g1->Draw();
    g2->Draw("same");
}
void DrawAveWaveform(const char* name="16chbase-HV2200-SPE-half-D-A7A4-4average")
{
    int CHid1 = 2;
    int CHid2 = 4;
    char buff[1024];
    //const char *path = "/mnt/f/R10754/KT0881/A7A8";
    const char *path = "/mnt/f/R10754/KT0881/A7light_A4crosstalk";
    //const char *path = "/mnt/f/R10754/KT0881/A4light_A7crosstalk";
    sprintf(buff,"%s/%s.root",path,name);
    if(gSystem->AccessPathName(buff))  
    {
        cout<<"Error!! The File "<<buff<<" doesn't exist"<<endl;
      return;
    }
    TFile *f1 = new TFile(buff,"read");
    
    TGraph *g1 = (TGraph*)f1->Get("CH2average_waveform");
    TGraph *g2 = (TGraph*)f1->Get("CH4average_waveform");

    g1->SetLineWidth(3);
    g2->SetLineWidth(3);

    g1->GetXaxis()->SetRangeUser(46,56);
    g2->GetXaxis()->SetRangeUser(46,56);
    new TCanvas();
    g1->Draw();
    sprintf(buff,"%s/%sWFsignal.png",path,name);
    gPad->SaveAs(buff);

    new TCanvas();
    g2->Draw();
    sprintf(buff,"%s/%sWFct.png",path,name);
    gPad->SaveAs(buff);

    g2->SetLineColor(2);
    new TCanvas();
    g1->Draw();
    g2->Draw("same");
    sprintf(buff,"%s/%sWFtogether.png",path,name);
    gPad->SaveAs(buff);
}