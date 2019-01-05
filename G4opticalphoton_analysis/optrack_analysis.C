#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>
#include <TRandom.h>

void optrack_analysis()
{
    void DrawMyGraph(TGraph *datagraph, char *xtitle, char *ytitle, float MSize=1, int MStyle =28, int MColor=1, int LColor=1, float LWidth=1, int LStyle=1, int FColor=16);
    TFile *f1 = new TFile("CRTest.root","READ");
    TTree *t1 = (TTree*)f1->Get("Run");


    int counts = 0;

    std::vector<int>* fID;
    std::vector<double>* fGt;
    std::vector<double>* fL;
    std::vector<int>* fBounce;
    std::vector<double>* fWaveL;


    std::vector<double>* ftR;
    std::vector<double>* ftL;
    std::vector<int>* fIDR;
    std::vector<int>* fIDL;

    fID = new std::vector<int>;
    fGt = new std::vector<double>;
    fL = new std::vector<double>;
    fWaveL = new std::vector<double>;
    fBounce = new std::vector<int>;
    ftR = new std::vector<double>;
    ftL = new std::vector<double>;
    fIDL = new std::vector<int>;
    fIDR = new std::vector<int>;


    t1->SetBranchAddress("PmtR.t",&ftR);
    t1->SetBranchAddress("PmtL.t",&ftL);
    t1->SetBranchAddress("PmtR.trackID",&fIDR);
    t1->SetBranchAddress("PmtL.trackID",&fIDL);
    t1->SetBranchAddress("op.Gt",&fGt);
    t1->SetBranchAddress("op.L",&fL);
    t1->SetBranchAddress("op.WaveL",&fWaveL);
    t1->SetBranchAddress("op.ID",&fID);
    t1->SetBranchAddress("op.Bounce",&fBounce);

    t1->SetBranchAddress("op.crkov",&counts);

    int N = 0;

    int ID = 0,B = 0,binN = 0;
    double L = 0;// path length (cm)
    double t = 0; // hit time (ns)
    double wvL=0; // waveLength (nm)
    bool flag=0;
    int temp=0;
    int detected = 0; //val=0,undetected; val=1,detected by L; val=2,detected by R;

    //binN = ID+1;
    //TH1D *HB = new TH1D;

    //TH1D *HL = new TH1D("HL","length",200,0,0.4);



    //t1->Draw("op.L>>HL","op.ID==binN");

    ofstream outputdata("PhotonsInform.dat",ios::trunc);

    t1->GetEntry(0);
    N = counts;

    cout<<"the photon quantity is "<<N<<endl;

    TH1D *HB = new TH1D("HB","bounce",N+1,0,N+1);

    t1->Draw("op.Bounce>>HB");

    for(int i=0; i<N; i++){
        ID = i+1;
        binN = ID+1;
        B = HB->GetBinContent(i+1);
        //cout<<B<<endl;
        //return;
        //continue;
        /*
           temp=fBounce->size();
           for(int j=0;j<temp;j++){
           if((*fBounce)[j]==ID) B++;
           }
         */
        temp=fL->size();
        flag = 1;
        L = 0;
        t = 0;
        for(int k = 0;k<temp;k++)
        {
            if((*fID)[k]==ID) 
            {
                if(flag)    {wvL=(*fWaveL)[k];flag=0;}
                if (t <= (*fGt)[k]) t = (*fGt)[k];
                L+=(*fL)[k];
            }
        }

        detected = 0;
        //t = L*1.58/3e8/100*1e9;

        temp=fIDL->size();
        for(int s=0;s<temp;s++){
            if ((*fIDL)[s]==ID) {
                detected =1;
                //t = (*ftL)[s];
                //cout<<"hit time = "<<t<<endl;
            }
        }

        temp=fIDR->size();
        for(int q=0;q<temp;q++){
            if ((*fIDR)[q]==ID) {
                detected =2;
                //t = (*ftR)[q];
            }
        }
        outputdata<<ID<<"\t"<<wvL<<"\t"<<B<<"\t"<<L<<"\t"<<t<<"\t"<<detected<<endl;



    }





}
