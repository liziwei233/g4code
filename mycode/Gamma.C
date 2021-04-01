void Gamma(double r,double x0, double beta)
{
    TF1 *f1 = new TF1("f1","[0]*TMath::GammaDist(x,[1],[2],[3])",0,10);
    f1->SetParameters(1e3,r,x0,beta);
    TH1F *h = new TH1F("h","",1e3,0,10);
    h->FillRandom("f1",1e4);
    h->Draw();
    /*
    for(int i=0;i<1e4;i++)
    {

    double r = f1->GetRandom();

    }
    */

}