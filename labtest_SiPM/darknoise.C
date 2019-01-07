void darknoise(int ch=2,int Nfiles=5,const char* mark="trace"){

    TFile* f1 = new TFile("darknoise.root","RECREATE");
    TTree* t1 = new TTree("t1","");
    TTree* t2 = new TTree("t2","");
    float T1,A1;
    vector<float> T,A;
    vector<int>* id;
    id = new vector<int>;
    int DN=0;
    float DNR=0.0;
    float th=-20e-3;
    //float ta2,A2;
    //float tb1,B1;
    char comma;
    char sentence[1024];
    char histname[1024];
    char buff[1024];

    //TCanvas* c1=new TCanvas("c1","c1",800,600);
    //TCanvas* c2=new TCanvas("c2","c2",800,600);
    //TCanvas* c3=new TCanvas("c3","c3",800,600);

    int cFile=0;
    int cPoint=0;
    t1->Branch("T1",&T1,"T1/F");
    //t1->Branch("ta2",&ta2,"ta2/F");
    //t1->Branch("tb1",&tb1,"tb1/F");
    t1->Branch("A1",&A1,"A1/F");
    //t1->Branch("A2",&A2,"A2/F");
    //t1->Branch("B1",&B1,"B1/F");
    t2->Branch("DNR",&DNR,"DNR/F");
    t2->Branch("id",id);
    //TH1F *hA1[200];
    //TH1F *hA2[200];
    //TH1F *hB1[200];

    ifstream input1,input2,input3;
    for(int k=1;k<Nfiles;k++){
        sprintf(buff,"C%d%s%05d.txt",ch,mark,k);
        input1.open(buff);


        if(!input1){
            cout<<buff<<" file can't be found!"<<endl;
            return;
        }	
        //while(input1&&(!input1.eof())&&counter<3){
        cPoint=0;  //initial Point Counter

        cout<<buff<<" Read File : "<<buff<<endl;
        
        vector<float>().swap(T);
        vector<float>().swap(A);
        while(input1&&(!input1.eof())){
            //sprintf(histname,"hmcp_%d",counter);
            //hA1[counter]= new TH1F(histname,";T (25 ps);Amp (V)",20002,0,20002);//25ps/point 
            //sprintf(histname,"hB1_%d",counter);
            //hA2[counter]= new TH1F(histname,";T (25 ps);Amp (V)",20002,0,20002);//25ps/point
            //sprintf(histname,"hB2_%d",counter);
            //hB1[counter]= new TH1F(histname,";T (25 ps);Amp (V)",20002,0,20002);//25ps/point

            //hA1[counter]->Reset();
            if(cPoint<6)
            {    input1>>sentence;
            //cout<<sentence<<endl;
            }
            //input2>>sentence;
            //input3>>sentence;
            //for(int i=0;i<20002;i++){
            else	
            {
                input1>>T1>>comma>>A1;
                T.push_back(T1);
                A.push_back(A1);
          //      cout<<T1<<"\t"<<A1<<endl;
            }
            //input2>>ta2>>comma>>A2;
            //input3>>tb1>>comma>>B1;
            //hA1[counter]->SetBinContent(i,A1);
            //hA2[counter]->SetBinContent(i,A2);
            //hB1[counter]->SetBinContent(i,B1);
            t1->Fill();
            cPoint++;
        }
        input1.close();
        DN=0;
        vector<int>().swap(*id);
        for(int i=5;i<T.size();i++){
            if(A.at(i)<th&&A.at(i-1)<th&&A.at(i-2)<th&&A.at(i-3)>th&&A.at(i-4)>th&&A.at(i-5)>th)
                //Set 2.5ns deadtime 
           {     
            id->push_back(i);
             if(*(id->end()-1)-*(id->end()-2)>50)   DN++;
            }
        }


        //DNR=DN*(1./(T.size()*50e-12))/1e6;
        DNR=DN*(1./10e-6)/1e6;
        t2->Fill();
        //t1->Fill();
        //c1->cd();

        //hA1[counter]->SetLineColor(2);
        //hA1[counter]->Draw();

        //hA2[counter]->SetLineColor(6);
        //hA2[counter]->Draw("same");

        //hB1[counter]->SetLineColor(9);
        //hB1[counter]->Draw("same");
        cout<<"The Dark Noise Ratio is:"<<DNR<<"M cps"<<endl;
        cFile++;

        }
        //input2.close();
        //input3.close();

        
        t1->Write();
        t2->Write();
        TH1F* ht=new TH1F("ht","ht",2000,0,100);
        t2->Draw("DNR>>ht");
        cout<<"DNR's mean value is "<<ht->GetMean()<<"M cps"<<endl;
        
        /*
           if(f1->IsOpen()){
           for(int s=0;s<counter;s++){
           hA1[s]->Write();
           hA2[s]->Write();
           hB1[s]->Write();
           }
           t1->Write();

           }
           */
    }

    /*
    bool Smooth(vector<float> A,int pos,float th){
        TH1F *h1 = new TH1F("h1","",20,1,20);
        for(int i=1;i<20;i++)
        {
            h1->SetBinContent(i,-1*A.at(pos-10+i));
        }
        h1->Smooth();
        h1->FindFirstBinAbove(-1*th);

    }*/
