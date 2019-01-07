void txt2wave(){
	TFile* f1 = new TFile("savewave.root","RECREATE");
	TTree* t1 = new TTree("t1","");
	float ta1,A1;
	float ta2,A2;
	float tb1,B1;
	char comma;
	char sentence[1024];
	char histname[1024];
	char buff[1024];
	
	TCanvas* c1=new TCanvas("c1","c1",800,600);
	TCanvas* c2=new TCanvas("c2","c2",800,600);
	TCanvas* c3=new TCanvas("c3","c3",800,600);
	
	int counter=0;
	t1->Branch("ta1",&ta1,"ta1/F");
	t1->Branch("ta2",&ta2,"ta2/F");
	t1->Branch("tb1",&tb1,"tb1/F");
	t1->Branch("A1",&A1,"A1/F");
	t1->Branch("A2",&A2,"A2/F");
	t1->Branch("B1",&B1,"B1/F");
	
	TH1F *hA1[200];
	TH1F *hA2[200];
	TH1F *hB1[200];
	
	ifstream input1,input2,input3;
	for(int k=0;k<2;k++){
		sprintf(buff,"C2--A1A2B1B2--0000%d.txt",k);
		input1.open(buff);
		sprintf(buff,"C3--A1A2B1B2--0000%d.txt",k);
		input2.open(buff);
		sprintf(buff,"C4--A1A2B1B2--0000%d.txt",k);
		input3.open(buff);
		
	if(!input1){
		cout<<buff<<" file can't be found!"<<endl;
		return;
	}	
	//while(input1&&(!input1.eof())&&counter<3){
	while(input1&&(!input1.eof())){
			sprintf(histname,"hmcp_%d",counter);
			hA1[counter]= new TH1F(histname,";T (25 ps);Amp (V)",20002,0,20002);//25ps/point 
			sprintf(histname,"hB1_%d",counter);
			hA2[counter]= new TH1F(histname,";T (25 ps);Amp (V)",20002,0,20002);//25ps/point
			sprintf(histname,"hB2_%d",counter);
			hB1[counter]= new TH1F(histname,";T (25 ps);Amp (V)",20002,0,20002);//25ps/point
			
			//hA1[counter]->Reset();
			input1>>sentence;
			input2>>sentence;
			input3>>sentence;
			for(int i=0;i<20002;i++){
			input1>>ta1>>comma>>A1;
			input2>>ta2>>comma>>A2;
			input3>>tb1>>comma>>B1;
			hA1[counter]->SetBinContent(i,A1);
			hA2[counter]->SetBinContent(i,A2);
			hB1[counter]->SetBinContent(i,B1);
			t1->Fill();
			}
			
		
		t1->Fill();
		c1->cd();
		
		hA1[counter]->SetLineColor(2);
		//hA1[counter]->Draw();
		
		hA2[counter]->SetLineColor(6);
		//hA2[counter]->Draw("same");
		
		hB1[counter]->SetLineColor(9);
		//hB1[counter]->Draw("same");
		cout<<"The Wave ID is:"<<counter<<endl;
		counter++;
		
	}
	input1.close();
	input2.close();
	input3.close();
	}
	
	
	
	
	if(f1->IsOpen()){
		for(int s=0;s<counter;s++){
			hA1[s]->Write();
			hA2[s]->Write();
			hB1[s]->Write();
		}
		t1->Write();
		
	}

}