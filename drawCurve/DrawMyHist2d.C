TH2F* DrawMyHist2d(const char* name="hist",double x1=0,double x2=1,double y1=0,double y2=1,const char *xtitle="",const char *ytitle="", Color_t LColor=1, Width_t LWidth=1.5)
{
	TH2F *datahist= new TH2F(name,"",200,x1,x2,200,y1,y2);
	TCanvas *c1 =new TCanvas("c1","c1",800,600);
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
	 datahist->GetXaxis()->SetTitleFont( 42 );
     datahist->GetYaxis()->SetTitleFont( 42 );
	 
	 Color_t TitleColor=1;
     datahist->GetXaxis()->SetTitleColor( TitleColor);
     datahist->GetYaxis()->SetTitleColor( TitleColor );
     datahist->GetXaxis()->SetTitleSize(0.06);
     datahist->GetYaxis()->SetTitleSize(0.06);
     datahist->GetXaxis()->SetTitleOffset(1.0);
     datahist->GetYaxis()->SetTitleOffset(1.0);
	 //datahist->GetXaxis()->SetBorderSize(5);
     datahist->GetXaxis()->SetNdivisions(505);
     datahist->GetYaxis()->SetNdivisions(510);
     datahist->GetXaxis()->CenterTitle();
     datahist->GetYaxis()->CenterTitle();
	 return datahist;
}
/*
void DrawMyHist2d(TH2 *datahist, char *xtitle, char *ytitle, Color_t LColor, Width_t LWidth){
	
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
     datahist->GetXaxis()->SetTitleColor( TitleColor);
     datahist->GetYaxis()->SetTitleColor( TitleColor );
     datahist->GetXaxis()->SetTitleSize(0.06);
     datahist->GetYaxis()->SetTitleSize(0.06);
     datahist->GetXaxis()->SetTitleOffset(1.0);
     datahist->GetYaxis()->SetTitleOffset(1.0);
     //datahist->GetXaxis()->SetBorderSize(5);
     datahist->GetXaxis()->SetNdivisions(510);
     datahist->GetYaxis()->SetNdivisions(510);
     datahist->GetXaxis()->CenterTitle();
     datahist->GetYaxis()->CenterTitle();
}*/
