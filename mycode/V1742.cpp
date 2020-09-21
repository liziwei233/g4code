#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <stdio.h>
#include "Include/DrawMyClass.h"
const char path[1024]="/mnt/d/OneDrive - mail.ustc.edu.cn/Experiment/V1742code/binarytest";
//const char path[1024]="/mnt/d/OneDrive - mail.ustc.edu.cn/Experiment/V1742code/test";
void drawwaveform(float* x, float* y, int Npoint);
void readbinary(const char* filename);
void readtxt(const char* filename);
void Channel(char* fullpath);
int ReadBinary(vector<float> * wave);
bool GetNextEvent();
//char buff[1024];
const int Nsample=1024;
const float Tstep = 0.200; //unit: ns;
//float waveform_y[Nsample];
//float waveform_x[Nsample];

ifstream infile;
vector<float> wave;
vector<float> waveform_x;
vector<float> waveform_y;
float amp;
int event_number;

int main(){
    char filename[1024]="wave_1.dat";
    sprintf(buff,"%s/%s",path,filename);
    cout << "The file path: " ;
    cout << buff << endl;
    Channel(buff);
    while(!GetNextEvent())
    {
        
        cout<<"Event: "<<event_number<<endl;
    }
    //readtxt(buff);
    return 0;
}
void Channel(char* fullpath){
    event_number = 0;
    infile.open(fullpath,ifstream::binary);
    if(!infile){
        cout<<"File Error!"<<endl;
        return;
    }
}
int ReadBinary(vector<float> * wave)
{
     
    
    amp=0;

    if(infile.eof()){
        return 1;
    }
    
    for(int i=0;i<1024;i++)
        {
            infile.read((char*)&amp,sizeof(float));
            wave->push_back(amp);
            //cout<<wave->at(i)<<endl;
        }

        return 0;
}
bool GetNextEvent()
{
    wave.clear();
    waveform_x.clear();
    waveform_y.clear();
    if(ReadBinary(&wave))
        {cerr << "Unknown format " << endl << endl;
        infile.close();
        return 1;}

    for(int i = 0; i < 1024; ++i)
    {
        waveform_x.push_back(i*0.2);
        waveform_y.push_back(wave.at(i));
        //cout<<"waveform_x= "<<waveform_x.at(i)<<", waveform_y= "<<waveform_y.at(i)<<endl;
    }

    event_number++;
    //cout<<" reading Event "<<event_number<<"..."<<endl;
    return 0;
}

/*
void readbinary(const char* filename){
    
    sprintf(buff,"%s.dat",filename);
    cout<<"Reading File:  "<<buff<<endl;
    ifstream infile(buff,ios::binary);
    
    if(!infile){
        cout<<"File error!"<<endl;
        return;
    }
    float amp;
    int counter=0;
    while(infile){
        for(int i=0;i<Nsample;i++)
        {
            infile.read((char*)&amp,sizeof(float));
            waveform_y[i]=amp;
            waveform_x[i]=i*Tstep;
        //cout<<"waveform_x= "<<waveform_x[i]<<", waveform_y= "<<waveform_y[i]<<endl;
        }
        //drawwaveform(waveform_x,waveform_y,Nsample);
        counter++;
        cout<<"the number of waveforms: "<<counter<<endl;
        //break;
    }
    infile.close();
}

void readtxt(const char* filename){
    
    sprintf(buff,"%s.txt",filename);
    cout<<"Reading File:  "<<buff<<endl;
    ifstream infile(buff,ios::in);
    
    if(!infile){
        cout<<"File error!"<<endl;
        return;
    }
    double amp;
    int counter=0;
    while(infile){
        for(int i=0;i<Nsample;i++)
        {
            infile>>amp;
            //infile.read((char*)&amp,sizeof(double));
            waveform_y[i]=amp;
            waveform_x[i]=i*Tstep;
        cout<<"waveform_x= "<<waveform_x[i]<<", waveform_y= "<<waveform_y[i]<<endl;
        }
        //drawwaveform(waveform_x,waveform_y,Nsample);
        counter++;
        cout<<"the number of waveforms: "<<counter<<endl;
        break;
    }
    infile.close();
}

void drawwaveform(float* x, float* y, int Npoint)
{
    TCanvas *c = cdC(0);
    TGraph *g = new TGraph(Npoint,x,y);
    g->Draw("AP");
}
*/