#include <TH1.h>
#include <TF1.h>
#include <stdio.h>

char buff[1024];
int counter = 0;

float cal_parallelRes(float R1=50, float R2=1.1)
{
    float val=0;
    val= (R1*R2)/(R1+R2);
    cout<<"R1= "<<R1<<"\t"<<",R2= "<<R2<<endl;
    cout<<"R1 || R2 = "<<val<<endl;
    return val;
}