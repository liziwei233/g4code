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
float val_parallelRes(float R0=50, float val=1.1)
{
    float Rx=0;
    Rx=val*R0/(R0-val);
    cout<<"R0= "<<R0<<"\t"<<",Rx= "<<Rx<<endl;
    cout<<"R0 || Rx = "<<val<<endl;
    return Rx;
}


double TraceWidth(double TraceToPlaneDistance = 1.)
{
    // TraceToPlaneDistance - unit:mm
    double CharacteristicImpedance=50;
    double Er=4.05;
    double TraceHeight= 0.035;// 1OZ=35um or 1.4mil
    double TraceWidth=0;
    
    TraceWidth = ((5.98*TraceToPlaneDistance)/exp(CharacteristicImpedance/(60/sqrt(Er*(1-exp(-1.55*(0.00002+TraceToPlaneDistance)/TraceToPlaneDistance)))))-TraceHeight)/0.8;
    cout<<TraceWidth<<endl;
    return TraceWidth;

}
double CharacteristicImpedance(double TraceWidth=2.4871,double TraceToPlaneDistance = 1.5)
{
    // TraceWidth, TraceToPlaneDistance, unit is mm
    double CharacteristicImpedance = 50;
    double Er=4.05;
    double TraceHeight= 0.035;// 1OZ=35um or 1.4mil
    double Impendance = 0;
    Impendance = (60/sqrt(Er*(1-exp(-1.55*(0.00002+TraceToPlaneDistance)/TraceToPlaneDistance))))*log(5.98*TraceToPlaneDistance/(0.8*TraceWidth+TraceHeight));
    cout<<Impendance<<endl;
    return Impendance;
}