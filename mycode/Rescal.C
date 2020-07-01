#include <TH1.h>
#include <TF1.h>
#include <stdio.h>

char buff[1024];
int counter = 0;

float cal_parallelRes(float R1=50, float R2=1.1)
{
    float val=0;
    val= (R1*R2)/(R1+R2);
    //cout<<"R1= "<<R1<<"\t"<<",R2= "<<R2<<endl;
    //cout<<"R1 || R2 = "<<val<<endl;
    return val;

}
float val_parallelRes(float R1=50, float Rpara=1.1)
{
    float R2=0;
    R2=TMath::Abs(Rpara*R1/(R1-Rpara));
    //cout<<"R1= "<<R1<<"\t"<<",R2= "<<R2<<endl;
    //cout<<"R1 || R2 = "<<Rpara<<endl;
    return R2;
}


void baseRes(double HV1=300, double HV2=800, double HV3=800, double HV4=200,double I=100,double MCP1 = 200, double MCP2=200)
{
    //MCP res unit: Mohm
    // HV unit is V
    //I unit: uA
    double HV = HV1+HV2+HV3+HV4;
    double R1,R2,R3,R4;
    R1 = HV1/I;
    R2 = HV2/I;
    R3 = HV3/I;
    R4 = HV4/I;
    if(MCP1>0)
    R2 = val_parallelRes(R2,MCP1);
    if(MCP2>0)
    R3 = val_parallelRes(R3,MCP2);
    cout<<" [ *Usage* : Scanf the value of divided HV you set. ] "<<endl;
    cout<<">> HV set to: "<<HV<<"V"<<endl;
    cout<<">> I set to: "<<I<<"uA"<<endl;
    cout<<"-->> The resistor value is: "<<R1<<"M, "<<R2<<"M, "<<R3<<"M, "<<R4<<"M"<<endl;
}
void baseHV(double R1=3, double R2=8.333, double R3=8.333, double R4=2,double I=100,double MCP1 = 200, double MCP2=200)
{
    //MCP res unit: Mohm
    // HV unit is V
    //I unit: uA
    double HV1,HV2,HV3,HV4;
    if(MCP1>0)
    R2 = cal_parallelRes(R2,MCP1);
    if(MCP2>0)
    R3 = cal_parallelRes(R3,MCP2);
    
    HV1 = R1*I;
    HV2 = R2*I;
    HV3 = R3*I;
    HV4 = R4*I;
    double HV = HV1+HV2+HV3+HV4;
    cout<<" [ *Usage* : Scanf the value of resistors on base. ] "<<endl;
    cout<<">> The Res value actually is: "<<R1<<"M, "<<R2<<"M, "<<R3<<"M, "<<R4<<"M"<<endl;
    cout<<">> I set to: "<<I<<"uA"<<endl;
    cout<<"-->> Divided HV is: "<<HV1<<"V, "<<HV2<<"V, "<<HV3<<"V, "<<HV4<<"V"<<endl;
    cout<<"-->> HV set to: "<<HV<<"V"<<endl;
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