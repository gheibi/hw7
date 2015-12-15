#include<iostream>
#include<cmath>
using namespace std;

void Arenstrof(double*p,const double mass,double x,double y,double xp,double yp);
void maximum (double *norm,double & MAX);

int main(){

  const double a21=1.0/5.0;
  const double a31=3.0/40.0,       a32=9.0/40.0;
  const double a41=44.0/45.0,      a42=-56.0/15.0,      a43=32/9;
  const double a51=19372.0/6561.0, a52=-25360.0/2187.0, a53=64448.0/6561.0, a54=-212.0/729.0;
  const double a61=9017.0/3168.0,  a62=-355.0/33.0,     a63=46732.0/5247.0, a64=49.0/176.0,  a65=-5103.0/18656.0;
  const double a71=35.0/384.0,     a72=0.0,             a73=500.0/1113.0,   a74=125.0/192.0, a75=-2187.0/6784.0,    a76=11.0/84.0;
//-------------------------------------------------------------------------------------------------------------------------  
  const double b14=35.0/384.0,    b24=0.0,            b34=500.0/1113.0,  b44=125.0/192.0,b54=-2187.0/6784.0,  b64=11.0/84.0;
  const double b15=5179.0/51600.0,b25=0.0,            b35=7571.0/16695.0,b45=393.0/640.0,b55=-92097.0/339200.0,b65=187.0/2100.0, b75=1.0/40.0;
  
  double k1[4], k2[4], k3[4], k4[4], k5[4], k6[4],k7[4], RK4_var[4], RK5_var[4], norm[4];

  double x=0.994;                    // first variation x
  double y=0;                       // second variation y 
  double xp=0;                     // third variation   x'=xp
  double yp=-2.00158510637908;    // fourth variation   y'=yp
  double const mass=0.012277471; // mass = mass ratio
  double const  Tol=1e-5;
  double dt=1e-4;
  double MAX=0.0;
  double const  T=17.065216560157;
  double q=0.1;
     RK4_var[0]=x;
     RK4_var[1]=y;
     RK4_var[2]=xp;
     RK4_var[3]=yp;
                   RK5_var[0]=x;
                   RK5_var[1]=y;
                   RK5_var[2]=xp;
                   RK5_var[3]=yp;
   cout<<0.0<<'\t'<<dt<<'\t'<<x<<'\t'<<y<<'\t'<<xp<<'\t'<<yp<<endl;
for(double t=dt; t<=T; t+=dt){
  
 Arenstrof (k1,mass, x, y, xp, yp);
 Arenstrof (k2,mass, x+dt*(a21*k1[0]),   y+dt*(a21*k1[1]),  xp+dt*(a21*k1[2]), yp+dt*(a21*k1[3]));
 Arenstrof (k3,mass, x+dt*(a31*k1[0]+a32*k2[0]),   y+dt*(a31*k1[1]+a32*k2[1]),  xp+dt*(a31*k1[2]+a32*k2[2]), yp+dt*(a31*k1[3]+a32*k2[3]));
 Arenstrof (k4,mass, x+dt*(a41*k1[0]+a42*k2[0]+a43*k3[0]),    y+dt*(a41*k1[1]+a42*k2[1]+a43*k3[1]),   xp+dt*(a41*k1[2]+a42*k2[2]+a43*k3[2]), yp+dt*(a41*k1[3]+a42*k2[3]+a43*k3[3]));
 Arenstrof (k5,mass, x+dt*(a51*k1[0]+a52*k2[0]+a53*k3[0]+a54*k4[0]),    y+dt*(a51*k1[1]+a52*k2[1]+a53*k3[1]+a54*k4[1]),   xp+dt*(a51*k1[2]+a52*k2[2]+a53*k3[2]+a54*k4[2]), yp+dt*(a51*k1[3]+a52*k2[3]+a53*k3[3]+a54*k4[3]));
 Arenstrof (k6,mass, x+dt*(a61*k1[0]+a62*k2[0]+a63*k3[0]+a64*k4[0]+a65*k5[0]),    y+dt*(a61*k1[1]+a62*k2[1]+a63*k3[1]+a64*k4[1]+a65*k5[1]),   xp+dt*(a61*k1[2]+a62*k2[2]+a63*k3[2]+a64*k4[2]+a65*k5[2]), yp+dt*(a61*k1[3]+a62*k2[3]+a63*k3[3]+a64*k4[3]+a65*k5[3]));
 Arenstrof (k7,mass, x+dt*(a71*k1[0]+a72*k2[0]+a73*k3[0]+a74*k4[0]+a75*k5[0]+a76*k6[0]),    y+dt*(a71*k1[1]+a72*k2[1]+a73*k3[1]+a74*k4[1]+a75*k5[1]+a76*k6[1]),   xp+dt*(a71*k1[2]+a72*k2[2]+a73*k3[2]+a74*k4[2]+a75*k5[2]+a76*k6[2]), yp+dt*(a71*k1[3]+a72*k2[3]+a73*k3[3]+a74*k4[3]+a75*k5[3]+a76*k6[3]));

   for(int i=0; i<=3; i++){
    RK4_var[i]+=dt*(b14*k1[i]+b24*k2[i]+b34*k3[i]+b44*k4[i]+b54*k5[i]+b64*k6[i]);
    RK5_var[i]+=dt*(b15*k1[i]+b25*k2[i]+b35*k3[i]+b45*k4[i]+b55*k5[i]+b65*k6[i]+b75*k7[i]);
 
    norm[i]=abs(RK5_var[i]-RK4_var[i]);
   }
    maximum (norm, MAX);
     x=RK4_var[0];
     y=RK4_var[1];
     xp=RK4_var[2];
     yp=RK4_var[3];

  cout<<t<<'\t'<<dt<<'\t'<<x<<'\t'<<y<<'\t'<<xp<<'\t'<<yp<<endl;
   
   dt*=q*pow((Tol/MAX),0.2);
     MAX=0.0; 
}
  return 0;
}

void Arenstrof(double*p ,const double mass,double x,double y,double xp,double yp){
double r=sqrt(pow(x+mass,2)+pow(y,2));
double s=sqrt(pow(x-1.0+mass,2)+pow(y,2));
p[0]=xp;     
p[1]=yp;    
p[2]=x+2.0*yp-((1.0-mass)*(x+mass)/pow(r,3))-(mass/pow(s,3))*(x-1.0+mass);
p[3]=y-2.0*xp-((1.0-mass)*y/pow(r,3))-(mass*y/pow(s,3));
}

void maximum(double *norm,double & MAX){
 for(int i=0; i<=3; i++)
   if (norm[i]>MAX) MAX=norm[i];
}
