#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
using namespace std;
//se definen los parametros con los que va a trabajar el resto del codigo
double a=0.9;
double M=1.0;
double E=0.9;
double Lz=-1.0;
double C=4.0;
//se definen funciones que se utilizaran dentro de otras en el avance de cada coordenada
double delta(double r){
 return r*r-2*M*r+a*a;
}
double rhocua(double r, double theta){
 return r*r+a*a*cos(theta)*cos(theta);
}
double sigma2(double r, double theta){
 return (r*r+a*a)*(r*r+a*a)-a*a*delta(r)*sin(theta)*sin(theta);
}
//se definen funciones de evolucion en cada coordenada
double rpun(double r,double theta){
 double numerador=((r*r+a*a)*E-a*Lz)*((r*r+a*a)*E-a*Lz)-delta(r)*C;
 double denominador=rhocua(r,theta)*rhocua(r,theta);
// if(denominador==0||numerador<0){
// return 0;//por errores por sqrt y 0 en denom, pero ya se vio donde se dan estos errores, de 1 a -1 esta bien en Lz
// }
return -sqrt(numerador/denominador);
}
double thetapun(double r,double theta){
 double numerador=C-(a*E*sin(theta)-Lz*(1/sin(theta)))*(a*E*sin(theta)-Lz*(1/sin(theta)));
 double denominador=rhocua(r,theta)*rhocua(r,theta);
// if(denominador==0||numerador<0){
// return 0;//por errores por sqrt y 0 en denom
// }
return sqrt(numerador/denominador);
}
double phipun(double r, double theta){
 double numerador=(2*a*M*r*E+(rhocua(r,theta)-2*M*r)*Lz*(1/sin(theta))*(1/sin(theta)));
 double denominador=(delta(r)*rhocua(r,theta));
// if(denominador==0||numerador<0){
// return 0;//por errores por sqrt y 0 en denom
// }
return numerador/denominador;
}
double tpun(double r, double theta){
 double numerador=(sigma2(r,theta)*E-2*a*M*r*Lz);
 double denominador=delta(r);
// if(denominador==0){
// return 0;
// }
return numerador/denominador;
}
//se define una estructura que se va a usar para guardar y manejar datos de manera mas facil
struct estado{
double r;
double theta;
double t;
double phi;
};
//se define una funcion con el metodo de runge kutta de 4to orden
estado kuta4(estado s,double h){
double k1r=h*rpun(s.r,s.theta);
double k1theta=h*thetapun(s.r,s.theta);
double k1t=h*tpun(s.r,s.theta);
double k1phi=h*phipun(s.r,s.theta);
double k2r=h*rpun(s.r+0.5*k1r,s.theta+0.5*k1theta);
double k2theta=h*thetapun(s.r+0.5*k1r,s.theta+0.5*k1theta);
double k2t=h*tpun(s.r+0.5*k1r,s.theta+0.5*k1theta);
double k2phi=h*phipun(s.r+0.5*k1r,s.theta+0.5*k1theta);
double k3r=h*rpun(s.r+0.5*k2r,s.theta+0.5*k2theta);
double k3theta=h*thetapun(s.r+0.5*k2r,s.theta+0.5*k2theta);
double k3t=h*tpun(s.r+0.5*k2r,s.theta+0.5*k2theta);
double k3phi=h*phipun(s.r+0.5*k2r,s.theta+0.5*k2theta);
double k4r=h*rpun(s.r+k3r,s.theta+k3theta);
double k4theta =h*thetapun(s.r+k3r,s.theta+k3theta);
double k4t=h*tpun(s.r+k3r,s.theta+k3theta);
double k4phi=h*phipun(s.r+k3r,s.theta+k3theta);
estado s_next;
  s_next.r=s.r+(k1r+2*k2r+2*k3r+k4r)/6.0;
  s_next.theta=s.theta+(k1theta+2*k2theta+2*k3theta+k4theta)/6.0;
  s_next.t=s.t+(k1t+2*k2t+2*k3t+k4t)/6.0;
  s_next.phi=s.phi+(k1phi+2*k2phi+2*k3phi+k4phi)/6.0;
return s_next;
}
int main(){
double h=0.001;
int pasos=10000;
estado s={10,M_PI/2-0.2,0.0,0.0};//condiciones iniciales en cada coordenada
ofstream datos("fotongeneral.dat");//se guardan los datos
for(int i=0;i<pasos;i++){
 datos<<s.t<<","<<s.r<<","<<s.theta<<","<<s.phi<<endl;
 s=kuta4(s,h);
}
datos.close();
return 0;
}
