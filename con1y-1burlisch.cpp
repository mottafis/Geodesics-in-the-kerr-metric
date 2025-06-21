#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
using namespace std;
double a = 0.9;//se definen los parametros para describir el agujero negro y la particula
double M = 1.0;//masa agujero negro, todo esta normalizado
double E = 0.9;
double C = 4.0;
//se definen funciones para despues definir la evolucion conrespecto al parametro tao de las coordenadas, estas coordenadas antes de convertirlas en el codigo de las graficas son las de Boyer–Lindquist, con estas se define la metrica
double delta(double r) {
return r*r-2*M*r+a*a;
}
double rhocua(double r, double theta) {
return r*r+a*a*cos(theta)*cos(theta);
}
double sigma2(double r, double theta) {
return (r*r+a*a)*(r*r+a*a)-a*a*delta(r)*sin(theta)*sin(theta);
}
//ecuaciones de evolucion segun chandrasekar para un foton
double rpun(double r, double theta, double Lz) {
double numerador= ((r*r+a*a)*E-a*Lz)*((r*r+a*a)*E-a*Lz)-delta(r)*C;
double denominador= rhocua(r,theta)*rhocua(r,theta);
return -sqrt(numerador/denominador);
}
double thetapun(double r, double theta, double Lz) {
double numerador= C-(a*E*sin(theta)-Lz*(1/sin(theta)))*(a*E*sin(theta)-Lz*(1/sin(theta)));
double denominador= rhocua(r,theta)*rhocua(r,theta);
return sqrt(numerador/denominador);
}
double phipun(double r, double theta, double Lz) {
double numerador= (2*a*M*r*E+(rhocua(r,theta)-2*M*r)*Lz*(1/sin(theta))*(1/sin(theta)));
double denominador= (delta(r)*rhocua(r,theta));
return numerador/denominador;
}
double tpun(double r, double theta, double Lz) {
double numerador=(sigma2(r, theta)*E-2*a*M*r*Lz);
double denominador= delta(r);
return numerador/denominador;
}
//se define una estructura para facilitar el manejo de cada coordenada
struct estado {
double r;
double theta;
double t;
double phi;
};
//se define una funcion con las derivadas
estado f(double tao, const estado &y, double Lz) {
estado derivadas;
derivadas.r = rpun(y.r, y.theta, Lz);
derivadas.theta = thetapun(y.r, y.theta, Lz);
derivadas.t = tpun(y.r, y.theta, Lz);
derivadas.phi = phipun(y.r, y.theta, Lz);
return derivadas;
}
//se ejecuta en una funcion el metodo del punto medio segun las ecuaciones de...
estado puntomediomod(const estado &y0, double tao, double H, int n, double Lz) {
double h= H/n;
vector<estado> z(n+1);
z[0]= y0;
z[1]= z[0];
estado k= f(tao, z[0], Lz);
z[1].r += h*k.r;
z[1].theta += h*k.theta;
z[1].t += h*k.t;
z[1].phi += h*k.phi;
for (int i=1; i<n;i++) {
z[i+1].r= z[i-1].r+2*h*f(tao+i*h,z[i],Lz).r;
z[i+1].theta= z[i-1].theta + 2 * h * f(tao + i * h, z[i], Lz).theta;
z[i+1].t= z[i-1].t + 2 * h * f(tao + i*h,z[i],Lz).t;
z[i+1].phi= z[i-1].phi+2*h*f(tao+i*h, z[i], Lz).phi;
}
k= f(tao + H, z[n], Lz);
estado y_fin;
y_fin.r= 0.5*(z[n].r+z[n-1].r+h*k.r);
y_fin.theta= 0.5*(z[n].theta+z[n-1].theta+h*k.theta);
y_fin.t= 0.5*(z[n].t + z[n-1].t + h*k.t);
y_fin.phi = 0.5*(z[n].phi+z[n-1].phi+h*k.phi);
return y_fin;
}
//se hace una funcion de extrapolacion de richardson
estado extrapolacion(const vector<estado> &estimaciones,int p){
vector<estado> temp=estimaciones;//el ciclo en j limita la estimacion hasta p 
for(int j=1;j<=p;j++){
  for(int i=0;i<=p-j;i++){
   double factor=pow(2,2*j);//se hace con cada coordenada
   temp[i].r=(factor*temp[i+1].r-temp[i].r)/(factor-1);//formula de extrapolacion
   temp[i].theta=(factor*temp[i+1].theta-temp[i].theta)/(factor-1);
   temp[i].t=(factor*temp[i+1].t-temp[i].t)/(factor-1);
   temp[i].phi=(factor*temp[i+1].phi-temp[i].phi)/(factor-1);
  }
 }
return temp[0];//se imprime la de mayor precision porque se va reemplazando en cada iteracion
}
//se define una funcion que llame a las otras dos que son las base de Burlisch y se retorna al final la extrapolacion
estado burlisch(estado y0, double x, double H, int pmax, double Lz) {
vector<estado> estimaciones(pmax + 1);
for (int p = 1; p <= pmax; p++) {
int n = 2 * p;
estimaciones[p] = puntomediomod(y0, x, H, n, Lz);
}
return extrapolacion(estimaciones,pmax);
}
//se hace para que salgan los dos archivos de datos al mismo tiempo
void ejecutarSimulacion(double Lz, const string& nombreArchivo) {
estado s = {10.0, M_PI / 2 - 0.2, 0.0, 0.0};//se definen las condiciones iniciales
double h = 0.001;
int pasos = 10000;
ofstream datos(nombreArchivo);//se guardan los archivos separados por , y con los dos extremos, Lz=1 y Lz=-1
for (int i = 0; i < pasos; i++) {
datos << s.t << "," << s.r << "," << s.theta << "," << s.phi << endl;
s = burlisch(s, s.t, h, 8, Lz);
}
datos.close();
}
int main() {
//se ejecuta la ultima funcion que ejecuta las demas y guarda los datos en dos archivos dependiendo de Lz.
ejecutarSimulacion(1.0, "fotongeneralbur_1.dat");
ejecutarSimulacion(-1.0, "fotongeneralbur_-1.dat");
return 0;
}
