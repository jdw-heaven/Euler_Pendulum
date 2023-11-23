#ifndef __MY_EULER_H__
#define __MY_EULER_H__
/*The above is to prevent the header file from being included multiple times,
it can be omitted, it is best to have, any name, and ensure that it is unique*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
/*The following macro definition is optional, here we define the parameter we need*/
#define PI 3.14159265358979
#define E 1e-3
#define Bg 1.20 //we resume that the neodymium magnet's remanence is the averger value of 1.18-1.22T which is the N35's parameter
//#define mu0 4*PI*0.0000001 //the unit is T*m/A  to make double available, we chose to use Bg/mu0 which we called Bm.
#define Bm 3e4/PI
#define g 9.80 //N/kg I got it through experiment
#define R 1.45 //cm 
#define h 1.946 //cm 
#define r 0.475 //cm
#define l 7.6 //cm
#define d 0.387 //cm
#define m 0.0443 //kg

/*The following is the definition of functions*/
/*the functions to calculate the mgnetic field*/
//the component of x axis
double m_bx(double x, double y, double z, double w, double th){
    double K = pow( (x-R*cos(th))*(x-R*cos(th))+(y-R*sin(th))*(y-R*sin(th))+(z-w)*(z-w) , 1.5);
    return (z-w)*R*cos(th)/K;
}
double m_bx_th(double x, double y, double z, double w){
	int i=1,j,k,n=5;
	double T[5],a=0.1,b=2*PI,s=0.0; 
	T[0]=0.5*(b-a)*(m_bx(x,y,z,w,a)+m_bx(x,y,z,w,b)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_bx(x,y,z,w,a+(2*k-1)*(b-a)/pow(2,j)); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[4]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[4]-T[0])>E;i++) 
	{ 
		T[0]=T[4]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[4]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
	return T[4];
}
double m_Bx(double x, double y, double z){
    int i=1,j,k,n=5;
	double T[5],a=0.0,b=h,s=0.0; 
	T[0]=0.5*(b-a)*(m_bx_th(x,y,z,a)+m_bx_th(x,y,z,b)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_bx_th(x,y,z,a+(2*k-1)*(b-a)/pow(2,j)); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[4]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[4]-T[0])>E;i++) 
	{ 
		T[0]=T[4]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[4]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
	return T[4];
}
//the component of y axis
double m_by(double x, double y, double z, double w, double th){
    double K = pow( (x-R*cos(th))*(x-R*cos(th))+(y-R*sin(th))*(y-R*sin(th))+(z-w)*(z-w) , 1.5);
    //printf("%10.4lf%10.4lf%10.4lf%10.4lf%10.4lf%10.4lf\n", x, y, z, w, th, (z-w)*R*sin(th)/K);
    return (z-w)*R*sin(th)/K;
}
double m_by_th(double x, double y, double z, double w){
	int i=1,j,k,n=5;
	double T[5],a=0.1,b=2*PI,s=0.0; 
	T[0]=0.5*(b-a)*(m_by(x,y,z,w,a)+m_by(x,y,z,w,b)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_by(x,y,z,w,a+(2*k-1)*(b-a)/pow(2,j)); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[4]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[4]-T[0])>E;i++) 
	{ 
		T[0]=T[4]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[4]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
	return T[4];
}
double m_By(double x, double y, double z){
    int i=1,j,k,n=5;
	double T[5],a=0.0,b=h,s=0.0; 
	T[0]=0.5*(b-a)*(m_by_th(x,y,z,a)+m_by_th(x,y,z,b)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_by_th(x,y,z,a+(2*k-1)*(b-a)/pow(2,j)); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[4]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[4]-T[0])>E;i++) 
	{ 
		T[0]=T[4]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[4]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
	return T[4];
}
//the component of z axis
double m_bz(double x, double y, double z, double w, double th){
    double K = pow( (x-R*cos(th))*(x-R*cos(th))+(y-R*sin(th))*(y-R*sin(th))+(z-w)*(z-w) , 1.5);
    return ( -(y-R*sin(th))*R*sin(th)-(x-R*cos(th))*R*cos(th) )/K;
}
double m_bz_th(double x, double y, double z, double w){
	int i=1,j,k,n=5;
	double T[5],a=0.1,b=2*PI,s=0.0; 
	T[0]=0.5*(b-a)*(m_bz(x,y,z,w,a)+m_bz(x,y,z,w,b)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_bz(x,y,z,w,a+(2*k-1)*(b-a)/pow(2,j)); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[4]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[4]-T[0])>E;i++) 
	{ 
		T[0]=T[4]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[4]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
	return T[4];
}
double m_Bz(double x, double y, double z){
    int i=1,j,k,n=5;
	double T[5],a=0.0,b=h,s=0.0; 
	T[0]=0.5*(b-a)*(m_bz_th(x,y,z,a)+m_bz_th(x,y,z,b)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_bz_th(x,y,z,a+(2*k-1)*(b-a)/pow(2,j)); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[4]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[4]-T[0])>E;i++) 
	{ 
		T[0]=T[4]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[4]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
	return T[4];
}
//fx
double m_fx(double b, double u, double al){
    return ( -r*sin(b)*cos(al)*m_Bz(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)) \
    + r*sin(b)*sin(al)*m_By(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)));
}
double m_fx_b(double u, double al){
	int i=1,j,k,n=5;
	double T[5],a=0.0,b=2*PI,s=0.0; 
	T[0]=0.5*(b-a)*(m_fx(a,u,al)+m_fx(b,u,al)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_fx(a+(2*k-1)*(b-a)/pow(2,j),u,al); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[4]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[4]-T[0])>E/100;i++) 
	{ 
		T[0]=T[4]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[4]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
}
double m_Fx(double al){
	int i=1,j,k,n=5;
	double T[5],a=0.0,b=l,s=0.0; 
	T[0]=0.5*(b-a)*(m_fx_b(a,al)+m_fx_b(b,al)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_fx_b(a+(2*k-1)*(b-a)/pow(2,j),al); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[4]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[4]-T[0])>E/100;i++) 
	{ 
		T[0]=T[4]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[4]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
	return Bm*T[4];
}
//fy
double m_fy(double b, double u, double al){
    return ( -r*sin(b)*sin(al)*m_Bx(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)) \
    + r*cos(b)*m_Bz(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)));
}
double m_fy_b(double u, double al){
	int i=1,j,k,n=5;
	double T[5],a=0.0,b=2*PI,s=0.0; 
	T[0]=0.5*(b-a)*(m_fy(a,u,al)+m_fy(b,u,al)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_fy(a+(2*k-1)*(b-a)/pow(2,j),u,al); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[4]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[4]-T[0])>E/100;i++) 
	{ 
		T[0]=T[4]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[4]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
	return T[4];
}
double m_Fy(double al){
	int i=1,j,k,n=5;
	double T[5],a=0.0,b=l,s=0.0; 
	T[0]=0.5*(b-a)*(m_fy_b(a,al)+m_fy_b(b,al)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_fy_b(a+(2*k-1)*(b-a)/pow(2,j),al); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[4]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[4]-T[0])>E/100;i++) 
	{ 
		T[0]=T[4]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[4]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
	return Bm*T[4];
}
//fz
double m_fz(double b, double u, double al){
    return ( -r*cos(b)*m_By(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)) \
    + r*sin(b)*cos(al)*m_Bx(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)));
}
double m_fz_b(double u, double al){
	int i=1,j,k,n=5;
	double T[5],a=0.0,b=2*PI,s=0.0; 
	T[0]=0.5*(b-a)*(m_fz(a,u,al)+m_fz(b,u,al)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_fz(a+(2*k-1)*(b-a)/pow(2,j),u,al); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[4]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[4]-T[0])>E/100;i++) 
	{ 
		T[0]=T[4]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[4]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
	return T[4];
}
double m_Fz(double al){
	int i=1,j,k,n=5;
	double T[5],a=0.0,b=l,s=0.0; 
	T[0]=0.5*(b-a)*(m_fz_b(a,al)+m_fz_b(b,al)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_fz_b(a+(2*k-1)*(b-a)/pow(2,j),al); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[4]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[4]-T[0])>E/100;i++) 
	{ 
		T[0]=T[4]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[4]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
	return Bm*T[4];
}
#endif