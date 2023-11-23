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
//The inner layer of the double integral
double m_Bx_th_add(double x, double y, double z, double w, int k, double th1, double t){
    double a = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        a += m_bx(x, y, z, w, th1+(2*i-1)*t/pow(2,k-1));
    }
    return a*t/(k>2? pow(2,k-2):1);
}
double m_Bx_th(double x, double y, double z, double w){
    double answer = 0;
    double th1 = 3*E, th2 = 2*PI;
    double t = th2-th1;
    int M = 40;     //the initial M, if it's not enough, we'll give you an alert and you need to enlarge M
    double *Px;
    Px = (double *)malloc(M*M*sizeof(double));

    Px[0] = ( m_bx(x, y, z, w,th1)+m_bx(x, y, z, w, th2) )*t/2;
    int i = 1, j = 1;
    while(i){
        Px[i*M+0] = (Px[(i-1)*M+0]+m_Bx_th_add(x,y,z,w,i+1,th1,t))/2;
        while(j){
            Px[i*M+j] = Px[i*M+j-1]+(Px[i*M+j-1]-Px[(i-1)*M+j-1])/(pow(4,j)-1);
            if(j == i){
                j = 1;
                break;
            }else{
                j++;
            }
        }
        double e = fabs(Px[i*M+i]-Px[(i-1)*M+i-1]);
        if(e<E&&i!=1){
            answer = Px[i*M+i];
            break;
        }else if(i == M-1){
            printf("x=%lf, y=%lf, z=%lf\n", x, y, z);
            printf("Please enlarge Px!!!\n");
        }else{
            i++;
        }
    }

    free(Px);
    return answer;
}
//the magnet's x component
double m_Bx_add(double x, double y, double z, int k, double w1, double t){
    double a = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        a += m_Bx_th(x, y, z, w1+(2*i-1)*t/pow(2,k-1));
    }
    return a*t/(k>2? pow(2,k-2):1);
}
double m_Bx(double x, double y, double z){
    double answer = 0;
    double w1 = 0, w2 = h;
    double t = w2-w1;
    int M = 40;
    double *Hx;
    Hx = (double *)malloc(M*M*sizeof(double));

    Hx[0] = ( m_Bx_th(x, y, z, w1)+m_Bx_th(x, y, z, w2) )*t/2;
    int i = 1, j = 1;
    while(i){
        Hx[i*M+0] = (Hx[(i-1)*M+0]+m_Bx_add(x,y,z,i+1,w1,t))/2;
        while(j){
            Hx[i*M+j] = Hx[i*M+j-1]+(Hx[i*M+j-1]-Hx[(i-1)*M+j-1])/(pow(4,j)-1);
            if(j == i){
                j = 1;
                break;
            }else{
                j++;
            }
        }
        double e = fabs(Hx[i*M+i]-Hx[(i-1)*M+i-1]);
        if(e<E&&i!=1){
            answer = Bg*Hx[i*M+i]/(4*PI);
            break;
        }else if(i == M-1){
            printf("x=%lf, y=%lf, z=%lf\n", x, y, z);
            printf("Please enlarge Hx!!!\n");
        }else{
            i++;
        }
    }

    free(Hx);
    return answer;
}
//the component of y axis
double m_by(double x, double y, double z, double w, double th){
    double K = pow( (x-R*cos(th))*(x-R*cos(th))+(y-R*sin(th))*(y-R*sin(th))+(z-w)*(z-w) , 1.5);
    //printf("%10.4lf%10.4lf%10.4lf%10.4lf%10.4lf%10.4lf\n", x, y, z, w, th, (z-w)*R*sin(th)/K);
    return (z-w)*R*sin(th)/K;
}
//The inner layer of the double integral
double m_By_th_add(double x, double y, double z, double w, int k, double th1, double t){
    double a = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        a += m_by(x, y, z, w, th1+(2*i-1)*t/pow(2,k-1));
    }
    return a*t/(k>2? pow(2,k-2):1);
}
double m_By_th(double x, double y, double z, double w){
    double answer = 0;
    double th1 = 3*E, th2 = 2*PI;
    double t = th2-th1;
    int M = 40;     //the initial M, if it's not enough, we'll give you an alert and you need to enlarge M
    double *Py;
    Py = (double *)malloc(M*M*sizeof(double));

    Py[0] = ( m_by(x, y, z, w,th1)+m_by(x, y, z, w, th2) )*t/2;
    int i = 1, j = 1;
    while(i){
        Py[i*M+0] = (Py[(i-1)*M+0]+m_By_th_add(x,y,z,w,i+1,th1,t))/2;
        while(j){
            Py[i*M+j] = Py[i*M+j-1]+(Py[i*M+j-1]-Py[(i-1)*M+j-1])/(pow(4,j)-1);
            if(j == i){
                j = 1;
                break;
            }else{
                j++;
            }
        }
        //printf("%d,%.8lf\n",i, Py[i*M+i]);
        double e = fabs(Py[i*M+i]-Py[(i-1)*M+i-1]);
        //printf("%.10lf\n\n", e);
        if(e<E&&i!=1){//for the special case            
            answer = Py[i*M+i];
            break;
        }else if(i == M-1){
            printf("x=%lf, y=%lf, z=%lf\n", x, y, z);
            printf("Please enlarge Py!!!\n");
        }else{
            i++;
        }
    }

    free(Py);
    return answer;
}
//the magnet's y component
double m_By_add(double x, double y, double z, int k, double w1, double t){
    double a = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        a += m_By_th(x, y, z, w1+(2*i-1)*t/pow(2,k-1));
    }
    return a*t/(k>2? pow(2,k-2):1);
}
double m_By(double x, double y, double z){
    double answer = 0;
    double w1 = 0, w2 = h;
    double t = w2-w1;
    int M = 40;
    double *Hy;
    Hy = (double *)malloc(M*M*sizeof(double));

    Hy[0] = ( m_By_th(x, y, z, w1)+m_By_th(x, y, z, w2) )*t/2;
    int i = 1, j = 1;
    while(i){
        //printf("%d\n",i);
        Hy[i*M+0] = (Hy[(i-1)*M+0]+m_By_add(x,y,z,i+1,w1,t))/2;
        while(j){
            Hy[i*M+j] = Hy[i*M+j-1]+(Hy[i*M+j-1]-Hy[(i-1)*M+j-1])/(pow(4,j)-1);
            if(j == i){
                j = 1;
                break;
            }else{
                j++;
            }
        }
        double e = fabs(Hy[i*M+i]-Hy[(i-1)*M+i-1]);
        if(e<E&&i!=1){
            answer = Bg*Hy[i*M+i]/(4*PI);
            break;
        }else if(i == M-1){
            printf("x=%lf, y=%lf, z=%lf\n", x, y, z);
            printf("Please enlarge Hy!!!\n");
        }else{
            i++;
        }
    }

    free(Hy);
    return answer;
}

//the component of z axis
double m_bz(double x, double y, double z, double w, double th){
    double K = pow( (x-R*cos(th))*(x-R*cos(th))+(y-R*sin(th))*(y-R*sin(th))+(z-w)*(z-w) , 1.5);
    return ( -(y-R*sin(th))*R*sin(th)-(x-R*cos(th))*R*cos(th) )/K;
}
//The inner layer of the double integral
double m_Bz_th_add(double x, double y, double z, double w, int k, double th1, double t){
    double a = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        a += m_bz(x, y, z, w, th1+(2*i-1)*t/pow(2,k-1));
    }
    return a*t/(k>2? pow(2,k-2):1);
}
double m_Bz_th(double x, double y, double z, double w){
    double answer = 0;
    double th1 = 3*E, th2 = 2*PI;
    double t = th2-th1;
    int M = 40;     //the initial M, if it's not enough, we'll give you an alert and you need to enlarge M
    double *Pz;
    Pz = (double *)malloc(M*M*sizeof(double));

    Pz[0] = ( m_bz(x, y, z, w,th1)+m_bz(x, y, z, w, th2) )*t/2;
    int i = 1, j = 1;
    while(i){
        Pz[i*M+0] = (Pz[(i-1)*M+0]+m_Bz_th_add(x,y,z,w,i+1,th1,t))/2;
        while(j){
            Pz[i*M+j] = Pz[i*M+j-1]+(Pz[i*M+j-1]-Pz[(i-1)*M+j-1])/(pow(4,j)-1);
            if(j == i){
                j = 1;
                break;
            }else{
                j++;
            }
        }
        double e = fabs(Pz[i*M+i]-Pz[(i-1)*M+i-1]);
        if(e<E&&i!=1){
            answer = Pz[i*M+i];
            break;
        }else if(i == M-1){
            printf("x=%lf, y=%lf, z=%lf\n", x, y, z);
            printf("Please enlarge Pz!!!\n");
        }else{
            i++;
        }
    }

    free(Pz);
    return answer;
}
//the magnet's z component
double m_Bz_add(double x, double y, double z, int k, double w1, double t){
    double a = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        a += m_Bz_th(x, y, z, w1+(2*i-1)*t/pow(2,k-1));
    }
    return a*t/(k>2? pow(2,k-2):1);
}
double m_Bz(double x, double y, double z){
    double answer = 0;
    double w1 = 0, w2 = h;
    double t = w2-w1;
    int M = 40;
    double *Hz;
    Hz = (double *)malloc(M*M*sizeof(double));

    Hz[0] = ( m_Bz_th(x, y, z, w1)+m_Bz_th(x, y, z, w2) )*t/2;
    int i = 1, j = 1;
    while(i){
        Hz[i*M+0] = (Hz[(i-1)*M+0]+m_Bz_add(x,y,z,i+1,w1,t))/2;
        while(j){
            Hz[i*M+j] = Hz[i*M+j-1]+(Hz[i*M+j-1]-Hz[(i-1)*M+j-1])/(pow(4,j)-1);
            if(j == i){
                j = 1;
                break;
            }else{
                j++;
            }
        }
        double e = fabs(Hz[i*M+i]-Hz[(i-1)*M+i-1]);
        if(e<E&&i!=1){
            answer = Bg*Hz[i*M+i]/(4*PI);
            break;
        }else if(i == M-1){
            printf("x=%lf, y=%lf, z=%lf\n", x, y, z);
            printf("Please enlarge Hz!!!\n");
        }else{
            i++;
        }
    }

    free(Hz);
    return answer;
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