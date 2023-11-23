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
#define d 1.004 //cm
#define m 0.0443 //kg
#define N 20

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
    double th1 = 0.01, th2 = 2*PI;
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
    double w1 = 1e-3, w2 = h;
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
    double th1 = 0.01, th2 = 2*PI;
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
    double w1 = 1e-3, w2 = h;
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
    double th1 = 0.01, th2 = 2*PI;
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
    double w1 = 1e-3, w2 = h;
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
}/*
//Get Gauss Integration's weights
void Weigths(double *weights){
    FILE *fp;
    fp = fopen("weights.txt", "r");
    for(int i =0; i < N; i++){
        fscanf(fp, "%lf", &weights[i]);
        //check
        //printf("weights[%d] is %.16lf\n", i, weights[i]);
    }
}
//Get Gauss Integratin's abscissa
void Abscissa(double *abscissa){
    FILE *fp;
    fp = fopen("abscissa.txt", "r");
    for(int i =0; i < N; i++){
        fscanf(fp, "%lf", &abscissa[i]);
        //check
        //printf("abscissa[%d] is %.16lf\n", i, abscissa[i]);
    }
}
//Fx
double m_fx(double b, double u, double al, double *weights, double *abscissa){
    return Bm*( -r*sin(b)*cos(al)*m_Bz(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)) \
    + r*sin(b)*sin(al)*m_By(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)));
}
double m_Fx_th(double u, double al, double *weights, double *abscissa){
    double ans = 0;
    double a = 0, b = 2*PI;
    for(int i = 0; i < N; i++){
        ans += weights[i]*m_fx(((b-a)*abscissa[i]+b+a)/2, u, al, weights, abscissa)*(b-a)/2;
    }
    //check
    //printf("Fx_th(%.4lf, %.4lf, %.4lf) is %.8lf\n", x, y, z, ans);
    return ans;
}
double m_Fx(double al, double *weights, double *abscissa){
    double ans = 0;
    double a = 0, b = 2*l;
    for(int i = 0; i < N; i++){
        ans += weights[i]*m_Fx_th(((b-a)*abscissa[i]+b+a)/2, al, weights, abscissa)*(b-a)/2;
    }
    //check
    //printf("Fx(%.4lf, %.4lf, %.4lf) is %.8lf\n", x, y, z, ans);
    return ans;
}
//Fy
double m_fy(double b, double u, double al, double *weights, double *abscissa){
    return Bm*( -r*sin(b)*sin(al)*m_Bx(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)) \
    + r*cos(b)*m_Bz(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)));
}
double m_Fy_th(double u, double al, double *weights, double *abscissa){
    double ans = 0;
    double a = 0, b = 2*PI;
    for(int i = 0; i < N; i++){
        ans += weights[i]*m_fy(((b-a)*abscissa[i]+b+a)/2, u, al, weights, abscissa)*(b-a)/2;
    }
    //check
    //printf("Fy_th(%.4lf, %.4lf, %.4lf) is %.8lf\n", x, y, z, ans);
    return ans;
}
double m_Fy(double al, double *weights, double *abscissa){
    double ans = 0;
    double a = 0, b = 2*l;
    for(int i = 0; i < N; i++){
        ans += weights[i]*m_Fy_th(((b-a)*abscissa[i]+b+a)/2, al, weights, abscissa)*(b-a)/2;
    }
    //check
    //printf("Fy(%.4lf, %.4lf, %.4lf) is %.8lf\n", x, y, z, ans);
    return ans;
}
//Fz
double m_fz(double b, double u, double al, double *weights, double *abscissa){
    return Bm*( -r*cos(b)*m_By(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)) \
    + r*sin(b)*cos(al)*m_Bx(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)) );
}
double m_Fz_th(double u, double al, double *weights, double *abscissa){
    double ans = 0;
    double a = 0, b = 2*PI;
    for(int i = 0; i < N; i++){
        ans += weights[i]*m_fz(((b-a)*abscissa[i]+b+a)/2, u, al, weights, abscissa)*(b-a)/2;
    }
    //check
    //printf("Fz_th(%.4lf, %.4lf, %.4lf) is %.8lf\n", x, y, z, ans);
    return ans;
}
double m_Fz(double al, double *weights, double *abscissa){
    double ans = 0;
    double a = 0, b = 2*l;
    for(int i = 0; i < N; i++){
        ans += weights[i]*m_Fz_th(((b-a)*abscissa[i]+b+a)/2, al, weights, abscissa)*(b-a)/2;
    }
    //check
    //printf("Fz(%.4lf, %.4lf, %.4lf) is %.8lf\n", x, y, z, ans);
    return ans;
}*/

//Fx
float m_fx(float b, float u, float a){
    return ( -r*sin(b)*cos(a)*m_Bz(-r*sin(b),r*cos(b)*cos(a)+r*sin(a)-r*cos(a),r*cos(b)*sin(a)-d-u*cos(a)-r*sin(a)) \
    + r*sin(b)*sin(a)*m_By(-r*sin(b),r*cos(b)*cos(a)+r*sin(a)-r*cos(a),r*cos(b)*sin(a)-d-u*cos(a)-r*sin(a)));
}
float m_Fx_b_add(int k, float b1, float t, float u, float a){
    float ans = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        ans += m_fx(b1+(2*i-1)*t/pow(2,k-1),u,a);
    }
    return ans*t/(k>2? pow(2,k-2):1);
}
float m_Fx_b(float u, float a){
    float answer = 0;
    float b1 = 0, b2 = 2*PI;
    float t = b2-b1;
    int M = 40;
    float *Px;
    Px = (float *)malloc(M*M*sizeof(float));
    Px[0] = ( m_fx(b1,u,a)+m_fx(b2,u,a) )*t/2;
    int i = 1, j = 1;
    while(i){
        Px[i*M+0] = (Px[(i-1)*M+0]+m_Fx_b_add(i+1,b1,t,u,a))/2;
        while(j){
            Px[i*M+j] = Px[i*M+j-1]+(Px[i*M+j-1]-Px[(i-1)*M+j-1])/(pow(4,j)-1);
            if(j==i){
                j=i;
                break;
            }else{
                j++;
            }
        }
        float e = fabs(Px[i*M+i]-Px[(i-1)*M+i-1]);
        if( e<E/1000&&i!=1 ){
            answer = Px[i*M+i];
            break;
        }else if(i==M-1){
            printf("u = %lf, a = %lf.\n", u, a);
            printf("Please enlarge Px!");
        }else{
            i++;
        }
    }
    free(Px);
    return answer;
}
float m_Fx_add(int k, float u1, float t, float a){
    float ans = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        ans += m_Fx_b(u1+(2*i-1)*t/pow(2,k-1),a);
    }
    return ans*t/(k>2? pow(2,k-2):1);
}
float m_Fx(float a){
    float answer = 0;
    float u1 = 0, u2 = 2*l;
    float t = u2-u1;
    int M = 40;
    float *Hx;
    Hx = (float *)malloc(M*M*sizeof(float));
    Hx[0] = ( m_Fx_b(u1,a)+m_Fx_b(u2,a) )*t/2;
    int i = 1, j = 1;
    while(i){
        Hx[i*M+0] = (Hx[(i-1)*M+0]+m_Fx_add(i+1,u1,t,a))/2;
        while(j){
            Hx[i*M+j] = Hx[i*M+j-1]+(Hx[i*M+j-1]-Hx[(i-1)*M+j-1])/(pow(4,j)-1);
            if(j==i){
                j=i;
                break;
            }else{
                j++;
            }
        }
        float e = fabs(Hx[i*M+i]-Hx[(i-1)*M+i-1]);
        if( e<E/1000&&i!=1 ){
            answer = Bm*Hx[i*M+i];
            break;
        }else if(i==M-1){
            printf("a = %lf.\n", a);
            printf("Please enlarge Hx!");
        }else{
            i++;
        }
    }
    free(Hx);
    return answer;
}

//Fy
float m_fy(float b, float u, float a){
    return ( -r*sin(b)*sin(a)*m_Bx(-r*sin(b),r*cos(b)*cos(a)+r*sin(a)-r*cos(a),r*cos(b)*sin(a)-d-u*cos(a)-r*sin(a)) \
    + r*cos(b)*m_Bz(-r*sin(b),r*cos(b)*cos(a)+r*sin(a)-r*cos(a),r*cos(b)*sin(a)-d-u*cos(a)-r*sin(a)));
}
float m_Fy_b_add(int k, float b1, float t, float u, float a){
    float ans = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        ans += m_fy(b1+(2*i-1)*t/pow(2,k-1),u,a);
    }
    return ans*t/(k>2? pow(2,k-2):1);
}
float m_Fy_b(float u, float a){
    float answer = 0;
    float b1 = 0, b2 = 2*PI;
    float t = b2-b1;
    int M = 40;
    float *Py;
    Py = (float *)malloc(M*M*sizeof(float));
    Py[0] = ( m_fy(b1,u,a)+m_fy(b2,u,a) )*t/2;
    int i = 1, j = 1;
    while(i){
        Py[i*M+0] = (Py[(i-1)*M+0]+m_Fy_b_add(i+1,b1,t,u,a))/2;
        while(j){
            Py[i*M+j] = Py[i*M+j-1]+(Py[i*M+j-1]-Py[(i-1)*M+j-1])/(pow(4,j)-1);
            if(j==i){
                j=i;
                break;
            }else{
                j++;
            }
        }
        float e = fabs(Py[i*M+i]-Py[(i-1)*M+i-1]);
        if( e<E/1000&&i!=1 ){
            answer = Py[i*M+i];
            break;
        }else if(i==M-1){
            printf("u = %lf, a = %lf.\n", u, a);
            printf("Please enlarge Py!");
        }else{
            i++;
        }
    }
    free(Py);
    return answer;
}
float m_Fy_add(int k, float u1, float t, float a){
    float ans = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        ans += m_Fy_b(u1+(2*i-1)*t/pow(2,k-1),a);
    }
    return ans*t/(k>2? pow(2,k-2):1);
}
float m_Fy(float a){
    float answer = 0;
    float u1 = 0, u2 = 2*l;
    float t = u2-u1;
    int M = 40;
    float *Hy;
    Hy = (float *)malloc(M*M*sizeof(float));
    Hy[0] = ( m_Fy_b(u1,a)+m_Fy_b(u2,a) )*t/2;
    int i = 1, j = 1;
    while(i){
        Hy[i*M+0] = (Hy[(i-1)*M+0]+m_Fy_add(i+1,u1,t,a))/2;
        while(j){
            Hy[i*M+j] = Hy[i*M+j-1]+(Hy[i*M+j-1]-Hy[(i-1)*M+j-1])/(pow(4,j)-1);
            if(j==i){
                j=i;
                break;
            }else{
                j++;
            }
        }
        float e = fabs(Hy[i*M+i]-Hy[(i-1)*M+i-1]);
        if( e<E/1000&&i!=1 ){
            answer = Bm*Hy[i*M+i];
            break;
        }else if(i==M-1){
            printf("a = %lf.\n", a);
            printf("Please enlarge Hy!");
        }else{
            i++;
        }
    }
    free(Hy);
    return answer;
}

//Fz
float m_fz(float b, float u, float a){
    return ( -r*cos(b)*m_By(-r*sin(b),r*cos(b)*cos(a)+r*sin(a)-r*cos(a),r*cos(b)*sin(a)-d-u*cos(a)-r*sin(a)) \
    + r*sin(b)*cos(a)*m_Bx(-r*sin(b),r*cos(b)*cos(a)+r*sin(a)-r*cos(a),r*cos(b)*sin(a)-d-u*cos(a)-r*sin(a)) );
}
float m_Fz_b_add(int k, float b1, float t, float u, float a){
    float ans = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        ans += m_fz(b1+(2*i-1)*t/pow(2,k-1),u,a);
    }
    return ans*t/(k>2? pow(2,k-2):1);
}
float m_Fz_b(float u, float a){
    float answer = 0;
    float b1 = 0, b2 = 2*PI;
    float t = b2-b1;
    int M = 40;
    float *Pz;
    Pz = (float *)malloc(M*M*sizeof(float));
    Pz[0] = ( m_fz(b1,u,a)+m_fz(b2,u,a) )*t/2;
    int i = 1, j = 1;
    while(i){
        Pz[i*M+0] = (Pz[(i-1)*M+0]+m_Fz_b_add(i+1,b1,t,u,a))/2;
        while(j){
            Pz[i*M+j] = Pz[i*M+j-1]+(Pz[i*M+j-1]-Pz[(i-1)*M+j-1])/(pow(4,j)-1);
            if(j==i){
                j=i;
                break;
            }else{
                j++;
            }
        }
        float e = fabs(Pz[i*M+i]-Pz[(i-1)*M+i-1]);
        if( e<E/1000&&i!=1 ){
            answer = Pz[i*M+i];
            break;
        }else if(i==M-1){
            printf("u = %lf, a = %lf.\n", u, a);
            printf("Please enlarge Pz!");
        }else{
            i++;
        }
    }
    free(Pz);
    return answer;
}
float m_Fz_add(int k, float u1, float t, float a){
    float ans = 0;
    for(int i = 1; i <= (k>2? pow(2,k-2):1); i++){
        ans += m_Fz_b(u1+(2*i-1)*t/pow(2,k-1),a);
    }
    return ans*t/(k>2? pow(2,k-2):1);
}
float m_Fz(float a){
    float answer = 0;
    float u1 = 0, u2 = 2*l;
    float t = u2-u1;
    int M = 40;
    float *Hz;
    Hz = (float *)malloc(M*M*sizeof(float));
    Hz[0] = ( m_Fz_b(u1,a)+m_Fz_b(u2,a) )*t/2;
    int i = 1, j = 1;
    while(i){
        Hz[i*M+0] = (Hz[(i-1)*M+0]+m_Fz_add(i+1,u1,t,a))/2;
        while(j){
            Hz[i*M+j] = Hz[i*M+j-1]+(Hz[i*M+j-1]-Hz[(i-1)*M+j-1])/(pow(4,j)-1);
            if(j==i){
                j=i;
                break;
            }else{
                j++;
            }
        }
        //printf("i = %d,%lf\n", i, Hz[i*M+i]);
        float e = fabs(Hz[i*M+i]-Hz[(i-1)*M+i-1]);
        if( e<E/1000&&i!=1 ){
            answer = Bm*Hz[i*M+i];
            break;
        }else if(i==M-1){
            printf("a = %lf.\n", a);
            printf("Please enlarge Hz!");
        }else{
            i++;
        }
    }
    free(Hz);
    return answer;
}
/*
//fx
double m_fx(double b, double u, double al){
    return Bm*( -r*sin(b)*cos(al)*m_Bz(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)) \
    + r*sin(b)*sin(al)*m_By(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)));
}
double m_fx_b(double a, double b, double u, double al)
{
    double c = a+(b-a)/2;
    double p1 = (m_fx(a,u,al)+4*m_fx(c,u,al)+m_fx(b,u,al))*(b-a)/6;
    double p2 = (m_fx(a,u,al)+4*m_fx((a+c)/2,u,al)+2*m_fx(c,u,al)+4*m_fx((c+b)/2,u,al)+m_fx(b,u,al))*(b-a)/12;
    printf("p1 = %.8lf, p2 = %.8lf.\n", p1, p2);
    if(fabs(p1-p2)<=15*E)
        return p2+(p2-p1)/15;
    return m_fx_b(a,c,u,al)+m_fx_b(c,b,u,al);
}
double m_fx_u(double a, double b, double al)
{
    double ab = 0, bb = 2*PI;
    double c = a+(b-a)/2;
    double p1 = (m_fx_b(ab,bb,a,al)+4*m_fx_b(ab,bb,c,al)+m_fx_b(ab,bb,b,al))*(b-a)/6;
    double p2 = (m_fx_b(ab,bb,a,al)+4*m_fx_b(ab,bb,(a+c)/2,al)+2*m_fx_b(ab,bb,c,al)+4*m_fx_b(ab,bb,(c+b)/2,al)+m_fx_b(ab,bb,b,al))*(b-a)/12;
    if(fabs(p1-p2)<=15*E)
        return p2+(p2-p1)/15;
    return m_fx_u(a,c,al)+m_fx_u(c,b,al);
}
double m_Fx(double al){
    double au = 0, bu = 2*l;
    return m_fx_u(au,bu,al);
}

//fy
double m_fy(double b, double u, double al){
    return Bm*( -r*sin(b)*sin(al)*m_Bx(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)) \
    + r*cos(b)*m_Bz(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)));
}
double m_fy_b(double a, double b, double u, double al)
{
    double c = a+(b-a)/2;
    double p1 = (m_fy(a,u,al)+4*m_fy(c,u,al)+m_fy(b,u,al))*(b-a)/6;
    double p2 = (m_fy(a,u,al)+4*m_fy((a+c)/2,u,al)+2*m_fy(c,u,al)+4*m_fy((c+b)/2,u,al)+m_fy(b,u,al))*(b-a)/12;
    if(fabs(p1-p2)<=15*E)
        return p2+(p2-p1)/15;
    return m_fy_b(a,c,u,al)+m_fy_b(c,b,u,al);
}
double m_fy_u(double a, double b, double al)
{
    double ab = 0, bb = 2*PI;
    double c = a+(b-a)/2;
    double p1 = (m_fy_b(ab,bb,a,al)+4*m_fy_b(ab,bb,c,al)+m_fy_b(ab,bb,b,al))*(b-a)/6;
    double p2 = (m_fy_b(ab,bb,a,al)+4*m_fy_b(ab,bb,(a+c)/2,al)+2*m_fy_b(ab,bb,c,al)+4*m_fy_b(ab,bb,(c+b)/2,al)+m_fy_b(ab,bb,b,al))*(b-a)/12;
    if(fabs(p1-p2)<=15*E)
        return p2+(p2-p1)/15;
    return m_fy_u(a,c,al)+m_fy_u(c,b,al);
}
double m_Fy(double al){
    double au = 0, bu = 2*l;
    return m_fy_u(au,bu,al);
}
//fz
double m_fz(double b, double u, double al){
    return Bm*( -r*cos(b)*m_By(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)) \
    + r*sin(b)*cos(al)*m_Bx(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)));
}
double m_fz_b(double a, double b, double u, double al)
{
    double c = a+(b-a)/2;
    double p1 = (m_fz(a,u,al)+4*m_fz(c,u,al)+m_fz(b,u,al))*(b-a)/6;
    double p2 = (m_fz(a,u,al)+4*m_fz((a+c)/2,u,al)+2*m_fz(c,u,al)+4*m_fz((c+b)/2,u,al)+m_fz(b,u,al))*(b-a)/12;
    if(fabs(p1-p2)<=15*E)
        return p2+(p2-p1)/15;
    return m_fz_b(a,c,u,al)+m_fz_b(c,b,u,al);
}
double m_fz_u(double a, double b, double al)
{
    double ab = 0, bb = 2*PI;
    double c = a+(b-a)/2;
    double p1 = (m_fz_b(ab,bb,a,al)+4*m_fz_b(ab,bb,c,al)+m_fz_b(ab,bb,b,al))*(b-a)/6;
    double p2 = (m_fz_b(ab,bb,a,al)+4*m_fz_b(ab,bb,(a+c)/2,al)+2*m_fz_b(ab,bb,c,al)+4*m_fz_b(ab,bb,(c+b)/2,al)+m_fz_b(ab,bb,b,al))*(b-a)/12;
    if(fabs(p1-p2)<=15*E)
        return p2+(p2-p1)/15;
    return m_fz_u(a,c,al)+m_fz_u(c,b,al);
}
double m_Fz(double al){
    double au = 0, bu = 2*l;
    return m_fz_u(au,bu,al);
}*/
/*
//fx
double m_fx(double b, double u, double al){
    return ( -r*sin(b)*cos(al)*m_Bz(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)) \
    + r*sin(b)*sin(al)*m_By(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)));
}
double m_fx_b(double u, double al){
	int i=1,j,k,n=8;
	double T[8],a=0.0,b=2*PI,s=0.0; 
	T[0]=0.5*(b-a)*(m_fx(a,u,al)+m_fx(b,u,al)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_fx(a+(2*k-1)*(b-a)/pow(2,j),u,al); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[7]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[7]-T[0])>E/100;i++) 
	{ 
		T[0]=T[7]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[7]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
}
double m_Fx(double al){
	int i=1,j,k,n=8;
	double T[8],a=0.0,b=l,s=0.0; 
	T[0]=0.5*(b-a)*(m_fx_b(a,al)+m_fx_b(b,al)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_fx_b(a+(2*k-1)*(b-a)/pow(2,j),al); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[7]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[7]-T[0])>E/100;i++) 
	{ 
		T[0]=T[7]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[7]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
	return Bm*T[7];
}
//fy
double m_fy(double b, double u, double al){
    return ( -r*sin(b)*sin(al)*m_Bx(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)) \
    + r*cos(b)*m_Bz(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)));
}
double m_fy_b(double u, double al){
	int i=1,j,k,n=8;
	double T[8],a=0.0,b=2*PI,s=0.0; 
	T[0]=0.5*(b-a)*(m_fy(a,u,al)+m_fy(b,u,al)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_fy(a+(2*k-1)*(b-a)/pow(2,j),u,al); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[7]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[7]-T[0])>E/100;i++) 
	{ 
		T[0]=T[7]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[7]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
	return T[7];
}
double m_Fy(double al){
	int i=1,j,k,n=8;
	double T[8],a=0.0,b=l,s=0.0; 
	T[0]=0.5*(b-a)*(m_fy_b(a,al)+m_fy_b(b,al)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_fy_b(a+(2*k-1)*(b-a)/pow(2,j),al); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[7]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[7]-T[0])>E/100;i++) 
	{ 
		T[0]=T[7]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[7]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
	return Bm*T[7];
}
//fz
double m_fz(double b, double u, double al){
    return ( -r*cos(b)*m_By(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)) \
    + r*sin(b)*cos(al)*m_Bx(-r*sin(b),r*cos(b)*cos(al)+r*sin(al)-r*cos(al),r*cos(b)*sin(al)-d-u*cos(al)-r*sin(al)));
}
double m_fz_b(double u, double al){
	int i=1,j,k,n=8;
	double T[8],a=0.0,b=2*PI,s=0.0; 
	T[0]=0.5*(b-a)*(m_fz(a,u,al)+m_fz(b,u,al)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_fz(a+(2*k-1)*(b-a)/pow(2,j),u,al); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[7]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[7]-T[0])>E/100;i++) 
	{ 
		T[0]=T[7]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[7]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
	return T[7];
}
double m_Fz(double al){
	int i=1,j,k,n=8;
	double T[8],a=0.0,b=l,s=0.0; 
	T[0]=0.5*(b-a)*(m_fz_b(a,al)+m_fz_b(b,al)); 
	for(j=1;j<n-1;j++) 
	{ 
		for(k=1;k<=pow(2,j-1);k++) 
			s+=m_fz_b(a+(2*k-1)*(b-a)/pow(2,j),al); 
		T[j]=0.5*(T[j-1]+(b-a)*s/pow(2,j-1)); 
		s=0.0; 
	}
	T[7]=(4*T[1]-T[0])/(double)3; 
	for(;fabs(T[7]-T[0])>E/100;i++) 
	{ 
		T[0]=T[7]; 
		for(j=1;j<n-1-i;j++) 
			T[j]=(pow(4,i)*T[j+1]-T[j])/(pow(4,i)-1); 
		T[7]=(pow(4,i+1)*T[1]-T[0])/(pow(4,i+1)-1); 
	}
	return Bm*T[7];
}*/
#endif