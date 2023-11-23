#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"my_headfile.h"


void main(void)
{
	FILE *fp;
    fp = fopen("test_f_b.txt", "w");
	//for(double b=0; b<=2*PI; b+=0.5){
		for(double a=0; a<=PI/2; a+=0.5){
			for(double u=0; u<=2*l; u+=0.2){
				fprintf(fp, "x,%lf\n", m_Fx_b(u,a));
				fprintf(fp, "y,%lf\n", m_Fy_b(u,a));
				fprintf(fp, "z,%lf\n", m_Fz_b(u,a));
			}		
		}
	//}
}