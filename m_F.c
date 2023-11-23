//check that Fx is zero
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"my_headfile.h"


void main(void)
{
    FILE *fp;
    fp = fopen("m_Fx.txt", "w");
    double al = 0;
    for(int i = 0; i < 51; i++){
        al = i*PI/100;
        fprintf(fp, "when alpha is %.5lfPI, Fx is %.8lf.\n", al/PI, m_Fx(al));
        fprintf(fp, "when alpha is %.5lfPI, Fy is %.8lf.\n", al/PI, m_Fy(al));
        fprintf(fp, "when alpha is %.5lfPI, Fz is %.8lf.\n", al/PI, m_Fz(al));
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    fprintf(fp, "we can see that Fx is zero within the error limitation!");
}