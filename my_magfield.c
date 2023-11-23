#include"my_headfile.h"

#define A 0.1

void main(void){
    FILE *fp;
    fp = fopen("my_magnetic.txt", "w");
    fprintf(fp, "          x           y           z            Bx            By            Bz\n");
    for(double x = -3.00; x <= 3.01; x += A){
        for(double y = -3.00; y <= 3.01; y += A){
            for(double z = -16.00; z < 0.0000; z += A*10){
                fprintf(fp, "%11.4lf,%11.4lf,%11.4lf,%13.8lf,%13.8lf,%13.8lf\n", x, y, z,\
                m_Bx(x, y, z), m_By(x, y, z), m_Bz(x, y, z));
            }
        }
    }
}