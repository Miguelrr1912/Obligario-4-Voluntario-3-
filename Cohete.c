#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 6.65e-11
#define MT 5.9736e24
#define ML 0.07349e24
#define d 3.844e8
#define omega 2.6617e-6
#define RT 6.378160e6
#define RL 1.7374e6
#define h 1 //minuto
#define hs (h*60) //segundos

typedef struct {
    double r,phi;
    double pr,pphi;
} cohete;

//hago las funciones para evaluar las derivadas a apartir de f
double fr(double pr)
{
    return pr;
}

double fphi(double pphi,double r)
{
    return pphi/(r*r);
}

double fpr(double r, double phi, double pr, double pphi,double t)
{  
    double rprima, mu, delta;
    delta=G*MT/(d*d*d);
    mu=ML/MT;
    rprima=sqrt(1+r*r - 2*r*cos(phi-omega*t));
    return pphi*pphi/(r*r*r)-delta*(1/(r*r)+mu*1/(rprima*rprima*rprima)*(r-cos(phi-omega*t)));
}

double fpphi(double r, double phi, double pr, double pphi,double t)
{
    double delta,mu, rprima;
    delta=G*MT/(d*d*d);
    mu=ML/MT;
    rprima=sqrt(1+r*r - 2*r*cos(phi-omega*t));
    return -delta*mu*r/(rprima*rprima*rprima)*sin(phi-omega*t);
}

void runge_kuttap(cohete *cohete, FILE *file,double t)
{
    double k[4][4];  

    // Evaluo k1
    k[0][0] = hs * fr(cohete->pr);
    k[1][0] = hs * fphi(cohete->pphi, cohete->r);
    k[2][0] = hs * fpr(cohete->r, cohete->phi, cohete->pr, cohete->pphi, t);
    k[3][0] = hs * fpphi(cohete->r, cohete->phi, cohete->pr, cohete->pphi, t);

    // Evaluo k2
    k[0][1] = hs * fr(cohete->pr + k[2][0]/2);
    k[1][1] = hs * fphi(cohete->pphi + k[3][0]/2, cohete->r + k[0][0]/2);
    k[2][1] = hs * fpr(cohete->r + k[0][0]/2, cohete->phi + k[1][0]/2, cohete->pr + k[2][0]/2, cohete->pphi + k[3][0]/2, t+hs/2);
    k[3][1] = hs * fpphi(cohete->r + k[0][0]/2, cohete->phi + k[1][0]/2, cohete->pr + k[2][0]/2, cohete->pphi + k[3][0]/2, t+hs/2);

    // Evaluo k3
    k[0][2] = hs * fr(cohete->pr + k[2][1]/2);
    k[1][2] = hs * fphi(cohete->pphi + k[3][1]/2, cohete->r + k[0][1]/2);
    k[2][2] = hs * fpr(cohete->r + k[0][1]/2, cohete->phi + k[1][1]/2, cohete->pr + k[2][1]/2, cohete->pphi + k[3][1]/2, t+hs/2);
    k[3][2] = hs * fpphi(cohete->r + k[0][1]/2, cohete->phi + k[1][1]/2, cohete->pr + k[2][1]/2, cohete->pphi + k[3][1]/2, t+hs/2);

    // Evaluo k4
    k[0][3] = hs * fr(cohete->pr + k[2][2]);
    k[1][3] = hs * fphi(cohete->pphi + k[3][2], cohete->r + k[0][2]);
    k[2][3] = hs * fpr(cohete->r + k[0][2], cohete->phi + k[1][2], cohete->pr + k[2][2], cohete->pphi + k[3][2], t+hs);
    k[3][3] = hs * fpphi(cohete->r + k[0][2], cohete->phi + k[1][2], cohete->pr + k[2][2], cohete->pphi + k[3][2], t+hs);

    // Actualizo las variables 
    cohete->r   += (k[0][0] + 2*k[0][1] + 2*k[0][2] + k[0][3]) / 6;
    cohete->phi      += (k[1][0] + 2*k[1][1] + 2*k[1][2] + k[1][3]) / 6;
    cohete->pr  += (k[2][0] + 2*k[2][1] + 2*k[2][2] + k[2][3]) / 6;
    cohete->pphi     += (k[3][0] + 2*k[3][1] + 2*k[3][2] + k[3][3]) / 6;

    // Guardo los resultados
    fprintf(file, "%lf, %lf\n", cohete->r*cos(cohete->phi), cohete->r*sin(cohete->phi));
}



int main (void)
{
    cohete cohete; 
    double xluna, yluna, t, LAT;
    FILE *filecohete, *fileluna;
    int i;
    
    // creo un archivo para el cohete y otro para la luna
    filecohete = fopen("cohete.txt", "w"); 
    fileluna = fopen("luna.txt", "w");

    // Inicializo las variables del cohete;
    LAT=3.1415/2; //
    cohete.r=RT/d;
    cohete.phi=3.1415/2+0.15;
    cohete.pr=sqrt(2*G*MT/RT)*cos(LAT-cohete.phi)/d;
    cohete.pphi=cohete.r*sqrt(2*G*MT/RT)*sin(LAT-cohete.phi)/(d*d);

    
    //coloco coordenadas de la luna (Al principio alineado en el eje y)
    xluna=0;
    yluna=1;
   
    // Imprimo la posición inicial de la luna y del cohete
    fprintf(fileluna, "%lf, %lf\n", xluna, yluna);
    fprintf(filecohete, "%lf, %lf\n", cohete.r*cos(LAT), cohete.r*sin(LAT));

    for(i=0; i<7000; i++)
    {   
        t=i*h*60; // tiempo en segundos
        runge_kuttap(&cohete, filecohete, t);

        // Actualizo la posición de la luna
        xluna =-1*sin(omega*t);
        yluna =1*cos(omega*t);
        fprintf(fileluna, "%lf, %lf\n", xluna, yluna);

       // if(sqrt(cohete.r*cohete.r - 2*cohete.r*cos(cohete.phi-omega*t)) < RL/RT)
       // {
       //     printf("El cohete colisiona con la luna en el tiempo: %lf minutos\n", t/60);
       //     break;
       // }
    }

    return 0;
}