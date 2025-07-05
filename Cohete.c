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

void runge_kuttap(cohete *cohete, FILE *file,double t,double hs)
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

void runge_kuttaadaptado(cohete *cohetei, FILE *file,FILE *aux,double *t, double *ha, double epsmax)           
{   
    //El que va a ir con h/2
    cohete coheteh2, coheteh;
    double error[4],errormax, hmax,s,saux;
    int i;
    
    coheteh=*cohetei; //Guardo el cohete original 
    coheteh2=*cohetei; //Guardo el cohete original

    runge_kuttap(&coheteh,aux,*t,*ha);//que se mande a un archivo

    //Ahora lo hago igual pero con h/2
     //Apunto al cohete que va a ir con h/2

    runge_kuttap(&coheteh2,aux,*t,*ha/2); //Lo evalúo con h/2
     //Ahora añado el otro h/2 
    runge_kuttap(&coheteh2,aux, *t+*ha/2,*ha/2); //Lo evalúo con el otro h/2

    //Ahora puedo calcular el error
    error[0]=16*fabs(coheteh2.r - coheteh.r)/15; //Error en r 
    error[1]=16*fabs(coheteh2.phi - coheteh.phi)/15; //Error en phi
    error[2]=16*fabs(coheteh2.pr - coheteh.pr)/15; //Error en pr 
    error[3]=16*fabs(coheteh2.pphi - coheteh.pphi)/15; //Error en pphi

    //Veo cual es el error máximo 
    errormax=error[0];
    for(i=1; i<4;i++)
    {
        if(error[i]>errormax)
            errormax=error[i];
    }
    
    s=pow(errormax/epsmax,0.2);

    if(s< 1e-8) s=1e-8; //Si el error es mas pequeño que ese valor lo sustituyo.

    hmax=*ha/s;//Calculo el hmax

    if(s>2)      
   {*ha=hmax;
     return;} //si s mayor que 2 sustituto
    else if(s<2)
    {*t+=*ha;
     *cohetei=coheteh;
     fprintf(file, "%lf, %lf\n", cohetei->r*cos(cohetei->phi), cohetei->r*sin(cohetei->phi));
    }

    
    if(*ha<hmax) *ha=2*(*ha); //Si el h es menor que hmax lo sustituyo por 2*h
    
    return;
}

double Hamiltoniano(cohete coheteh,double t)
{  
    // Defino las cosas
   double delta,mu,rprima,H;
   delta=G*MT/(d*d);
   mu=ML/MT;
    rprima=sqrt(1+coheteh.r*coheteh.r - 2*coheteh.r*cos(coheteh.phi-omega*t));
    H=coheteh.pr*coheteh.pr*0.5+coheteh.pphi*coheteh.pphi*0.5/(coheteh.r*coheteh.r)-delta*(1/coheteh.r+mu/rprima); //t=0 para el Hamiltoniano
    return H;
}

int main (void)
{
    cohete cohete; 
    double xluna, yluna, t, thetha,Hprima;
    FILE *filecohete, *fileluna,*filecohetead, *filelunaad, *fileaux, *fileHprima;
    int i;
    
    // creo un archivo para el cohete y otro para la luna
    filecohete = fopen("cohete.txt", "w"); 
    fileluna = fopen("luna.txt", "w");
    filecohetead = fopen("coheteaux.txt", "w"); 
    filelunaad = fopen("lunaaux.txt", "w");
    fileaux = fopen("aux.txt", "w");
    fileHprima = fopen("Hprima.txt", "w"); // Archivo auxiliar para runge_kuttaadaptado

    // Inicializo las variables del cohete;
    thetha=2.2; //
    cohete.r=RT/d;
    cohete.phi=2.1;
    cohete.pr=sqrt(2*G*MT/RT)*cos(thetha-cohete.phi)/d; // Velocidad radial inicial
    cohete.pphi=cohete.r*sqrt(2*G*MT/RT)*sin(thetha-cohete.phi)/(d*d);

    
    //coloco coordenadas de la luna (Al principio alineado en el eje y)
    xluna=0;
    yluna=1;
   
    // Imprimo la posición inicial de la luna y del cohete
    fprintf(fileluna, "%lf, %lf\n", xluna, yluna);
    fprintf(filecohete, "%lf, %lf\n", cohete.r*cos(cohete.phi), cohete.r*sin(cohete.phi));
    fprintf(filelunaad, "%lf, %lf\n", xluna, yluna);
    fprintf(filecohetead, "%lf, %lf\n", cohete.r*cos(cohete.phi), cohete.r*sin(cohete.phi));

    for(i=0; i<7000; i++)
    {   
        t=i*h*60; // tiempo en segundos
        runge_kuttap(&cohete, filecohete, t,h*60); // paso de tiempo en segundos

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

    // Inicializo las variables del cohete(otra vez);
    thetha=2.2; //
    cohete.r=RT/d;
    cohete.phi=2.1;
    cohete.pr=sqrt(2*G*MT/RT)*cos(thetha-cohete.phi)/d;
    cohete.pphi=cohete.r*sqrt(2*G*MT/RT)*sin(thetha-cohete.phi)/(d*d);

    
    //coloco coordenadas de la luna (Al principio alineado en el eje y)
    xluna=0;
    yluna=1;


    t = 0; // Reinicio el tiempo para el método adaptado
    double tmax= 7000 * h * 60; // tiempo máximo en segundos
    double hadaptado = h * 60; // paso de tiempo adaptado en segundos
    while(t<tmax)
    { runge_kuttaadaptado(&cohete, filecohetead,fileaux, &t, &hadaptado, 1e-6); // paso de tiempo adaptado en segundos
      // Actualizo la posición de la luna
      xluna =-1*sin(omega*t);
      yluna =1*cos(omega*t);
      fprintf(filelunaad, "%lf, %lf\n", xluna, yluna);
      Hprima =Hamiltoniano(cohete, t)-cohete.pphi*omega; // Calculo el Hamiltoniano
      fprintf(fileHprima, "%.12lf\n", Hprima); // Guardo el Hamiltoniano


    }

    // Cierro los archivos
    fclose(filecohete);
    fclose(fileluna);
    fclose(filecohetead);
    fclose(filelunaad);


    return 0;
}