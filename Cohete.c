#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define h 0.01
#define g 9.81

typedef struct {
    double thetha,phi;
    double thethap,phip;
} cuerpos;

//hago las funciones para evaluar las derivadas a apartir de f
double ftheta(double theta)
{
    return theta;
}

double fphi(double phi)
{
    return phi;
}

double fthethap(double thetha, double phi, double thethap, double phip)
{
    return (g*sin(phi)*cos(-phi+thetha) - 2*g*sin(thetha) - pow(phip, 2)*cos(-phi+thetha)*sin(-phi+thetha) - pow(phip, 2)*sin(-phi + thetha))/(2 - pow(cos(-phi+thetha), 2));
}

double fphip(double thetha, double phi, double thethap, double phip)
{
    return (g*sin(thetha)*cos(-phi+thetha) - g*sin(phi) + 0.5*pow(phip, 2)*cos(-phi+thetha)*sin(-phi+thetha) + pow(phip, 2)*sin(-phi+thetha))/(1 - 0.5*pow(cos(-phi+thetha), 2));
}


void runge_kuttap(cuerpos *cuerpo, FILE *file)
{
    double k[4][4];  

    // Evaluo k1
    k[0][0] = h * ftheta(cuerpo->thethap);
    k[1][0] = h * fphi(cuerpo->phip);
    k[2][0] = h * fthethap(cuerpo->thetha, cuerpo->phi, cuerpo->thethap, cuerpo->phip);
    k[3][0] = h * fphip(cuerpo->thetha, cuerpo->phi, cuerpo->thethap, cuerpo->phip);

    // Evaluo k2
    k[0][1] = h * ftheta(cuerpo->thethap + k[2][0]/2);
    k[1][1] = h * fphi(cuerpo->phip + k[3][0]/2);
    k[2][1] = h * fthethap(cuerpo->thetha + k[0][0]/2, cuerpo->phi + k[1][0]/2, cuerpo->thethap + k[2][0]/2, cuerpo->phip + k[3][0]/2);
    k[3][1] = h * fphip(cuerpo->thetha + k[0][0]/2, cuerpo->phi + k[1][0]/2, cuerpo->thethap + k[2][0]/2, cuerpo->phip + k[3][0]/2);

    // Evaluo k3
    k[0][2] = h * ftheta(cuerpo->thethap + k[2][1]/2);
    k[1][2] = h * fphi(cuerpo->phip + k[3][1]/2);
    k[2][2] = h * fthethap(cuerpo->thetha + k[0][1]/2, cuerpo->phi + k[1][1]/2, cuerpo->thethap + k[2][1]/2, cuerpo->phip + k[3][1]/2);
    k[3][2] = h * fphip(cuerpo->thetha + k[0][1]/2, cuerpo->phi + k[1][1]/2, cuerpo->thethap + k[2][1]/2, cuerpo->phip + k[3][1]/2);

    // Evaluo k4
    k[0][3] = h * ftheta(cuerpo->thethap + k[2][2]);
    k[1][3] = h * fphi(cuerpo->phip + k[3][2]);
    k[2][3] = h * fthethap(cuerpo->thetha + k[0][2], cuerpo->phi + k[1][2], cuerpo->thethap + k[2][2], cuerpo->phip + k[3][2]);
    k[3][3] = h * fphip(cuerpo->thetha + k[0][2], cuerpo->phi + k[1][2], cuerpo->thethap + k[2][2], cuerpo->phip + k[3][2]);

    // Actualizo las variables 
    cuerpo->thetha   += (k[0][0] + 2*k[0][1] + 2*k[0][2] + k[0][3]) / 6;
    cuerpo->phi      += (k[1][0] + 2*k[1][1] + 2*k[1][2] + k[1][3]) / 6;
    cuerpo->thethap  += (k[2][0] + 2*k[2][1] + 2*k[2][2] + k[2][3]) / 6;
    cuerpo->phip     += (k[3][0] + 2*k[3][1] + 2*k[3][2] + k[3][3]) / 6;

    // Guardo los resultados
    fprintf(file, "%lf, %lf, %lf, %lf\n", cuerpo->thetha, cuerpo->phi, cuerpo->thethap, cuerpo->phip);
}



int main (void)
{
    double desv, desvsum=0, coeflyu=0; 
    double energias[] = {1, 3, 5, 10, 15};
    int num_energias = 5,j,E;
    cuerpos cuerpo;
    cuerpos cuerpoperturbado;

    FILE *file, *file2, *file3, *file4;
    file = fopen("Pendulo_dobleA2.txt", "w");
    if (file == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    file2 = fopen("Pendulo_doblePerturbado.txt", "w");
    if (file2 == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    file3 = fopen("Lyupanov.txt", "w");
    if (file3 == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    
    //Lo hago para 5 energías diferentes y paralelizo
    #pragma omp parallel for private(j)
    for(j=0; j<num_energias;j++)
    {
    E=energias[j];
    cuerpo.thetha=0.15;
    cuerpo.phi=0.08;
    cuerpo.thethap=sqrt(2*(E+9.8*(2*cos(cuerpo.thetha)+cos(cuerpo.phi))-3*9.8));
    //cuerpo.thethap=0.0;
    cuerpo.phip=0.0;

    desv=0.007;
    cuerpoperturbado.thetha=cuerpo.thetha;
    cuerpoperturbado.phi=cuerpo.phi+desv;    
    cuerpoperturbado.thethap=cuerpo.thethap;
    cuerpoperturbado.phip=cuerpo.phip;
    

    //Voy a guardar los resultados en un archivo
    #pragma omp critica 
    { fprintf(file, "%lf, %lf, %lf, %lf\n", cuerpo.thetha, cuerpo.phi, cuerpo.thethap, cuerpo.phip);
    fprintf(file2, "%lf, %lf, %lf, %lf\n", cuerpoperturbado.thetha, cuerpoperturbado.phi, cuerpoperturbado.thethap, cuerpoperturbado.phip);
    }

    for(double i=0; i<100; i += h)
    {
        runge_kuttap(&cuerpo, file);
        runge_kuttap(&cuerpoperturbado,file2);

        //Calculo la desviación 
        desvsum=sqrt(pow((cuerpo.thetha-cuerpoperturbado.thetha),2)+
                     pow((cuerpo.phi-cuerpoperturbado.phi),2)+
                     pow((cuerpo.thethap-cuerpoperturbado.thethap),2)+
                     pow((cuerpo.phip-cuerpoperturbado.phip),2));

        

        //Obtengo el coeficiente de lyupanov;
        coeflyu+=log(fabs(desvsum)/fabs(desv));

        //cuerpoperturbado.thetha=cuerpo.thetha+desv*(cuerpoperturbado.thetha-cuerpo.thetha)/(abs(cuerpoperturbado.thetha-cuerpo.thetha));
        #pragma omp critica
        fprintf(file3,"%lf\n",coeflyu/(h+i+h)); //Guardo el coeficiente de lyupanov en el archivo
    }
    coeflyu=coeflyu/(100*h);
    printf("%f\n",coeflyu);

    coeflyu=0; //Reinicio el coeficiente de lyupanov para la siguiente energía
    desvsum=0; //Reinicio la suma de desviaciones para la siguiente energía
    //Hago entre una energia y otra una separación en el archivo 
    fprintf(file3, "\n");
   
}

    fclose(file);
    fclose(file2);
   
    

    return 0;
}