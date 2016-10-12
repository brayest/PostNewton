#include <stdlib.h> /* Librerias Adicionales */
#include <math.h>  /* Librerias de Matematica */
#include </home/brayest/LibreriasProgra/derivadas.h> /*Derivadas*/
#include </home/brayest/LibreriasProgra/metodos.h> /*Rutina de intergracion*/


void PostNewtonNcuerpos ( double x, double **y, double **dvdt, double *M, int N) {
  double G,c;
  int z,k,w,v,r, coordenadas, relaciones;
  
  coordenadas = 3;
  G = 6.6726e-11;
  c = 299792458;  
  
  // x12,x13,x23,...,etc
  relaciones = N*(N-1)/2;

  
  double **Xab, *Rab;
  
  Xab = malloc( relaciones *sizeof(double));
  Rab = malloc( relaciones *sizeof(double));
  for ( z = 0; z < relaciones; z++) {
    Xab[z] = (double *)malloc(coordenadas *sizeof(double));
  }
  
  z = 0;
  k = 0;
  w = 0;
  while ( k < N-1 ) {
//    printf("%d		%d	%d	%d	\n",k,z,N,w);
//    printf("\n");
    if (z >= N-1)
      r = k;
    else r = z;
    for ( w =r+1; w < N; w++) {
//      printf("%d		%d	%d	%d	\n",k,z,N,w);
//      printf("\n");
      for (v = 0; v < coordenadas; v++)
        Xab[z][v] = y[k][v] - y[w][v];
      z +=1;
    }
//    printf("%d		%d	%d	%d	\n",k,z,N,w);
//    printf("\n");
    k++;
  }    
    
    
      
/*  z = 0;
  for ( z = 0; z < relaciones; z++) {
    printf("X%d		",z);
    for ( v = 0; v < coordenadas; v++) {
      printf("%.3f	",Xab[z][v]);
    }
    printf("\n");
  }*/
  
  
  for ( z = 0; z < relaciones; z++) {
    Rab[z]  = sqrt( pow(Xab[z][0],2) + pow(Xab[z][1],2) + pow(Xab[z][2],2) );
//    printf("%.3f \n",Rab[z]);
  }
  
  double accel = 0;
  
  for ( z = 0; z < N; z++) {
    for ( v = 0; v < 3; v++) {
      dvdt[z][v] = y[z][v+3];
    }
    for ( v = 3; v < 6; v++) {
      for ( w = 0; w < relaciones; w++) {
        if ( w != z) {
          accel = accel + (G*M[w]*Xab[w][v])/pow(Rab[w],3);
          printf("%d	\n",w);
        }
      dvdt[z][v] = -accel;
      accel = 0;
      }	
    }
    
  }
  
  

    
}









main (int argc, char *argv[]) {
  
  // Definición numero de cuerpos
  int cuerpos = 2;
  int variables = 6;
  int i,j,k;
  
  // constantes
  double G,c;  
  G = 6.6726e-11;
  c = 299792458;
  
  // Definicion de variables  respecto a cantidad de cuerpos
  double *masas, **iniciales, **movimiento;
  
  // Vectores 
  masas = malloc(cuerpos  *sizeof(double)); 
  iniciales = malloc(cuerpos *sizeof(double));
  movimiento = malloc(cuerpos *sizeof(double));
  
  
  // Matrices : Con la siguiente forma
  //
  // Masa	x	y	z	vx	vy	vz
  //  1		*	*	*	*	*	*
  //  2		*	*	*	*	*	*
  //  3		*	*	*	*	*	*
  
  for ( i = 0; i < cuerpos; i++) {
    iniciales[i] = (double *)malloc( variables *sizeof(double));
    movimiento[i] = (double *)malloc( variables *sizeof(double));
  }
  
  masas[0] = 1.989e30;
  masas[1] = 5.872e24;

  iniciales[0][0] = 0;
  iniciales[0][1] = 0;
  iniciales[0][2] = 0;
  iniciales[0][3] = 0;
  iniciales[0][5] = 0;
  iniciales[1][1] = 0;
  iniciales[1][2] = 0;
  iniciales[1][3] = 0;
  iniciales[1][5] = 0;
  iniciales[1][4] = 28.73e3; //Vy2
  iniciales[0][4] = -(masas[1]/masas[0])*iniciales[1][4];  // Vy1 para conservar momentum
  iniciales[1][0] = G*masas[0]/pow(iniciales[1][4],2); // X2 dependiente de velocidad para Orbita circular ( E = 0 )
  
  
  
  for (i = 0; i < cuerpos; i++) {
    for ( j = 0; j < variables; j++) {
      movimiento[i][j] = iniciales[i][j];
    }
  }
  
  
  
  // Variables de Control
  double dt, tiempo, Nit;
  int print;
  
  dt = strtod(argv[1], NULL);
  print = atoi(argv[2]);
  Nit = atoi(argv[3]);
  tiempo = 0;
  
  FILE *out;
  out = fopen("SolucionNPN.dat","w");
  
  for ( i = 0; i < cuerpos; i++) {
    printf("%d	",i);
    for ( j = 0; j < variables; j++) {
      printf("%.3f	",iniciales[i][j]);
    }
    printf("\n");
  }
  printf("\n"); 
  
  while ( i <= Nit ) {
    if ( i % print == 0 ) {
      for ( j = 0; j < cuerpos; j++) {
        for ( k = 0; k < variables/2; k++) {
          fprintf(out, "%.3f	",movimiento[j][k]);
        }
      }
    fprintf(out,"\n");
    }
    RK4MATRIX(movimiento,iniciales,dt,tiempo,PostNewtonNcuerpos,cuerpos,masas);
    
    for ( j = 0; j < cuerpos; j++)
      for ( k = 0; k < variables; k++)
        iniciales[j][k] = movimiento[j][k];
    i++;
    tiempo += dt;
  }
  fclose(out);
  free(masas);
  free(iniciales);
  free(movimiento);
  
}

  
