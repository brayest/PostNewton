#include <stdlib.h> /* Librerias Adicionales */
#include <math.h>  /* Librerias de Matematica */
#include </home/brayest/LibreriasProgra/derivadas.h> /*Derivadas*/
#include </home/brayest/LibreriasProgra/metodos.h> /*Rutina de intergracion*/


void PostNewtonNcuerposV2 ( double x, double **y, double **dvdt, double *M, int N) {
  double G,c;
  int z,k,w,v,r, coordenadas, relaciones;
  
  coordenadas = 3;
  G = 6.6726e-11;
  c = 299792458;  
  
  double ***Xab,**Rab,***Vab;
  
  Xab = malloc( N * sizeof(double));
  Rab = malloc( N * sizeof(double));
  Vab = malloc( N * sizeof(double));
  for ( z = 0; z < N ; z++) {
    Xab[z] = (double **)malloc( N * sizeof(double));
    Rab[z] = (double *)malloc(N * sizeof(double));
    Vab[z] = (double **)malloc( N * sizeof(double));
    for ( v = 0; v < coordenadas; v++) {
      Xab[z][v] = (double *)malloc( coordenadas * sizeof(double));
      Vab[z][v] = (double *)malloc( coordenadas * sizeof(double));
    }
  }
    
  
  for ( z = 0; z < N; z++)
    for ( v = 0; v < N; v++) {
      for ( w = 0; w < coordenadas; w++) {
        Xab[z][v][w] = y[z][w] - y[v][w];
        Vab[z][v][w] = y[z][w+3] - y[v][w+3];
      }
      Rab[z][v] = sqrt( pow(Xab[z][v][0],2)+ pow(Xab[z][v][1],2) + pow(Xab[z][v][2],2) );
    }

  int N1,l;
  double *SUM, *CONST, *VAR;
  N1 = 10;
  SUM = malloc( N1 * sizeof(double));
  CONST = malloc( N1 * sizeof(double));
  VAR = malloc( coordenadas * sizeof(double));
    
  
  for ( z = 0; z < N; z++) {
    for ( v = 0; v < 3; v++)
      dvdt[z][v] = y[z][v+3];
    for ( v = 3; v < 6; v++) {
      for( k = 0; k < N; k++) {
        if ( k != z ) {
          if ( N > 2 ) {
            for ( w = 0; w < N; w++ ) {
              if ( w != z && w!= k ) {
                SUM[2] = SUM[2] + G*M[w]/Rab[k][w];
                SUM[3] = SUM[3] + G*M[w]/Rab[z][w];
                CONST[0] = 0;
                for ( l = 0; l < coordenadas; l++) {
                  CONST[0] = CONST[0] + Xab[z][k][l]*Xab[k][w][l]; // Xab*Xbc
                }
                SUM[4] = SUM[4] + G*M[w]/pow(Rab[k][w],3)*CONST[0];
                SUM[6] = SUM[6] + G*M[w]*Xab[k][w][v-3]/pow(Rab[k][w],3);
                
              }
            }
          }
          SUM[0] = SUM[0] + (G*M[k]*Xab[z][k][v-3])/pow(Rab[z][k],3);
          for ( l = 0; l < coordenadas; l++) {
            CONST[1] = CONST[1] + pow(y[z][l+3],2); // |Va|
            CONST[2] = CONST[2] + y[z][l+3]*y[k][l+3]; // Va*Vb
            CONST[3] = CONST[3] + pow(y[k][l+3],2); // |Vb|
            CONST[4] = CONST[4] + y[k][l+3]*(Xab[z][k][l]/Rab[z][k]); // |Vb*Nab|
            VAR[l] = 4*y[z][l+3] - 3*y[k][l+3]; // (4Va-3Vb)
            CONST[5] = CONST[5] + Xab[z][k][l]*VAR[l]; // Xab*(4Va-3Vb)
          }
          SUM[1] = SUM[1] + (G*M[k]*Xab[z][k][v-3])/pow(Rab[z][k],3)*( (4*G*M[k])/Rab[z][k]  
          + (5*G*M[z])/Rab[z][k] + SUM[2] + 4*SUM[3] - 0.5*SUM[4] - CONST[1] +4*CONST[2] - 2*CONST[3] - (3.0/2.0)*pow(CONST[4],2));
          SUM[5] = SUM[5] + (G*M[k]/Rab[z][k])*SUM[6];
          SUM[7] = SUM[7] + G*M[k]*CONST[5]*Vab[z][k][v-3]/pow(Rab[z][k],3);
        }
      }
      dvdt[z][v] = -SUM[0] + SUM[1]/pow(c,2) - (7.0/(2*pow(c,2)))*SUM[5] + SUM[7]/pow(c,2);
      for ( l = 0; l < N1; l++) {
        SUM[l] = 0;  
        CONST[l] = 0;
      }
    }
  }
  
  free(Xab);
  free(Rab);
  free(SUM);
  free(CONST);
  free(VAR);
  free(Vab);   
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
  movimiento = malloc(cuerpos  *sizeof(double));
  
  
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
//  masas[2] = 7.349e22;

  iniciales[0][0] = 0;
  iniciales[0][1] = 0;
  iniciales[0][2] = 0;
  iniciales[0][3] = 0;
  iniciales[0][5] = 0;
  iniciales[1][1] = 0;
  iniciales[1][2] = 0;
  iniciales[1][3] = 0;
  iniciales[1][5] = 0;
  iniciales[1][4] = c/4;
  iniciales[0][4] = -(masas[1]/masas[0])*iniciales[1][4];
  iniciales[1][0] = G*masas[0]/pow(iniciales[1][4],2);

// 	Tierra, Sol, Luna
/*
  iniciales[2][1] = 0;
  iniciales[2][2] = 0;
  iniciales[2][5] = 0;
  iniciales[1][4] = 28.73e3; //Vy2
  iniciales[2][4] = iniciales[1][4] + 1038.6;
  iniciales[0][4] = -(masas[1]/masas[0])*iniciales[1][4];  // Vy1 para conservar momentum
  iniciales[1][0] = G*masas[0]/pow(iniciales[1][4],2); // X2 dependiente de velocidad para Orbita circular ( E = 0 )
  iniciales[2][0] = G*(masas[1])/pow(iniciales[2][4]-iniciales[1][4],2)+iniciales[1][0]; // radio de luna
*/

  
  
  
  for (i = 0; i < cuerpos; i++) {
    for ( j = 0; j < variables; j++) {
      movimiento[i][j] = iniciales[i][j];
    }
  }
  
  
  
  // Variables de Control
  double dt, tiempo, Nit, *SUMAS,*AP;
  int print;
  
  SUMAS = malloc( cuerpos * sizeof(double));
  AP = malloc( 10 * sizeof(double));
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
  
//  PostNewtonNcuerposV2(tiempo, iniciales, movimiento, masas,cuerpos);
  i=0;  
  while ( i <= Nit ) {
    if ( i % print == 0 ) {
      fprintf(out,"%.10f	",tiempo);
      for ( j = 0; j < cuerpos; j++) {
        for ( k = 0; k < variables/2; k++) {
          fprintf(out, "%.15f	",movimiento[j][k]);
        }
      }
    fprintf(out,"\n");
    }
    RK4MATRIX(movimiento,iniciales,dt,tiempo,PostNewtonNcuerposV2,cuerpos,masas);
    
    for ( j = 0; j < cuerpos; j++)
      for ( k = 0; k < variables; k++)
        iniciales[j][k] = movimiento[j][k];
        
      
    // Cálculo de Energia
    for ( j = 0; j < cuerpos; j++) {
      for ( k = 0; k < variables; k++) {
        AP[0] = AP[0] + pow(iniciales[j][k+3],2); //Va^2
      }
      SUMAS[0] = SUMAS[0] + G*masas[j]*AP[0];      
    }
    
    i++;
    tiempo += dt;
  }
  
  fclose(out);
  free(masas);
  free(iniciales);
  free(movimiento);
  
}

  
