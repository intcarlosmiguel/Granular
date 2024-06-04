#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "bib/simulate.h"

int main(int argc, char *argv[]){
    int colunas = 20;
    int linhas = 23;
    double tempo_total = 100.;
    

    double dt = 1e-6;
    double atrito_particulas = atof(argv[1]);
    double atrito_retas = atof(argv[2]);
    bool rotacao = atof(argv[3]);
    double angulo = atof(argv[4]);
    int seed = atoi(argv[5]);
    double alpha = 2.5;
    
    seed += (int) atrito_particulas + atrito_retas;
    atrito_particulas /= 10;
    atrito_retas /=10;
    simulate(colunas,linhas,tempo_total,angulo,dt, atrito_particulas, atrito_retas,alpha,1562,rotacao);

}