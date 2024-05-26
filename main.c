#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "bib/simulate.h"

int main(int argc, char *argv[]){
    int colunas = 4;
    int linhas = 4;
    double tempo_total = 100.;
    
    double angulo = 45;

    double dt = 1e-6;
    double atrito_particulas = atof(argv[1]);
    double atrito_retas = atof(argv[2]);
    bool rotacao = atof(argv[3]);
    int seed = atoi(argv[4]);
    
    seed += (int) atrito_particulas + atrito_retas;
    atrito_particulas /= 10;
    atrito_retas /=10;
    simulate(colunas,linhas,tempo_total,angulo,dt, atrito_particulas, atrito_retas,seed,rotacao);

}