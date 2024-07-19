#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "bib/simulate.h"
#include <omp.h>

int main(int argc, char *argv[]){
    int colunas = 20;
    int linhas = 18;
    double tempo_total = 100.;
    
    double dt = 1e-6;
    double atrito_particulas = atof(argv[1]);
    double atrito_retas = atof(argv[2]);
    double angulo = atof(argv[3]);
    int seed = atoi(argv[4]);
    bool rotacao = atof(argv[5]);
    double alpha = 3.;
    
    seed += (int) atrito_particulas + atrito_retas;
    atrito_particulas /= 10;
    atrito_retas /=10;
    omp_set_num_threads(6);
    #pragma omp parallel for
    for ( int i = 0; i < 500; i++){
        simulate(colunas,linhas,tempo_total,angulo,dt, atrito_particulas, atrito_retas,alpha,seed+i,rotacao);
    }
    

}