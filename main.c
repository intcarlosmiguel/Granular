#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "bib/simulate.h"
#include <omp.h>

int main(int argc, char *argv[]){
    struct CONFIG config;
    config.colunas = 20;
    config.linhas = 40;
    config.tempo_total = 100.;
    config.dt = 1e-6;

    config.atrito_particulas = atof(argv[1]);
    config.atrito_retas = atof(argv[2]);
    config.angulo = atof(argv[3]);
    int seed = atoi(argv[4]);
    config.rotacao = atof(argv[5]);
    config.gamma = atof(argv[6]);
    config.abertura = atof(argv[7])/10;

    seed += (int) config.atrito_particulas + config.atrito_retas;
    config.atrito_particulas /= 10;
    config.atrito_retas /=10;
    omp_set_num_threads(14);

    int count = 0;

    #pragma omp parallel for
    for ( int i = 0; i < 500; i++){
        
        simulate(config, seed+i,i == 0);

        if((count+1)%100 == 0) printf("Simulação %d concluída\n", count+1);
        count++;
    }
    

}