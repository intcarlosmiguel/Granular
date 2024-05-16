#include <stdio.h>
#include "bib/simulate.h"

void main(){

    int colunas = 16;
    int linhas = 32;
    double tempo_total = 2.;
    
    double angulo = 45;

    double dt = 1e-6;
    simulate(colunas,linhas,tempo_total,angulo,dt);

}