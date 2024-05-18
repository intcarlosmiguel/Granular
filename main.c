#include <stdio.h>
#include "bib/simulate.h"

void main(){

    int colunas = 4;
    int linhas = 4;
    double tempo_total = 2.;
    
    double angulo = 45;

    double dt = 1e-6;
    simulate(colunas,linhas,tempo_total,angulo,dt);

}