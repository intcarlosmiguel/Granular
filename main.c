#include <stdio.h>
#include "bib/simulate.h"

void main(){

    int colunas = 2;
    int linhas = 2;

    double tempo_total = 1000;
    
    double angulo = 45;

    simulate(colunas,linhas,tempo_total,angulo);

}