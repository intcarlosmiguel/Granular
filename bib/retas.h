#pragma once
#include "vector.h"
#include "math.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "define.h"

void init_coef(struct reta* RETA,double x1,double y1,double x2, double y2){
    RETA->a = (x2 != x1)? (y2 - y1)/(x2 - x1): 1;
    RETA->b = (x2 != x1)? -1: 0;
    RETA->c = (x2 != x1)? y1 - RETA->a*x1:-1.0*RETA->a*x1 ;
    RETA->inicio.x = x1;
    RETA->inicio.y = y1;
    RETA->fim.x = x2;
    RETA->fim.y = y2;
}

double distance_ponto_reta(struct reta* RETA,struct VECTOR *posicao){
    return fabs(RETA->a*posicao->x + RETA->b*posicao->y + RETA->c )/sqrt(pow(RETA->a,2) + pow(RETA->b,2));
}
bool entre(struct reta *RETA,struct VECTOR *point){
    struct VECTOR AP, AB;
    double t, denom;

    // Calcula AB = B - A
    relative(&RETA->inicio, &RETA->fim, &AB);

    // Calcula AP = P - A
    relative(point, &RETA->fim, &AP);

    // Calcula o denominador (AB . AB)
    denom = dot(&AB, &AB);
    if (denom == 0) {
        printf("Erro: Divisão por zero ao calcular t, verifique os pontos A e B\n");
        exit(0);
    }

    // Calcula o numerador (AP . AB)
    double numer = dot(&AP, &AB);

    // Calcula t = (AP . AB) / (AB . AB)
    t = numer / denom;

    // Imprime o resultado
    return ((t >= 0) && (t<=1));
}