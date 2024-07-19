#pragma once
#include "vector.h"
#include "math.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "define.h"

void init_coef(struct reta* RETA,double x1,double y1,double x2, double y2){
    RETA->a = (x2 != x1)? (y2 - y1)/(x2 - x1): -1;
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

void entre(struct reta *RETA,struct particula *p,double* t){

    double B = RETA->fim.x - RETA->inicio.x;
    double D = RETA->fim.y - RETA->inicio.y;

    double A = RETA->inicio.x - p->posicao.x;
    double C = RETA->inicio.y - p->posicao.y;

    double A_ = pow(B,2) + pow(D,2);
    double B_ = 2*(B*A+D*C);
    double C_ = pow(A,2) + pow(C,2) - pow(p->raio,2);
    double discriminante = pow(B_,2) - 4*A_*C_;
    if(discriminante < 0){
        t[0] = -10;
        t[1] = -10;
    }
    else{
        t[0] = (-B_ + sqrt(discriminante))/(A_*2);
        t[1] = (-B_ - sqrt(discriminante))/(A_*2);

    }
}