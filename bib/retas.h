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
bool entre(struct reta *RETA,struct particula *p){

    if(RETA->b == 0) RETA->fim.y - p->posicao.y > 0;

    double alpha = RETA->fim.x - RETA->inicio.x;
    double beta = RETA->fim.y - RETA->inicio.y;

    double a = RETA->inicio.x - p->posicao.x;
    double b = RETA->inicio.y - p->posicao.y;
    double A = pow(alpha,2) + pow(beta,2);
    double B = 2*(alpha*a+beta*b);
    double C = pow(a,2) + pow(b,2) - pow(p->raio,2);
    double discriminante = pow(B,2) - 4*A*C;
    if(discriminante < 0) return false;
    double t1 = (-B + sqrt(discriminante))/A/2;
    double t2 = (-B - sqrt(discriminante))/A/2;
    bool b1 = ((t1 <= 1) &&(t1 >= 0 ));
    bool b2 = ((t2 <= 1) &&(t2 >= 0 ));
    return (b1 && b2);
}