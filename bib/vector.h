#pragma once
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "define.h"


void print(struct VECTOR *v){
    printf("(%.20f,%.20f)\n",v->x,v->y);
}
double dot(struct VECTOR *v1,struct VECTOR *v2){
    return v1->x*v2->x + v1->y*v2->y;
}

void relative(struct VECTOR *v1,struct VECTOR *v2, struct VECTOR *resultado){
    resultado->x = v1->x - v2->x;
    resultado->y = v1->y - v2->y;
}

void sum(struct VECTOR *v1,struct VECTOR *v2, struct VECTOR *resultado){
    resultado->x = v1->x + v2->x;
    resultado->y = v1->y + v2->y;
}

double norma(struct VECTOR *v){
    return sqrt(pow(v->x,2) + pow(v->y,2));
}

void adicionar(struct VECTOR *v,double c){
    v->x = v->x+c;
    v->y = v->y+c;
}

void mult(struct VECTOR *v,double c){
    v->x = c*v->x;
    v->y = c*v->y;
}

void rotate(struct VECTOR *v,struct VECTOR *rot,int check){
    switch (check){
    case 0:{
        rot->x = -v->x;
        rot->y = v->y;
        break;
    }
    case 1:{
        rot->x = v->x;
        rot->y = -v->y;
        break;
    }
    default:
        printf("Caso não encontrado!");
        break;
    }
}
void copiar(struct VECTOR *original,struct VECTOR *copia){
    copia->x = original->x;
    copia->y = original->y;
}


double distance_ponto_ponto(struct VECTOR *posicao1,struct VECTOR *posicao2){
    return sqrt(pow(posicao1->x - posicao2->x,2)+pow(posicao1->y - posicao2->y,2)); 
}