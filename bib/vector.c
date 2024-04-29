#include <math.h>
#include <stdio.h>
#include <stdlib.h>


struct VECTOR {
    double x;
    double y;
};



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

void reflect_vector(struct VECTOR *v, struct VECTOR *normal,struct VECTOR *ref) {
    double prod = dot(v, normal);
    mult(normal,2 * prod);

    relative(v, normal,ref);
    
    mult(normal,1/(2 * prod));
}
void copiar(struct VECTOR *original,struct VECTOR *copia){
    copia->x = original->x;
    copia->y = original->y;
}