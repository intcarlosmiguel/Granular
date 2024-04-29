#ifndef VECTOR_H
#define VECTOR_H

struct VECTOR {
    double x;
    double y;
};

void print(struct VECTOR *v);
double dot(struct VECTOR *v1,struct VECTOR *v2);
void relative(struct VECTOR *v1,struct VECTOR *v2, struct VECTOR *v);
void sum(struct VECTOR *v1,struct VECTOR *v2, struct VECTOR *resultado);
double norma(struct VECTOR *v);
void mult(struct VECTOR *v,double c);
void rotate(struct VECTOR *v,struct VECTOR *rot,int check);
void reflect_vector(struct VECTOR *v, struct VECTOR *normal,struct VECTOR *ref);
void copiar(struct VECTOR *original,struct VECTOR *copia);
void adicionar(struct VECTOR *v,double c);
#endif 
