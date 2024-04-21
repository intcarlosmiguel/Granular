#ifndef PARTICLE_H
#define PARTICLE_H

#include "vector.h"
#include <stdbool.h>

struct particula{

    struct VECTOR posicao;
    struct VECTOR velocidade;
    struct VECTOR Force;

    double angular;
    double aceleracao_angular;

    double massa;
    double raio;
    double atrito;
    double Young;
    double A;
    double gamma;
};

struct reta{
    double a;
    double b;
    double c;
    struct VECTOR inicio;
    struct VECTOR fim;
    double A;
    double atrito;
    double gamma;
};
void correct(struct reta *RETA,struct particula *p,struct VECTOR *DEFORMACAO);

double distance_ponto_reta(struct reta* RETA,struct VECTOR *posicao);
void intersecao_circulo_reta(struct reta *RETA,struct particula *p,double *deformacao,struct VECTOR *DEFORMACAO);
double maximo(double a,double b);
double minimo(double a,double b);
bool entre(struct reta *RETA,struct VECTOR *point);
void reflect_reta(struct VECTOR *point, struct reta *RETA,struct VECTOR *reflect);
bool acima(struct reta* RETA,struct VECTOR *posicao,struct VECTOR *CM);
void init_coef(struct reta* RETA,double x1,double y1,double x2, double y2);

double distance_ponto_ponto(struct VECTOR *posicao1,struct VECTOR *posicao2);
bool force_plano(struct particula *particula,struct reta *RETA,struct VECTOR* CORRECAO);
bool force(struct particula *particula1, struct particula *particula2,struct VECTOR * CORRECAO);
#endif

