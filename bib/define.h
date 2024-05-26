#pragma once
const double GRAVIDADE = 9.80665;
double PI;

struct GRID{
    int n_grids;
    int** celulas;
    int* ids;
    int nx;
    int ny;
    double largura;
    double x0,y0;
    /* data */
};


struct VECTOR {
    double x;
    double y;
};

struct particula{

    struct VECTOR posicao;
    struct VECTOR velocidade;
    struct VECTOR Force;
    struct VECTOR aceleracao;

    double angular;
    double aceleracao_angular;

    double massa;
    double Inertia;
    double raio;
    double atrito;
    double Young;
    double A;
    double gamma;
    double Force_Rot;

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