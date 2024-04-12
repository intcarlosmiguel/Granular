#include "vector.h"
#include "math.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

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

double distance_ponto_ponto(struct VECTOR *posicao1,struct VECTOR *posicao2){

    return sqrt(pow(posicao1->x - posicao2->x,2)+pow(posicao1->y - posicao2->y,2)); 
}

struct VECTOR find_tangente(struct VECTOR *velocidade1,struct VECTOR *velocidade2,struct VECTOR *NORMAL){

    struct VECTOR rot1,rot2;

    double** M = (double**) malloc(2*sizeof(double*));
    M[0] = (double*)calloc(2,sizeof(double));
    M[1] = (double*)calloc(2,sizeof(double));

    rotate(NORMAL,&rot1,0);
    rotate(NORMAL,&rot2,1);

    M[0][0] = dot(&rot1,velocidade1);
    M[0][1] = dot(&rot2,velocidade1);
    M[1][0] = dot(&rot1,velocidade2);
    M[1][1] = dot(&rot2,velocidade2);
    int a = (norma(velocidade1) > norma(velocidade2))? 0 : 1;
    bool b = (M[a][0] < M[a][1]);
    free(M[0]);
    free(M[1]);
    free(M);
    if(b) return rot2;
    else return rot1;
}

void intersecao_circulo_reta(struct reta *RETA,struct particula *p,double *deformacao,struct VECTOR *DEFORMACAO) {
    double discriminante;
    double x1, x2, y1, y2;

    // Determinando os coeficientes da equação quadrática resultante da substituição da reta na equação do círculo

    double a = RETA->a * RETA->a + RETA->b * RETA->b;

    double b = 2 * (RETA->a * (RETA->c - p->posicao.x) + RETA->b * (RETA->c - p->posicao.y));
    double c = p->posicao.x * p->posicao.x + p->posicao.y * p->posicao.y + RETA->c * RETA->c - 2 * (p->posicao.x * RETA->a + p->posicao.y * RETA->b) - p->raio * p->raio;

    // Calculando o discriminante
    discriminante = b * b - 4 * a * c;

    // Verificando as condições de interseção
    if (discriminante > 0) {
        // Duas interseções
        struct VECTOR ponto1;
        struct VECTOR ponto2;
        struct VECTOR Central;
        
        ponto1.x = (-b + sqrt(discriminante)) / (2 * a);
        ponto2.x = (-b - sqrt(discriminante)) / (2 * a);
        ponto1.y = (-RETA->a * x1 - RETA->c) / RETA->b;
        ponto2.y = (-RETA->a * x2 - RETA->c) / RETA->b;

        Central.x = 0.5*(ponto1.x+ponto2.x);
        Central.y = 0.5*(ponto1.y+ponto2.y);

        *deformacao = p->raio - distance_ponto_ponto(&Central,&p->posicao);
        relative(&p->posicao,&Central, DEFORMACAO);
    } else{
        printf("Deu problema!");
        exit(0);
    }
}

void force_plano(struct particula *particula,struct reta *RETA){
    struct VECTOR NORMAL,FORCE;
    FORCE.x = 0;
    FORCE.y = 0;
    double deformacao;
    double atrito = RETA->atrito;

    intersecao_circulo_reta(RETA,particula,&deformacao,&NORMAL);

    mult(&NORMAL,1/(particula->raio - deformacao));

    double dv = -dot(&NORMAL,&particula->velocidade);
    double forca_normal = 4/3*sqrt(particula->raio)*particula->Young*sqrt(deformacao)*(deformacao + particula->A*dv);

    double velocidade_tangencial = dv + particula->raio*particula->angular;

    if(forca_normal < 0) forca_normal = 0;

    double forca_tangencial = -particula->gamma *velocidade_tangencial;
    if(forca_tangencial < -atrito*forca_normal) forca_tangencial = -atrito*forca_normal;
    if(forca_tangencial > atrito*forca_normal) forca_tangencial = atrito*forca_normal;
    struct VECTOR vel;
    vel.x = 0;
    vel.y = 0;
    struct VECTOR ROTATE = find_tangente(&particula->velocidade, &vel,&NORMAL);
        
    mult(&ROTATE,forca_tangencial);
    mult(&NORMAL,forca_normal);
    sum(&NORMAL,&ROTATE,&FORCE);
    sum(&FORCE,&particula->Force,&FORCE);
}

void force(struct particula *particula1, struct particula *particula2){

    struct VECTOR NORMAL,FORCE;
    FORCE.x = 0;
    FORCE.y = 0;

    relative(&particula1->posicao,&particula2->posicao, &NORMAL);
    double dist = norma(&NORMAL);
    double deformacao = particula1->raio + particula2->raio - dist;
    
    if (deformacao > 0){

        double Young = (particula1->Young*particula2->Young)/(particula1->Young+particula2->Young);

        double Raio_effetivo = (particula2->raio < 1000*particula1->raio)? (particula1->raio*particula2->raio)/(particula1->raio+particula2->raio) :particula1->raio;

        double Massa =(particula2->massa < 1000*particula1->massa)?(particula1->massa*particula2->massa)/(particula1->massa+particula2->massa) : particula1->massa;

        double A = particula1->A + particula2->A;

        double atrito = (particula1->atrito > particula2->atrito)? particula1->atrito : particula2->atrito;

        double gamma = (particula1->gamma > particula2->gamma)? particula1->gamma : particula2->gamma;

        struct VECTOR VELOCIDADE,TANGENCIAL;
        
        relative(&particula1->velocidade,&particula2->velocidade, &VELOCIDADE);

        mult(&NORMAL,1/dist);

        double dv = -dot(&NORMAL,&VELOCIDADE);

        double velocidade_tangencial = dot(&VELOCIDADE,&NORMAL) + particula1->raio*particula1->angular+particula2->raio*particula2->angular;

        double forca_normal = 4/3*sqrt(Raio_effetivo)*Young*sqrt(deformacao)*(deformacao + A*dv);
        if(forca_normal < 0) forca_normal = 0;

        double forca_tangencial = -gamma*velocidade_tangencial;
        if(forca_tangencial < -atrito*forca_normal) forca_tangencial = -atrito*forca_normal;
        if(forca_tangencial > atrito*forca_normal) forca_tangencial = atrito*forca_normal;

        struct VECTOR ROTATE = find_tangente(&particula1->velocidade, &particula2->velocidade,&NORMAL);
        
        mult(&ROTATE,forca_tangencial);
        mult(&NORMAL,forca_normal);
        sum(&NORMAL,&ROTATE,&FORCE);
    }
    sum(&FORCE,&particula1->Force,&FORCE);
    sum(&FORCE,&particula2->Force,&FORCE);
}

double distance_ponto_reta(struct reta* RETA,struct VECTOR *posicao){
    return abs(RETA->a*posicao->x + RETA->b*posicao->y + RETA->c )/sqrt(RETA->a*RETA->a + RETA->b*RETA->b);
}   



double maximo(double a,double b){
    if(a > b) return a;
    else return b;
}
double minimo(double a,double b){
    if(a > b) return b;
    else return a;
}

bool entre(struct reta *RETA,struct VECTOR *point){
    bool a = ((point->x < maximo(RETA->inicio.x,RETA->fim.x)) && (point->x > minimo(RETA->inicio.x,RETA->fim.x)));
    bool b = ((point->y < maximo(RETA->inicio.y,RETA->fim.y)) && (point->y > minimo(RETA->inicio.y,RETA->fim.y)));
    return a && b;
}

void reflect_reta(struct VECTOR *point, struct reta *RETA,struct VECTOR *reflect) {
    double distance = distance_ponto_reta(RETA,point);
    double normalMagnitude = sqrt(RETA->a * RETA->a + RETA->b * RETA->b);
    
    // Calcula o ponto reflexo usando a fórmula: P' = P - 2 * distância * normal da reta
    reflect->x = point->x - 2 * distance * (RETA->a / normalMagnitude);
    reflect->y = point->y - 2 * distance * (RETA->a / normalMagnitude);
}

bool acima(struct reta* RETA,struct VECTOR *posicao){
    double plano = RETA->a*posicao->x + RETA->b*posicao->y + RETA->c;
    return plano > 0;
}

void init_coef(struct reta* RETA,double x1,double y1,double x2, double y2){
    RETA->a = (y2 - y1)/(x2 - x1);
    RETA->b = -1;
    RETA->c = y1 - RETA->a*x1;
    RETA->inicio.x = x1;
    RETA->inicio.y = y1;
    RETA->fim.x = x2;
    RETA->fim.y = y2;
}

