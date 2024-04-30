#include "vector.h"
#include "math.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

struct particula{

    struct VECTOR posicao;
    struct VECTOR velocidade;
    struct VECTOR Force;
    struct VECTOR aceleracao;

    double angular;
    double aceleracao_angular;

    double massa;
    double raio;
    double atrito;
    double Young;
    double A;
    double gamma;
    //double densidade;
    //double poisson;

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

void print_vector(struct VECTOR * vetor){
    printf("%f - %f\n",vetor->x,vetor->y);
}

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

    // Calculando o discriminante
    double B = 2*(RETA->a*(RETA->c - p->posicao.y) - p->posicao.x)/(1+pow(RETA->a,2));
    double C = (pow(RETA->c - p->posicao.y,2) + pow(p->posicao.x,2) - pow(p->raio,2) )/(1+pow(RETA->a,2));
    discriminante = pow(B,2) - 4*C;
    if(RETA->b ==0){
        discriminante = sqrt(pow(p->raio,2) - pow(-RETA->c - p->posicao.x,2));
        //printf("Teste: %f\n",p->posicao.y + sqrt(pow(p->raio,2) - pow(-RETA->c - p->posicao.x,2)));
    }
    // Verificando as condições de interseção
    if(fabs(discriminante) < 0.000001) discriminante = 0;
    if (discriminante >= 0) {
        // Duas interseções
        struct VECTOR ponto1;
        struct VECTOR ponto2;
        struct VECTOR Central;
        
        ponto1.x = (-B + sqrt(discriminante))/2;
        ponto2.x = (-B - sqrt(discriminante))/2;
        ponto1.y = RETA->a * ponto1.x + RETA->c;
        ponto2.y = RETA->a * ponto2.x + RETA->c;

        if(RETA->b == 0){
            ponto1.x = -RETA->c;
            ponto2.x = -RETA->c;
            ponto1.y = p->posicao.y + sqrt(pow(p->raio,2) - pow(-RETA->c - p->posicao.x,2));
            ponto2.y = p->posicao.y - sqrt(pow(p->raio,2) - pow(-RETA->c - p->posicao.x,2));
        }
        Central.x = 0.5*(ponto1.x+ponto2.x);
        Central.y = 0.5*(ponto1.y+ponto2.y);
        if(distance_ponto_ponto(&Central,&p->posicao) > p->raio){
            printf("%f\n",Central.x*RETA->a + Central.y*RETA->b + RETA->c);
            printf("Erro:  distancia maior que o raio!\n");
            exit(0);
        }
        double distance = distance_ponto_ponto(&Central,&p->posicao);
        *deformacao = p->raio - distance;

        relative(&p->posicao,&Central, DEFORMACAO);
    } else{
        printf("Deu problema na Colisão com reta!\n");
        printf("Particula: (%f,%f)\n",p->posicao.x,p->posicao.y);
        printf("RETA: (%f,%f)\n",RETA->inicio.x,RETA->inicio.y);
        printf("RETA: (%f,%f)\n",RETA->fim.x,RETA->fim.y);
        printf("RETA: (%f,%f,%f)\n",RETA->a,RETA->b,RETA->c);
        printf("%.20f\n",discriminante);
        printf("%f\n", B*B - 4*C);
        exit(0);
    }
}

bool force_plano(struct particula *particula,struct reta *RETA,struct VECTOR* FORCE){
    struct VECTOR NORMAL;
    FORCE->x = 0;
    FORCE->y = 0;
    double deformacao;
    double atrito = RETA->atrito;

    intersecao_circulo_reta(RETA,particula,&deformacao,&NORMAL);

    mult(&NORMAL,1./(particula->raio - deformacao));
    
    double dv = -dot(&NORMAL,&particula->velocidade);
    double forca_normal = 4/3*sqrt(particula->raio)*particula->Young*sqrt(deformacao)*(deformacao + 0.5*(particula->A+RETA->A)*dv);

    double velocidade_tangencial = -particula->velocidade.x*NORMAL.y + NORMAL.x*particula->velocidade.y + particula->raio*particula->angular;
    if(forca_normal < 0) forca_normal = 0;

    double forca_tangencial = -particula->gamma *velocidade_tangencial;
    if(forca_tangencial < -atrito*forca_normal) forca_tangencial = -atrito*forca_normal;
    if(forca_tangencial > atrito*forca_normal) forca_tangencial = atrito*forca_normal;

    struct VECTOR vel;
    vel.x = 0;
    vel.y = 0;
    struct VECTOR TANGENCIAL;
    TANGENCIAL.x = -NORMAL.y;
    TANGENCIAL.y = NORMAL.x;
    mult(&TANGENCIAL,forca_tangencial);
    mult(&NORMAL,forca_normal);
    sum(&NORMAL,&TANGENCIAL,FORCE);
    return false;
}

bool force(struct particula *particula1, struct particula *particula2,struct VECTOR * FORCE){

    struct VECTOR NORMAL;
    FORCE->x = 0;
    FORCE->y = 0;

    relative(&particula1->posicao,&particula2->posicao, &NORMAL);
    double dist = norma(&NORMAL);
    double deformacao = particula1->raio + particula2->raio - dist;
    bool done = true;
    if (deformacao > 0.001){
        done = false;
        double Young = (particula1->Young*particula2->Young)/(particula1->Young+particula2->Young);

        double Raio_effetivo =(particula1->raio*particula2->raio)/(particula1->raio+particula2->raio) ;

        double A = 0.5*(particula1->A + particula2->A);

        double atrito = (particula1->atrito > particula2->atrito)? particula1->atrito : particula2->atrito;

        double gamma = (particula1->gamma > particula2->gamma)? particula1->gamma : particula2->gamma;

        struct VECTOR VELOCIDADE;
        
        relative(&particula1->velocidade,&particula2->velocidade, &VELOCIDADE);
        mult(&NORMAL,1/dist);
        double dv = -dot(&NORMAL,&VELOCIDADE);

        double velocidade_tangencial = -VELOCIDADE.x*NORMAL.y + NORMAL.y*VELOCIDADE.x  + particula1->raio*particula1->angular+particula2->raio*particula2->angular;

        double forca_normal = 4/3*sqrt(Raio_effetivo)*Young*sqrt(deformacao)*(deformacao + A*dv);
        if(forca_normal < 0) forca_normal = 0;

        double forca_tangencial = -gamma*velocidade_tangencial;
        if(forca_tangencial < -atrito*forca_normal) forca_tangencial = -atrito*forca_normal;
        if(forca_tangencial > atrito*forca_normal) forca_tangencial = atrito*forca_normal;
        struct VECTOR TANGENCIAL;
        TANGENCIAL.x = -NORMAL.y;
        TANGENCIAL.y = NORMAL.x;
        FORCE->x = NORMAL.x*forca_normal + TANGENCIAL.x*forca_tangencial;
        FORCE->y = NORMAL.y*forca_normal + TANGENCIAL.y*forca_tangencial;
    }
    return done;

}

double distance_ponto_reta(struct reta* RETA,struct VECTOR *posicao){
    return fabs(RETA->a*posicao->x + RETA->b*posicao->y + RETA->c )/sqrt(pow(RETA->a,2) + pow(RETA->b,2));
}   


bool entre(struct reta *RETA,struct VECTOR *point){
    struct VECTOR A, B, P, AP, AB;
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

void reflect_reta(struct VECTOR *point, struct reta *RETA,struct VECTOR *reflect) {
    struct VECTOR NORMAL;
    NORMAL.x = RETA->a;
    NORMAL.y = RETA->b;
    mult(&NORMAL,1/norma(&NORMAL));
    if(RETA->fim.x == RETA->inicio.x){
        reflect->x = 2*fabs(RETA->inicio.x - reflect->x)+reflect->x;
        reflect->y = reflect->y;
    }
    else{
        double d = dot(point,&NORMAL) - RETA->c;
        mult(&NORMAL,-2*d);
        sum(point,&NORMAL,reflect);
    }
    //double distance = distance_ponto_reta(RETA,point);
    //double normalMagnitude = sqrt(RETA->a * RETA->a + RETA->b * RETA->b);
    
    // Calcula o ponto reflexo usando a fórmula: P' = P - 2 * distância * normal da reta
    
    
}

bool acima(struct reta* RETA,struct VECTOR *posicao,struct VECTOR *CM){
    double dir_CM = (RETA->a != 0)? RETA->a*CM->x + RETA->b*CM->y + RETA->c : CM->x - RETA->c;
    double dir_particula = (RETA->a != 0)? RETA->a*posicao->x + RETA->b*posicao->y + RETA->c : posicao->x - RETA->c;
    return dir_CM/dir_particula > 0;
}

void init_coef(struct reta* RETA,double x1,double y1,double x2, double y2){
    RETA->a = (x2 != x1)? (y2 - y1)/(x2 - x1): 1;
    RETA->b = (x2 != x1)? -1: 0;
    RETA->c = (x2 != x1)? y1 - RETA->a*x1:-1.0*RETA->a*x1 ;
    RETA->inicio.x = x1;
    RETA->inicio.y = y1;
    RETA->fim.x = x2;
    RETA->fim.y = y2;
}

