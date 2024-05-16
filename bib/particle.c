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
            printf("Erro:  distancia maior que o raio!\n");
            printf("%f\n",Central.x*RETA->a + Central.y*RETA->b + RETA->c);
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

double distance_ponto_reta(struct reta* RETA,struct VECTOR *posicao){
    return fabs(RETA->a*posicao->x + RETA->b*posicao->y + RETA->c )/sqrt(pow(RETA->a,2) + pow(RETA->b,2));
}   

void force_plano(struct particula *p,struct reta *RETA,struct VECTOR* FORCE){
    //printf("Entrou aqui!\n");
    struct VECTOR NORMAL;
    FORCE->x = 0;
    FORCE->y = 0;
    double deformacao,d;
    double atrito = RETA->atrito;
    intersecao_circulo_reta(RETA,p,&deformacao,&NORMAL);
    d = p->raio - distance_ponto_reta(RETA,&p->posicao);
    //printf("%f %f\n",d,deformacao);
    if(deformacao > 1e-5){

        double atrito = (p->atrito < RETA->atrito)? p->atrito : RETA->atrito;
        mult(&NORMAL,1./(p->raio - deformacao));
        double dv = -dot(&NORMAL,&p->velocidade);
        //printf("Reta: %e %f\n",deformacao,dv);
        //print(&p->posicao);
        double forca_normal = 4./3.*sqrt(p->raio)*p->Young*sqrt(deformacao)*(deformacao + 0.5*(p->A+RETA->A)*dv)/2;

        double velocidade_tangencial = -p->velocidade.x*NORMAL.y + NORMAL.x*p->velocidade.y + p->raio*p->angular;
        if(forca_normal < 0) forca_normal = 0;

        double forca_tangencial = -p->gamma *velocidade_tangencial;
        if(forca_tangencial < -atrito*forca_normal) forca_tangencial = -atrito*forca_normal;
        if(forca_tangencial > atrito*forca_normal) forca_tangencial = atrito*forca_normal;
        
        struct VECTOR TANGENCIAL;
        TANGENCIAL.x = -NORMAL.y;
        TANGENCIAL.y = NORMAL.x;
        FORCE->x = NORMAL.x*forca_normal + TANGENCIAL.x*forca_tangencial;
        FORCE->y = NORMAL.y*forca_normal + TANGENCIAL.y*forca_tangencial;
    }
}

void force(struct particula *particula1, struct particula *particula2,struct VECTOR * FORCE){

    struct VECTOR NORMAL;
    FORCE->x = 0;
    FORCE->y = 0;

    relative(&particula1->posicao,&particula2->posicao, &NORMAL);
    double dist = norma(&NORMAL);
    double deformacao = particula1->raio + particula2->raio - dist;
    if (deformacao > 1e-6){
        double Young = (particula1->Young*particula2->Young)/(particula1->Young+particula2->Young);

        double Raio_effetivo =(particula1->raio*particula2->raio)/(particula1->raio+particula2->raio) ;

        double A = 0.5*(particula1->A + particula2->A);

        double atrito = (particula1->atrito < particula2->atrito)? particula1->atrito : particula2->atrito;

        double gamma = (particula1->gamma < particula2->gamma)? particula1->gamma : particula2->gamma;

        struct VECTOR VELOCIDADE;
        
        relative(&particula1->velocidade,&particula2->velocidade, &VELOCIDADE);
        mult(&NORMAL,1./dist);
        double dv = -dot(&NORMAL,&VELOCIDADE);

        double velocidade_tangencial = -VELOCIDADE.x*NORMAL.y + NORMAL.x*VELOCIDADE.y  + particula1->raio*particula1->angular - particula2->raio*particula2->angular;

        double forca_normal = 4.0/3.0*sqrt(Raio_effetivo)*Young*sqrt(deformacao)*(deformacao + A*dv);
        if(forca_normal < 0) forca_normal = 0;

        double forca_tangencial = -gamma*velocidade_tangencial;

        if(forca_tangencial < -atrito*forca_normal) forca_tangencial = -atrito*forca_normal;
        if(forca_tangencial > atrito*forca_normal) forca_tangencial = atrito*forca_normal;

        struct VECTOR TANGENCIAL;
        TANGENCIAL.x = -NORMAL.y;
        TANGENCIAL.y = NORMAL.x;
        FORCE->x = NORMAL.x*forca_normal + TANGENCIAL.x*forca_tangencial;
        FORCE->y = NORMAL.y*forca_normal + TANGENCIAL.y*forca_tangencial;
        if(FORCE->x != 0){
            /* printf("Forca: %.20f %e %e\n",forca_normal,A,dv);
            print(FORCE); */
        }
    }
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

