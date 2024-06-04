#pragma once
#include "vector.h"
#include "retas.h"
#include "define.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>


void intersecao_circulo_reta(struct reta *RETA,struct particula *p,double *deformacao,struct VECTOR *DEFORMACAO) {
    double discriminante;

    // Calculando o discriminante
    double B = 2*(RETA->a*(RETA->c - p->posicao.y) - p->posicao.x)/(1+pow(RETA->a,2));
    double C = (pow(RETA->c - p->posicao.y,2) + pow(p->posicao.x,2) - pow(p->raio,2) )/(1+pow(RETA->a,2));
    discriminante = pow(B,2) - 4*C;
    if(RETA->b ==0){
        discriminante = sqrt(pow(p->raio,2) - pow(-RETA->c - p->posicao.x,2));
        //printf("Teste: %f\n",p->posicao.y + sqrt(pow(p->raio,2) - pow(-RETA->c - p->posicao.x,2)));
    }
    // Verificando as condições de interseção
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
        if(check)print(&Central);
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



double force_plano(struct particula *p,struct reta *RETA,struct VECTOR* FORCE){
    //printf("Entrou aqui!\n");
    struct VECTOR NORMAL;
    FORCE->x = 0;
    FORCE->y = 0;
    double deformacao;
    intersecao_circulo_reta(RETA,p,&deformacao,&NORMAL);
    //if(check)printf("%f\n",deformacao);
    if(deformacao > 0){

        double atrito = (p->atrito < RETA->atrito)? p->atrito : RETA->atrito;
        mult(&NORMAL,1./(p->raio - deformacao));
        double dv = -dot(&NORMAL,&p->velocidade);
        double forca_normal = 4./3.*sqrt(p->raio)*p->Young*sqrt(deformacao)*(deformacao + 0.5*(p->A+RETA->A)*dv)/2;

        double velocidade_tangencial = -p->velocidade.x*NORMAL.y + NORMAL.x*p->velocidade.y + p->raio*p->angular;
        if(forca_normal < 0) forca_normal = 0;

        double forca_tangencial = -p->gamma *velocidade_tangencial;
        if(forca_tangencial < -atrito*forca_normal) forca_tangencial = -atrito*forca_normal;
        if(forca_tangencial > atrito*forca_normal) forca_tangencial = atrito*forca_normal;

        if (isnan(forca_normal) || isnan(forca_tangencial)) {
            printf("Erro: forca_normal ou forca_tangencial é NaN na função force_plano\n");
            exit(1);
        }
        
        struct VECTOR TANGENCIAL;
        TANGENCIAL.x = -NORMAL.y;
        TANGENCIAL.y = NORMAL.x;
        FORCE->x = NORMAL.x*forca_normal + TANGENCIAL.x*forca_tangencial;
        FORCE->y = NORMAL.y*forca_normal + TANGENCIAL.y*forca_tangencial;
        return forca_tangencial;
    }
    return 0;
}

double force(struct particula *particula1, struct particula *particula2,struct VECTOR * FORCE){

    struct VECTOR NORMAL;
    FORCE->x = 0;
    FORCE->y = 0;

    relative(&particula1->posicao,&particula2->posicao, &NORMAL);
    double dist = norma(&NORMAL);
    double deformacao = particula1->raio + particula2->raio - dist;
    if (deformacao > 0){
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

        if (isnan(forca_normal) || isnan(forca_tangencial)) {
            printf("Erro: forca_normal ou forca_tangencial é NaN na função force\n");
            exit(1);
        }


        struct VECTOR TANGENCIAL;
        TANGENCIAL.x = -NORMAL.y;
        TANGENCIAL.y = NORMAL.x;
        FORCE->x = NORMAL.x*forca_normal + TANGENCIAL.x*forca_tangencial;
        FORCE->y = NORMAL.y*forca_normal + TANGENCIAL.y*forca_tangencial;
        return forca_tangencial;
    }
    return 0;
}
