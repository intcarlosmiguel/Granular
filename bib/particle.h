#pragma once
#include "vector.h"
#include "retas.h"
#include "define.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>


void intersecao_circulo_reta(struct reta *RETA,struct particula *p,double *deformacao,struct VECTOR *DEFORMACAO) {
    double alpha = RETA->fim.x - RETA->inicio.x;
    double beta = RETA->fim.y - RETA->inicio.y;

    double a = RETA->inicio.x - p->posicao.x;
    double b = RETA->inicio.y - p->posicao.y;
    double A = pow(alpha,2) + pow(beta,2);
    double B = 2*(alpha*a+beta*b);
    double C = pow(a,2) + pow(b,2) - pow(p->raio,2);
    double discriminante = pow(B,2) - 4*A*C;
    // Verificando as condições de interseção
    if (discriminante >= 0) {
        
        // Duas interseções
        struct VECTOR Central;

        double t1 = (-B + sqrt(discriminante))/(A*2);
        double t2 = (-B - sqrt(discriminante))/(A*2);
        bool b1 = ((t1 <= 1) &&(t1 >= 0 ));
        bool b2 = ((t2 <= 1) &&(t2 >= 0 ));
        if(b1+b2 == 2){
            Central.x =0.5*( RETA->inicio.x + t1*(RETA->fim.x - RETA->inicio.x) +RETA->inicio.x + t2*(RETA->fim.x - RETA->inicio.x ));
            Central.y =0.5*( RETA->inicio.y + t1*(RETA->fim.y - RETA->inicio.y) +RETA->inicio.y + t2*(RETA->fim.y - RETA->inicio.y ));
        }
        else{
            double t = (b1)? t1: t2;
            Central.x = RETA->inicio.x + t*(RETA->fim.x - RETA->inicio.x);
            Central.y = RETA->inicio.y + t*(RETA->fim.y - RETA->inicio.y);
        }
        double distance = distance_ponto_ponto(&Central,&p->posicao);
        //if(check) printf("%f :",distance);
        if(distance >= p->raio){
            *deformacao = fabs(RETA->fim.y - Central.y);
            DEFORMACAO->y =  (RETA->fim.y - Central.y)/ *deformacao;
            DEFORMACAO->x = 0;
            //mult(DEFORMACAO,1./ *deformacao);
        }
        else{
            *deformacao = p->raio - distance;
            relative(&p->posicao,&Central, DEFORMACAO);
            mult(DEFORMACAO,1./(p->raio - *deformacao));
        }
        

        //if(check)printf("%.30f %f\n",*deformacao,distance);
        //if(check) print(DEFORMACAO);
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



double force_plano(struct particula *p,struct reta *RETA,struct VECTOR* FORCE,struct VECTOR* Central){
    //printf("Entrou aqui!\n");
    struct VECTOR NORMAL;
    FORCE->x = 0;
    FORCE->y = 0;
    double deformacao = p->raio - distance_ponto_ponto(Central,&p->posicao);
    relative(&p->posicao,Central, &NORMAL);
    mult(&NORMAL,1./(p->raio - deformacao));
    //printf("%f %f\n",distance_ponto_ponto(Central,&p->posicao),p->raio);
    //exit(0);
    if(deformacao > 0){

        double atrito = (p->atrito < RETA->atrito)? p->atrito : RETA->atrito;
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
            printf("%f\n",deformacao);
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
