#pragma once
#include <stdbool.h>
#include <stdlib.h>
#include "particle.h"
#include "define.h"
#include "grid.h"
void integracao(struct particula *p,struct VECTOR *anterior,double dt){
    double valor;
    p->Force.y += -p->massa*GRAVIDADE;
    mult(&p->Force,1/p->massa);

    
    p->velocidade.x += dt*(p->Force.x + p->aceleracao.x);
    p->velocidade.y += dt*(p->Force.y + p->aceleracao.y);
    valor = (p->velocidade.x != 0)? 2*p->posicao.x - anterior->x + pow(dt,2)*(p->Force.x+ p->aceleracao.x) : p->posicao.x + pow(dt,2)*(p->Force.x+ p->aceleracao.x)/2;

    anterior->x = p->posicao.x;
    p->posicao.x = valor;

    valor= (p->velocidade.y != 0)? 2*p->posicao.y - anterior->y + (p->Force.y + p->aceleracao.y)*pow(dt,2) : p->posicao.y + pow(dt,2)*(p->Force.y + p->aceleracao.y)/2;
    
    anterior->y = p->posicao.y;
    p->posicao.y = valor;

    p->angular += (p->aceleracao_angular + p->Force_Rot/p->Inertia)*dt;

    p->aceleracao_angular = p->Force_Rot/p->Inertia;

    p->aceleracao.x = p->Force.x;
    p->aceleracao.y = p->Force.y;
    
    p->Force.x = 0;
    p->Force.y = 0;
    p->Force_Rot = 0;
}

struct particula* calc_force_par(int site,struct particula* p,int N,struct GRID *grid,bool rotacao){

    int id,viz,celula_vizinha,j,vizinho;
    struct VECTOR FORCE;
    FORCE.x = 0;
    FORCE.y = 0;
    double force_rotacao;

    id = grid->ids[site];

    if(id > 0){
        for ( viz = 0; viz < 9; viz++){

            celula_vizinha = get_vizinho(grid,id,viz);

            if((celula_vizinha > 0)){
                
                if(grid->celulas[celula_vizinha][0] > 0){

                    for ( j = 1; j <= grid->celulas[celula_vizinha][0]; j++){

                        vizinho = grid->celulas[celula_vizinha][j];
                        if(site >= vizinho) continue;
                        if(p[vizinho].posicao.y < 0.01) continue;
                        if(distance_ponto_ponto(&p[site].posicao,&p[vizinho].posicao)< p[site].raio + p[vizinho].raio){
                            force_rotacao = force(&p[site],&p[vizinho],&FORCE);
                            p[site].Force.x += FORCE.x;
                            p[site].Force.y += FORCE.y;
                            if(rotacao)p[site].Force_Rot += force_rotacao*p[site].raio;

                            p[vizinho].Force.x -= FORCE.x;
                            p[vizinho].Force.y -= FORCE.y;
                            if(rotacao)p[vizinho].Force_Rot -= force_rotacao*p[vizinho].raio;

                            FORCE.x = 0;
                            FORCE.y = 0;
                            force_rotacao = 0;
                        }
                    }
                    
                }
            }
        }
        
    }
    return p;
}
struct particula* calc_force_reta(struct particula* particulas,int site,struct reta *retas,int N,bool rotacao,int n_retas){
    int j;
    double force_rotacao;
    struct VECTOR FORCE;
    struct VECTOR Central;
    FORCE.x = 0;
    FORCE.y = 0;
    for ( j = 0; j < n_retas; j++){

        //entre(&retas[j],&particulas[site],t);
        double a = particulas[site].posicao.x - retas[j].inicio.x;
        double b = particulas[site].posicao.y - retas[j].inicio.y;
        double c = retas[j].fim.x - retas[j].inicio.x;
        double d = retas[j].fim.y - retas[j].inicio.y;

        double ratio = (a*c+b*d)/(c*c+d*d);

        if(ratio < 0.0){
            Central.x = retas[j].inicio.x;
            Central.y = retas[j].inicio.y;
        }
        else{
            if(ratio > 1.0){
                Central.x = retas[j].fim.x;
                Central.y = retas[j].fim.y;
            }
            else{
                Central.x = c*ratio + retas[j].inicio.x;
                Central.y = d*ratio + retas[j].inicio.y;
            }
        }

        force_rotacao = force_plano(&particulas[site],&retas[j],&FORCE,&Central);

        particulas[site].Force.x += FORCE.x;
        particulas[site].Force.y += FORCE.y;

        FORCE.x = 0;
        FORCE.y = 0;
        
        if(rotacao)particulas[site].Force_Rot += force_rotacao*particulas[site].raio;
        force_rotacao = 0;
    }
    return particulas;
}