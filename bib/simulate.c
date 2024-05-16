#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "particle.h"

const double GRAVIDADE = 9.80665;

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

    p->aceleracao.x = p->Force.x;
    p->aceleracao.y = p->Force.y;
    
    p->Force.x = 0;
    p->Force.y = 0;
}

struct particula* corrige_ponto(struct particula* particulas,int N,struct VECTOR* CORRECOES){

    int i,j;
    bool done = true;
    struct VECTOR FORCE;
    FORCE.x = 0;
    FORCE.y = 0;
    for ( i = 0; i < N; i++){
        
        if(particulas[i].posicao.y < 0.) continue;
        for ( j = i+1; j < N; j++){
            if(particulas[j].posicao.y < 0.) continue;
            if(distance_ponto_ponto(&particulas[i].posicao,&particulas[j].posicao)< particulas[i].raio + particulas[j].raio){

                done = done && force(&particulas[i],&particulas[j],&FORCE);
                
                if(FORCE.x != 0){
                    printf("Forca particula particula:");
                    print(&FORCE);
                }
                particulas[i].Force.x += FORCE.x;
                particulas[i].Force.y += FORCE.y;

                particulas[j].Force.x -= FORCE.x;
                particulas[j].Force.y -= FORCE.y;

            }
            
        }
    }
    return particulas;

}

struct particula* corrige_reta(struct particula* particulas,struct reta *retas,int N,struct VECTOR* CORRECOES){
    
    int i,j;
    bool done = true;
    struct VECTOR FORCE;
    struct VECTOR CM;
    CM.x = 150.0/2/1000;
    CM.y = 910.0/2/1000;
    FORCE.x = 0;
    FORCE.y = 0;
    for ( i = 0; i < N; i++){
        for ( j = 0; j < 6; j++){
            if(entre(&retas[j],&particulas[i].posicao)){
                
                
                if((distance_ponto_reta(&retas[j],&particulas[i].posicao) < particulas[i].raio)){
                    if((i == 1) &&(j == 2)) printf("Reta: %e %f\n",particulas[i].raio - distance_ponto_reta(&retas[j],&particulas[i].posicao),retas[j].b);
                    done = done && force_plano(&particulas[i],&retas[j],&FORCE);
                    if(FORCE.x != 0){
                        //printf("Forca Reta %d: %d\n",j,i);
                        //print(&FORCE);
                    }
                    particulas[i].Force.x += FORCE.x;
                    particulas[i].Force.y += FORCE.y;
                    FORCE.x = 0;
                    FORCE.y = 0;
                }
            }
        }
    }
    return particulas;
}

void simulate(int colunas,int linhas,double tempo_total,double angulo,double dt){
    int N = colunas*linhas,i,j;
    double t = 0.0;
    double PI = (4.0 * atan(1.0));
    struct particula* particulas = (struct particula*) malloc(N*sizeof(struct particula));
    struct VECTOR* anteriores = (struct VECTOR*) malloc(N*sizeof(struct VECTOR));
    struct VECTOR NORMAL;

    struct reta* retas = (struct reta*) malloc(6*sizeof(struct reta));
    double y0 = 98.0 + 154*tan(PI*angulo/180);
    //O4
    init_coef(&retas[0],0,0,0,98.0/1000);
    //O5
    init_coef(&retas[1],150./1000,0,150./1000,98.0/1000);
    //O0
    init_coef(&retas[2],0,98.0/1000,-154./1000,y0/1000);
    //O1
    init_coef(&retas[3],150./1000,98.0/1000,304./1000,y0/1000);
    //O2
    init_coef(&retas[4],-154./1000,y0/1000,-154./1000,910./1000);
    //O3
    init_coef(&retas[5],304./1000,y0/1000,304./1000,910./1000);

    for ( i = 0; i < 6; i++){
        retas[i].A = 0.01;
        retas[i].atrito = 0.;
        retas[i].gamma = 10;
    }
    for ( i = 0; i < N; i++){

        particulas[i].A = 0.01;
        particulas[i].aceleracao_angular = 0.;
        particulas[i].angular = 0.;
        particulas[i].atrito = 0.;
        particulas[i].gamma = 10.;
        particulas[i].raio = 7.5e-3;
        particulas[i].Young = 1e9;
        
        particulas[i].massa = 4./3.*PI*7860.*pow(particulas[i].raio,3.);
        particulas[i].Force.x = 0;
        particulas[i].Force.y = 0;

        particulas[i].velocidade.x = 0;
        particulas[i].velocidade.y = 0;

        particulas[i].aceleracao.x = 0;
        particulas[i].aceleracao.y = 0;
        
    }
    
    int c = 0;
    for ( i = 0; i < linhas; i++){
        for ( j = 0; j < colunas; j++){
            particulas[c].posicao.x = (-154+9+7 + (8+7.5)*j)/1000;
            particulas[c].posicao.y = (910/2 + (8+7.5)*i)/1000;
            anteriores[c].x = particulas[c].posicao.x;
            anteriores[c].y = particulas[c].posicao.y;
            c++;
        }
        
    }

    FILE *file = fopen("./example.txt","w");
    if (file == NULL) {
        printf("Erro ao abrir o arquivo.");
    }
    bool done = false;
    int count;
    struct VECTOR* CORRECOES = (struct VECTOR*)malloc(N*sizeof(struct VECTOR));
    while (t < tempo_total){

        
        count = 0;
        done = false;
        
        for ( i = 0; i < N; i++){
            if(particulas[i].posicao.y > -0.01){
                integracao(&particulas[i],&anteriores[i],dt);
                if((int)(t/dt)%20 == 0)printf("%f\t%d\t%f\t%f\t%f,%f\n",t,i,particulas[i].velocidade.x,particulas[i].velocidade.y,particulas[i].posicao.x,particulas[i].posicao.y);
                if((int)(t/dt)%500 == 0)fprintf(file,"%.4f %d %f %f\n",t,i,particulas[i].posicao.x,particulas[i].posicao.y);
            }
        }
        for ( i = 0; i < N; i++){
            CORRECOES[i].x = 0;
            CORRECOES[i].y = 0;
        }
        particulas = corrige_ponto(particulas,N,CORRECOES);

        particulas = corrige_reta(particulas,retas,N,CORRECOES);
        
        for ( i = 0; i < N; i++) if(particulas[i].posicao.y > -0.01) count++;
        //if(fabs(particulas[0].posicao.y) > -0.01) break;
        //fprintf(file,"%f\t%d\t%f\t%f\n",t,1,particulas[1].posicao.x,particulas[1].posicao.y);
        //printf("%f\n",t*1000);
        t += dt;
        if(count == 0) break;
    }
    free(particulas);
    free(anteriores);

}