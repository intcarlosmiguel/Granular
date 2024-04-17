#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "particle.h"
//#include "obstaculo.h"
double PI = 3.14159265359;
const double GRAVIDADE = 9.80665;

void update_position(struct particula *p,struct VECTOR *anterior,struct VECTOR *F,double dt){
    //p->velocidade.x += dt*p->Force.x;
    //p->velocidade.y += dt*p->Force.y - GRAVIDADE;
    anterior->x = (2 - p->Force.x)*p->posicao.x - (1 - p->Force.x)*anterior->x + dt*dt*p->Force.x/p->massa;
    anterior->y = (2 - p->Force.y - p->massa*GRAVIDADE)*p->posicao.y - (1 - p->Force.y- p->massa*GRAVIDADE)*anterior->y + (p->Force.y/p->massa - GRAVIDADE)*dt*dt;
    F->x = p->Force.x;
    F->y = p->Force.y;
    p->Force.x = 0;
    p->Force.y = 0;
}

void main(){
    int colunas = 20;
    int linhas = 20;
    int N = colunas*linhas,i,j;
    double tempo_total = 200, dt = 0.25,t = 0;
    struct particula* particulas = (struct particula*) malloc(N*sizeof(struct particula));
    struct VECTOR* anteriores = (struct VECTOR*) malloc(N*sizeof(struct VECTOR));
    struct VECTOR* F = (struct VECTOR*) malloc(N*sizeof(struct VECTOR));
    struct particula PAREDE;
    struct reta* retas = (struct reta*) malloc(6*sizeof(struct reta));
    double angulo = 30;
    double y0 = 9.8 + 154*tan(PI*angulo/180);

    //O4
    init_coef(&retas[0],0,0,0,9.8);
    //O5
    init_coef(&retas[1],150,0,150,9.8);
    //O0
    init_coef(&retas[2],0,9.8,-154,y0);
    //O1
    init_coef(&retas[3],150,9.8,304,y0);
    //O2
    init_coef(&retas[4],-154,y0,-154,910);
    //O3
    init_coef(&retas[5],304,y0,304,910);

    for ( i = 0; i < 6; i++){
        retas[i].A = 0.01;
        retas[i].atrito = 0.05;
        retas[i].gamma = 0.05;
    }

    for ( i = 0; i < N; i++){

        particulas[i].A = 0;
        particulas[i].aceleracao_angular = 0;
        particulas[i].angular = 0;
        particulas[i].atrito = 1;
        particulas[i].gamma = 1;
        particulas[i].massa = 1;
        particulas[i].raio = 7.5;
        particulas[i].Young = 1;

        particulas[i].Force.x = 0;
        particulas[i].Force.y = 0;
        particulas[i].velocidade.x = 0;
        particulas[i].velocidade.y = 0;
        
    }
    int c = 0;
    for ( i = 0; i < linhas; i++){
        for ( j = 0; j < colunas; j++){
            //printf("%d\n",-154+9 + 8*j);
            particulas[c].posicao.x = -154+9 + (8+7.5)*j;
            particulas[c].posicao.y = 910/2 + (8+7.5)*i;
            //printf("%f - %f\n",particulas[c].posicao.x,particulas[c].posicao.y);
            anteriores[c].x = particulas[c].posicao.x;
            anteriores[c].y = particulas[c].posicao.y;
            c++;
        }
        
    }
    for ( i = 0; i < N; i++){
        //printf("%d %f - %f\n",i,particulas[i].posicao.x,particulas[i].posicao.y);
        
    }

    while (t < tempo_total){
        for ( i = 0; i < N; i++){
            //if(i == 0)printf("Particula: %d - (%f,%f)\n",i,particulas[i].posicao.x,particulas[i].posicao.y);
            if(particulas[i].posicao.y<0) continue;
            for ( j = i+1; j < N; j++){
                if(particulas[j].posicao.y<0) continue;
                if(distance_ponto_ponto(&particulas[i].posicao,&particulas[j].posicao)<= particulas[i].raio + particulas[j].raio) {
                    if(i == 0)printf("Entrou %d\n",j);
                    force(&particulas[i],&particulas[j]);
                }
            }
            //if(i == 0)printf("Particula: %d - (%f,%f)\n",i,particulas[i].Force.x,particulas[i].Force.y);
            update_position(&particulas[i],&anteriores[i],&F[i],dt);
            //if(i == 0)printf("Particula: %d - (%f,%f)\n",i,anteriores[i].x,anteriores[i].y);
            for ( j = 0; j < 6; j++){
                if(entre(&retas[j],&particulas[i].posicao)){
                    if(!acima(&retas[j],&particulas[i].posicao)) reflect_reta(&particulas[i].posicao,&retas[j],&particulas[i].posicao);

                    if(distance_ponto_reta(&retas[j],&particulas[i].posicao) < particulas[i].raio)force_plano(&particulas[i],&retas[j]); 
                } 
            }
            
        }
        for ( i = 0; i < N; i++){
            particulas[i].posicao.x = anteriores[i].x;
            particulas[i].posicao.y = anteriores[i].y;
            particulas[i].velocidade.x = dt*F[i].x;
            particulas[i].velocidade.y = dt*F[i].y - GRAVIDADE;
        }

        t += dt;
        printf("%f\n",t);
    }

}