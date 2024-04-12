#include <stdio.h>
#include <stdlib.h>
#include "vector.h"
#include "particle.h"
//#include "obstaculo.h"

const double GRAVIDADE = 9.80665;

void update_position(struct particula *p,struct VECTOR *anterior,double dt){
    p->velocidade.x += dt*p->Force.x;
    p->velocidade.y += dt*p->Force.y - GRAVIDADE;
    p->posicao.x = (2 - p->Force.x)*p->posicao.x - (1 - p->Force.x)*anterior->x + dt*dt*p->Force.x/p->massa;
    p->posicao.y = (2 - p->Force.y)*p->posicao.y - (1 - p->Force.y)*anterior->y + p->Force.y/p->massa;
    p->Force.x = 0;
    p->Force.y = 0;
}

void main(){

    int N = 25,i,j;
    double tempo_total = 200, dt = 0.25,t = 0;
    struct particula* particulas = (struct particula*) malloc(N*sizeof(struct particula));
    struct VECTOR* anteriores = (struct VECTOR*) malloc(N*sizeof(struct VECTOR));
    struct particula PAREDE;
    struct reta reta1;
    init_coef(&reta1,0,0,1,1);
    for ( i = 0; i < N; i++){

        particulas[i].A = 0;
        particulas[i].aceleracao_angular = 0;
        particulas[i].angular = 0;
        particulas[i].atrito = 0.5;
        particulas[i].gamma = 0.5;
        particulas[i].massa = 0.5;
        particulas[i].raio = 0.5;
        particulas[i].Young = 0.5;

        particulas[i].Force.x = 0;
        particulas[i].Force.y = 0;
        particulas[i].posicao.x = 0;
        particulas[i].posicao.y = 0;
        particulas[i].velocidade.x = 0;
        particulas[i].velocidade.y = 0;
        anteriores[i].x = particulas[i].posicao.x;
        anteriores[i].y = particulas[i].posicao.y;

        reta1.A = 0.01;
        reta1.atrito = 0.5;
        reta1.gamma = 10;
    }
    
    while (t < tempo_total){
        for ( i = 0; i < N; i++){

            for ( j = i+1; j < N; j++) if(distance_ponto_ponto(&particulas[i].posicao,&particulas[j].posicao)<= particulas[i].raio + particulas[j].raio) force(&particulas[i],&particulas[j]);
            
            update_position(&particulas[i],&anteriores[i],dt);

            if(entre(&reta1,&particulas[i].posicao)){
                if(!acima(&reta1,&particulas[i].posicao)) reflect_reta(&particulas[i].posicao,&reta1,&particulas[i].posicao);

                if(distance_ponto_reta(&reta1,&particulas[i].posicao) < particulas[i].raio)force_plano(&particulas[i],&reta1); 
            }
            anteriores[i].x = particulas[i].posicao.x;
            anteriores[i].y = particulas[i].posicao.y;
        }
        t += dt;
    }
    

}