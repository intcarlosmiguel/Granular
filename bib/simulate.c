#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "particle.h"

const double GRAVIDADE = 9.80665;

void integracao(struct particula *p,struct VECTOR *anterior,double dt){
    double valor;
    p->Force.y = -p->massa*GRAVIDADE;
    mult(&p->Force,1/p->massa);
    p->velocidade.x += dt*(p->Force.x + p->aceleracao.x);
    p->velocidade.y += dt*(p->Force.y + p->aceleracao.y);
    valor = (p->velocidade.x != 0)? 2*p->posicao.x - anterior->x + pow(dt,2)*(p->Force.x+ p->aceleracao.x) : p->posicao.x + pow(dt,2)*(p->Force.x+ p->aceleracao.x);

    anterior->x = p->posicao.x;
    p->posicao.x = valor;

    valor= (p->velocidade.y != 0)? 2*p->posicao.y - anterior->y + (p->Force.y + p->aceleracao.y)*pow(dt,2) : p->posicao.y + pow(dt,2)*(p->Force.y + p->aceleracao.y);
    anterior->y = p->posicao.y;
    p->posicao.y = valor;

    copiar(&p->Force,&p->aceleracao);
    p->Force.x = 0;
    p->Force.y = 0;
}

bool corrige_ponto(struct particula* particulas,int N,struct VECTOR* CORRECOES){

    int i,j;
    bool done = true;
    struct VECTOR CORRECAO;
    CORRECAO.x = 0;
    CORRECAO.y = 0;
    for ( i = 0; i < N; i++){
        
        if(particulas[i].posicao.y < -10) continue;
        for ( j = i+1; j < N; j++){
            if(particulas[j].posicao.y < 0) continue;
            if(distance_ponto_ponto(&particulas[i].posicao,&particulas[j].posicao)< particulas[i].raio + particulas[j].raio){
                //printf("Ponto - Ponto : %f\n",particulas[i].raio + particulas[j].raio -distance_ponto_ponto(&particulas[i].posicao,&particulas[j].posicao) );
                done = done && force(&particulas[i],&particulas[j],&CORRECAO);
                sum(&CORRECOES[i],&CORRECAO,&CORRECOES[i]);
                mult(&CORRECAO,-1.0);
                sum(&CORRECOES[j],&CORRECAO,&CORRECOES[j]);
                CORRECAO.x = 0;
                CORRECAO.y = 0;
            }
        }
    }
    return done;

}

bool corrige_reta(struct particula* particulas,struct reta *retas,int N,struct VECTOR* CORRECOES){
    
    int i,j;
    bool done = true;
    struct VECTOR CORRECAO,NORMAL;
    struct VECTOR CM;
    CM.x = 150.0/2;
    CM.y = 910.0/2;
    CORRECAO.x = 0;
    CORRECAO.y = 0;
    for ( i = 0; i < N; i++){
        for ( j = 0; j < 6; j++){
            if(entre(&retas[j],&particulas[i].posicao)){
                
                if(!acima(&retas[j],&particulas[i].posicao,&CM)){
                    //printf("%d\n",j);
                    reflect_reta(&particulas[i].posicao,&retas[j],&particulas[i].posicao);
                }
                if((distance_ponto_reta(&retas[j],&particulas[i].posicao) < particulas[i].raio)){
                    done = done && force_plano(&particulas[i],&retas[j],&CORRECAO);
                    sum(&CORRECOES[i],&CORRECAO,&CORRECOES[i]);
                }
                //if(i == 0)print(&particulas[0].posicao);
            }
        }
    }
    return done;
}

void simulate(int colunas,int linhas,double tempo_total,double angulo){
    int N = colunas*linhas,i,j;
    double dt = 0.001,t = 0;
    double PI = (4.0 * atan(1.0));
    struct particula* particulas = (struct particula*) malloc(N*sizeof(struct particula));
    struct VECTOR* anteriores = (struct VECTOR*) malloc(N*sizeof(struct VECTOR));
    struct VECTOR NORMAL;

    struct reta* retas = (struct reta*) malloc(6*sizeof(struct reta));
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
        retas[i].atrito = 0.;
        retas[i].gamma = 0.05;
    }

    for ( i = 0; i < N; i++){

        particulas[i].A = 0.5;
        particulas[i].aceleracao_angular = 0.;
        particulas[i].angular = 0.;
        particulas[i].atrito = 0.;
        particulas[i].gamma = 0.5;
        particulas[i].massa = 0.1;
        particulas[i].raio = 7.5;
        particulas[i].Young = 1;
        
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
            particulas[c].posicao.x = -154+9+7 + (8+7.5)*j;
            particulas[c].posicao.y = 910/2 + (8+7.5)*i;
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
        
        for ( i = 0; i < N; i++)if(particulas[i].posicao.y > -20) integracao(&particulas[i],&anteriores[i],dt);
        while (!done){
            for ( i = 0; i < N; i++){
                CORRECOES[i].x = 0;
                CORRECOES[i].y = 0;
            }
            //print(&particulas[16].posicao);
            done = true;
            done = done && corrige_ponto(particulas,N,CORRECOES);

            done = done && corrige_reta(particulas,retas,N,CORRECOES);
            
            for ( i = 0; i < N; i++){
                sum(&particulas[i].posicao,&CORRECOES[i],&particulas[i].posicao);
            }
        }
        for ( i = 0; i < N; i++) if(particulas[i].posicao.y > -20) count++;
        for ( i = 0; i < N; i++) if(( (int)(t*1000)%100 == 0))fprintf(file,"%f\t%d\t%f\t%f\n",t,i,particulas[i].posicao.x,particulas[i].posicao.y);
        //fprintf(file,"%f\t%d\t%f\t%f\n",t,1,particulas[1].posicao.x,particulas[1].posicao.y);
        printf("%f\n",t*1000);
        t += dt;
        if(count == 0) break;
    }
    free(particulas);
    free(anteriores);

}