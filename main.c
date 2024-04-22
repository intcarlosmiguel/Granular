#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bib/vector.h"
#include "bib/particle.h"
const double GRAVIDADE = 9.80665;

void integracao(struct particula *p,struct VECTOR *anterior,double dt){
    double valor;
    p->velocidade.x += dt*p->Force.x/p->massa;
    p->velocidade.y += dt*p->Force.y/p->massa - GRAVIDADE*dt;
    valor =(p->velocidade.x != 0)? 2*p->posicao.x - anterior->x + pow(dt,2)*p->Force.x/p->massa : p->posicao.x + 1.5*pow(dt,2)*p->Force.x/p->massa;

    anterior->x = p->posicao.x;
    p->posicao.x = valor;

    valor= (p->velocidade.y != 0)? 2*p->posicao.y - anterior->y + (p->Force.y/p->massa - GRAVIDADE)*pow(dt,2) : p->posicao.y + 1.5*pow(dt,2)*p->Force.y/p->massa;

    anterior->y = p->posicao.y;
    p->posicao.y = valor;

    p->Force.x = 0;
    p->Force.y = 0;
}

bool corrige_ponto(struct particula* particulas,int N){
    int i,j;
    bool done = true;
    struct VECTOR CORRECAO;
    for ( i = 0; i < N; i++){
        if(particulas[i].posicao.y<-10) continue;
        for ( j = i+1; j < N; j++){
            if(particulas[j].posicao.y<0) continue;
            if(distance_ponto_ponto(&particulas[i].posicao,&particulas[j].posicao)<= particulas[i].raio + particulas[j].raio){

                done = done && force(&particulas[i],&particulas[j],&CORRECAO);

                adicionar(&CORRECAO,0.1);
                sum(&particulas[i].posicao,&CORRECAO,&particulas[i].posicao);
                mult(&CORRECAO,-1.0);
                sum(&particulas[j].posicao,&CORRECAO,&particulas[j].posicao);
                printf("Colisão da %d partícula com a partícula %d!\n",i,j);
            }
        }
    }
    return done;
}
bool corrige_reta(struct particula* particulas,struct reta *retas,int N){
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
            
            if((distance_ponto_reta(&retas[j],&particulas[i].posicao) < particulas[i].raio) && (acima(&retas[j],&particulas[i].posicao,&CM)) && (entre(&retas[j],&particulas[i].posicao))){
                NORMAL.x = retas[j].a;
                printf("Colisão da partícula %d com a parede %d %f %d!\n",i,j,distance_ponto_reta(&retas[j],&particulas[i].posicao),acima(&retas[j],&particulas[i].posicao,&CM));
                NORMAL.y = retas[j].b;
                mult(&NORMAL,1/norma(&NORMAL));

                done = done && force_plano(&particulas[i],&retas[j],&CORRECAO);

                /* double cos = dot(&NORMAL,&particulas[i].velocidade)/norma(&particulas[i].velocidade);

                particulas[i].velocidade.x = particulas[i].velocidade.x*sqrt(1 - pow(cos,2));
                particulas[i].velocidade.y = - particulas[i].velocidade.y*cos; */
            }
        }
    }
    for ( i = 0; i < N; i++){
        particulas[i].posicao.x = CORRECAO.x + particulas[i].posicao.x+0.1;
        particulas[i].posicao.y = CORRECAO.y + particulas[i].posicao.y+0.1;
    }

    return done;
}

void main(){
    double PI = (4.0 * atan(1.0));
    int colunas = 1;
    int linhas = 2;
    int N = colunas*linhas,i,j;
    double tempo_total = 1000, dt = 0.25/2,t = 0;

    struct particula* particulas = (struct particula*) malloc(N*sizeof(struct particula));
    struct VECTOR* anteriores = (struct VECTOR*) malloc(N*sizeof(struct VECTOR));
    struct VECTOR NORMAL;
    struct VECTOR CORRECAO;

    struct reta* retas = (struct reta*) malloc(6*sizeof(struct reta));
    double angulo = 30;
    double y0 = 9.8 + 154*tan(PI*angulo/180);
    printf("%f\n",y0);
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
        particulas[i].aceleracao_angular = 0;
        particulas[i].angular = 0;
        particulas[i].atrito = 0;
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
    while (t < tempo_total){
        count = 0;
        done = false;
        for ( i = 0; i < N; i++) integracao(&particulas[i],&anteriores[i],dt);
        while (!done){
            done = true;
            done = done && corrige_ponto(particulas,N);

            done = done && corrige_reta(particulas,retas,N);
        }
        for ( i = 0; i < N; i++) if(particulas[i].posicao.y > -20) count++;
        for ( i = 0; i < N; i++) fprintf(file,"%f\t%d\t%f\t%f\n",t,i,particulas[i].posicao.x,particulas[i].posicao.y);
        t += dt;
        if(count == 0) break;
        printf("%f\n",t);
    }
    free(particulas);
    free(anteriores);
}