#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "particle.h"
#include "grid.h"
double time;
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

struct particula* corrige_ponto(struct particula* particulas,int N,struct GRID *grid){

    int site,j;
    struct VECTOR FORCE;
    FORCE.x = 0;
    FORCE.y = 0;
    int id,viz,celula_vizinha,vizinho;
    bool a = false;
    for ( site = 0; site < N; site++){
        
        if(particulas[site].posicao.y < 0.01) continue;
        id = grid->ids[site];
        //printf("%d %d\n",site,id);
        if(id > 0){
            for ( viz = 0; viz < 9; viz++){

                celula_vizinha = get_vizinho(grid,id,viz);
                //if(site == 0)printf("Celula vizinha: %d\n",celula_vizinha);
                if((celula_vizinha > 0)){
                    //if(site == 0) if(celula_vizinha == 962) printf("%d\n",grid->celulas[celula_vizinha][0]);
                    if(grid->celulas[celula_vizinha][0] > 0){

                        for ( j = 1; j <= grid->celulas[celula_vizinha][0]; j++){

                            vizinho = grid->celulas[celula_vizinha][j];
                            if(site >= vizinho) continue;
                            if(particulas[vizinho].posicao.y < 0.01) continue;
                            if(distance_ponto_ponto(&particulas[site].posicao,&particulas[vizinho].posicao)< particulas[site].raio + particulas[vizinho].raio){
                                //if((site == 3)||(vizinho == 3)) printf("%d %d %f\n",site,vizinho,particulas[site].posicao.y);
                                force(&particulas[site],&particulas[vizinho],&FORCE);
                                particulas[site].Force.x += FORCE.x;
                                particulas[site].Force.y += FORCE.y;

                                particulas[vizinho].Force.x -= FORCE.x;
                                particulas[vizinho].Force.y -= FORCE.y;
                                FORCE.x = 0;
                                FORCE.y = 0;
                                a = true;
                            }
                        }
                        
                    }
                }
            }
            
        }
        
        /* for ( j = i+1; j < N; j++){
            if(particulas[j].posicao.y < 0.) continue;
            if(distance_ponto_ponto(&particulas[i].posicao,&particulas[j].posicao)< particulas[i].raio + particulas[j].raio){

                force(&particulas[i],&particulas[j],&FORCE);
                particulas[i].Force.x += FORCE.x;
                particulas[i].Force.y += FORCE.y;

                particulas[j].Force.x -= FORCE.x;
                particulas[j].Force.y -= FORCE.y;
                FORCE.x = 0;
                FORCE.y = 0;

            }
            
        } */
    }
    //if(a) exit(0);
    return particulas;

}

struct particula* corrige_reta(struct particula* particulas,struct reta *retas,int N){
    
    int i,j;
    struct VECTOR FORCE;
    struct VECTOR CM;
    CM.x = 150.0/2/1000;
    CM.y = 910.0/2/1000;
    FORCE.x = 0;
    FORCE.y = 0;
    for ( i = 0; i < N; i++){
        if(particulas[i].posicao.y <0) continue;
        for ( j = 0; j < 6; j++){
            if(entre(&retas[j],&particulas[i].posicao)){
                
                if(distance_ponto_reta(&retas[j],&particulas[i].posicao) < particulas[i].raio){
                    force_plano(&particulas[i],&retas[j],&FORCE);
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
    double L1 = 3.5*7.5e-3*2;
    double L2 = 154.e-3;
    double y0 = 98.0 + 154*tan(PI*angulo/180);

    //O4
    init_coef(&retas[0],0,0,0,98.0/1000);
    //O5
    init_coef(&retas[1],L1,0,L1,98.0/1000);
    //O0
    init_coef(&retas[2],0,98.0/1000,-154./1000,y0/1000);
    //O1
    init_coef(&retas[3],L1,98.0/1000,L1+L2,y0/1000);
    //O2
    init_coef(&retas[4],-154./1000,y0/1000,-154./1000,910./1000);
    //O3
    init_coef(&retas[5],L1+L2,y0/1000,L1+L2,910./1000);

    struct GRID grid;
    init_grid(&grid,-154./1000.,-2*7.5e-3,304./1000,910./1000,2*7.5e-3);
    grid.ids = (int*) calloc(N,sizeof(int));

    for (i = 0; i < N; i++)grid.ids[i] = -1;

    for ( i = 0; i < 6; i++){
        retas[i].A = 0.01;
        retas[i].atrito = 0.6;
        retas[i].gamma = 10;
    }
    for ( i = 0; i < N; i++){

        particulas[i].A = 0.01;
        particulas[i].aceleracao_angular = 0.;
        particulas[i].angular = 0.;
        particulas[i].atrito = 0.4;
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
            atualiza_celula(&grid,&particulas[c].posicao,grid.ids[c],c);
            anteriores[c].x = particulas[c].posicao.x;
            anteriores[c].y = particulas[c].posicao.y;
            c++;
        }
        
    }

    FILE *file = fopen("./example.txt","w");
    if (file == NULL) {
        printf("Erro ao abrir o arquivo.");
    }
    FILE *arquivo = fopen("./resultado.txt","w");
    if (file == NULL) {
        printf("Erro ao abrir o arquivo.");
    }
    bool done = false;
    int count = 0;
    bool* caiu = (bool*) calloc(N,sizeof(bool));
    double DT = 0;
    while (t < tempo_total){

        
        done = false;
        //printf("%d\n",grid.ids[1]);
        for ( i = 0; i < N; i++){
            if(particulas[i].posicao.y > 0.01){
                integracao(&particulas[i],&anteriores[i],dt);
                //if((int)(t/dt)%20 == 0)printf("%f\t%d\t%f\t%f\t%f,%f\n",t,i,particulas[i].velocidade.x,particulas[i].velocidade.y,particulas[i].posicao.x,particulas[i].posicao.y);
                if((int)(t/dt)%1000 == 0)fprintf(file,"%.4f %d %f %f\n",t,i,particulas[i].posicao.x,particulas[i].posicao.y);
                //printf("%d\n",i);
                if(grid.ids[i] > 0 )atualiza_celula(&grid,&particulas[i].posicao,grid.ids[i],i);
                
            }
            else{
                if(!caiu[i]){
                    //printf("Primeiro a sair:%d\n",i);
                    count++;
                    fprintf(arquivo,"%f %d\n",t - DT,count);
                    DT = t;
                    caiu[i] = true;
                }
            }
        }
        particulas = corrige_ponto(particulas,N,&grid);
        particulas = corrige_reta(particulas,retas,N);
        
        //if(fabs(particulas[0].posicao.y) > -0.01) break;
        //fprintf(file,"%f\t%d\t%f\t%f\n",t,1,particulas[1].posicao.x,particulas[1].posicao.y);
        //if((int)(t/dt)%1000 == 0)printf("%f\n",t);
        t += dt;
        time = t;
        if(count == N) break;
    }
    free(particulas);
    free(anteriores);
    fclose(arquivo);
    free(caiu);
    free(grid.ids);
}