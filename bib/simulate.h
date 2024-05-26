#pragma once
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <string.h>
#include <errno.h>
#include "config.h"
#include "particle.h"
#include "grid.h"
#include "define.h"
#include "mtwister.h"
#include "integrate.h"
void create_directory_if_not_exists(const char *folder_name) {
    struct stat st = {0};

    // Verificar se o diretório já existe
    if (stat(folder_name, &st) == -1) {
        // Diretório não existe, então tenta criar
        if (mkdir(folder_name, 0700) == -1) {
            perror("Erro ao criar a pasta");
        } else {
            printf("Pasta '%s' criada com sucesso!\n", folder_name);
        }
    } else {
        // Diretório já existe
        //printf("Pasta '%s' já existe, não foi criada novamente.\n", folder_name);
    }
}

struct particula* corrige_ponto(struct particula* particulas,int N,struct GRID *grid,bool rotacao){

    int site,j;
    struct VECTOR FORCE;
    FORCE.x = 0;
    FORCE.y = 0;
    int id,viz,celula_vizinha,vizinho;
    double force_rotacao = 0;
    for ( site = 0; site < N; site++){
        if(particulas[site].posicao.y <0.01) continue;
        particulas = calc_force_par(site,particulas,N,grid,rotacao);
    }
    return particulas;

}

struct particula* corrige_reta(struct particula* particulas,struct reta *retas,int N,bool rotacao){
    
    int i,j;
    struct VECTOR FORCE;
    FORCE.x = 0;
    FORCE.y = 0;
    double force_rotacao = 0;
    for ( i = 0; i < N; i++){
        if(particulas[i].posicao.y <0) continue;
        particulas = calc_force_reta(particulas,i,retas,N,rotacao);
    }
    return particulas;
}

void simulate(int colunas,int linhas,double tempo_total,double angulo,double dt, double atrito_particulas, double atrito_retas,int seed,bool rotacao){
    PI = (4.0 * atan(1.0));
    int N = colunas*linhas,i,j;
    double t = 0.0;
    struct particula* particulas = (struct particula*) malloc(N*sizeof(struct particula));
    struct VECTOR* anteriores = (struct VECTOR*) malloc(N*sizeof(struct VECTOR));

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

    init_values(N,6,atrito_retas,atrito_particulas,retas,particulas);

    init_genrand64(seed);
    
    int c = 0;
    for ( i = 0; i < linhas; i++){
        for ( j = 0; j < colunas; j++){
            particulas[c].posicao.x = (-154+9+7 + (2*7.5+0.2+genrand64_real1()/10)*j)/1000;
            particulas[c].posicao.y = (910/3 + (2*7.5+0.2+genrand64_real1()/10)*i)/1000;
            atualiza_celula(&grid,&particulas[c].posicao,grid.ids[c],c);
            anteriores[c].x = particulas[c].posicao.x;
            anteriores[c].y = particulas[c].posicao.y;
            c++;
        }
        
    }

    char folder_name[500]; 
    sprintf(folder_name, "./results/%d", (int) angulo);
    create_directory_if_not_exists(folder_name);

    char string[200];
    char example[200];
    sprintf(example, "./results/%d/example_%.2f_%.2f.dat",(int) angulo, atrito_particulas,atrito_retas);
    FILE *file = fopen(example,"r");
    bool arquivo_criado = false;
    if (file != NULL) {
        // Se o arquivo existe, fechamos e não fazemos nada
        fclose(file);
    } else {
        // Se o arquivo não existe, cria ele
        file = fopen(example, "w");
        arquivo_criado = true;
        // Fecha o arquivo
    }


    int count = 0;
    bool* caiu = (bool*) calloc(N,sizeof(bool));
    double** resultados = (double**) malloc(N*sizeof(double*));
    for ( i = 0; i < N; i++) resultados[i] = calloc(2,sizeof(double));
    double DT = 0;
    while (t < tempo_total){
        for ( i = 0; i < N; i++){
            if(particulas[i].posicao.y > 0.01){
                integracao(&particulas[i],&anteriores[i],dt);
                if(arquivo_criado)if((int)(t/dt)%5000 == 0)fprintf(file,"%.4f %d %f %f\n",t,i,particulas[i].posicao.x,particulas[i].posicao.y);
                if(grid.ids[i] > 0 )atualiza_celula(&grid,&particulas[i].posicao,grid.ids[i],i);
                
            }
            else{
                if(!caiu[i]){
                    count++;
                    resultados[count - 1][0] = t - DT;
                    resultados[count - 1][1] = count;
                    DT = t;
                    caiu[i] = true;
                }
            }
        }
        particulas = corrige_ponto(particulas,N,&grid,rotacao);
        particulas = corrige_reta(particulas,retas,N,rotacao);
        t += dt;
        if(count == N) break;
        if(t - DT > 2) break;
    }
    
    if(count == N) sprintf(string, "./results/%d/resultado_%.2f_%.2f.dat",(int) angulo, atrito_particulas,atrito_retas);
    else{
        printf("%d\n",seed);
        sprintf(string, "./results/%d/resultado_%.2f_%.2f_stop.dat",(int) angulo, atrito_particulas,atrito_retas);
    }
    FILE *arquivo = fopen(string,"a");
    if (arquivo == NULL) {
        printf("Erro ao abrir o arquivo.");
    }

    for ( i = 0; i < count; i++){
        fprintf(arquivo,"%f %d\n",resultados[i][0],(int)resultados[i][1]);
        free(resultados[i]);
    }
    free(resultados);
    free(particulas);
    free(anteriores);
    fclose(arquivo);
    free(caiu);
    free(grid.ids);
    for ( i = 0; i < grid.n_grids; i++){
        free(grid.celulas[i]);
    }
    free(grid.celulas);
}