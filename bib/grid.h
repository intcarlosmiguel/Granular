#pragma once
#include <stdbool.h>
#include <stdlib.h>
#include "particle.h"
#include "define.h"

int X(struct GRID *grid,int celula){
    return (int) celula%grid->nx;
}
int Y(struct GRID *grid,int celula){
    return (int) celula/grid->nx;
}

int get_vizinho(struct GRID *grid,int celula,int check){
    switch (check){
        case 0:{ // Direita
            if(X(grid,celula) != grid->nx - 1) return celula+1;
            else return -1;
            break;
        }
        case 1: // Cima
            if(Y(grid,celula) != grid->ny - 1) return celula+ grid->nx;
            else return -1;
            break;
        case 2: // Esquerda
            if(X(grid,celula) != 0) return celula-1;
            else return -1;
            break;
        case 3: // Baixo
            if(Y(grid,celula) != 0) return celula - grid->nx;
            else return -1;
            break;
        case 4:
            if(X(grid,celula) != grid->nx - 1) get_vizinho(grid, celula+1,1);
            else return -1;
            break;
        case 5:
            if(X(grid,celula) != 0) return get_vizinho(grid, celula-1,1);
            else return -1;
            break;
        case 6:
            if(X(grid,celula) != 0) return get_vizinho(grid, celula-1,3);
            else return -1;
            break;
        case 7:
            if(X(grid,celula) != grid->nx - 1) get_vizinho(grid, celula+1,3);
            else return -1;
            break;
        
        default:
            return celula;
            break;
    }
    //return 0;
} 

void init_grid(struct GRID *grid,double x0,double y0,double x,double y,double largura){

    grid->largura = largura;
    double xf = x0,yf = y0;
    grid->x0 = x0;
    grid->y0 = y0;
    grid->ny = 0;
    grid->nx = 0;
    while (xf < x){
        xf += largura;
        grid->nx++;
    }

    while (yf < y){
        yf += largura;
        grid->ny++;
    }

    grid->n_grids = grid->nx*grid->ny;
    grid->celulas = (int**) malloc(grid->n_grids*sizeof(int*));
    
    for (int i = 0; i < grid->n_grids; i++) grid->celulas[i] = (int*) calloc(1,sizeof(int));

}

void atualiza_celula(struct GRID *grid,struct VECTOR *posicao,int celula_anterior,int particula){
    int x = (posicao->x - grid->x0)/grid->largura;
    int y = (posicao->y - grid->y0)/grid->largura;
    int nova_celula = x + y*grid->nx;
    if(nova_celula < 0){
        int last = grid->celulas[celula_anterior][0];
        int i;
        for (i = 1; i <= last; i++)if(grid->celulas[celula_anterior][i] == particula) break;
        int aux =  grid->celulas[celula_anterior][last];
        grid->celulas[celula_anterior][last] = particula;
        grid->celulas[celula_anterior][i] = aux;
        grid->celulas[celula_anterior] = realloc(grid->celulas[celula_anterior],(grid->celulas[celula_anterior][0])*sizeof(int));
        grid->celulas[celula_anterior][0] -= 1;
        grid->ids[particula] = nova_celula;
    }
    else{

        if((nova_celula != celula_anterior) && (nova_celula <grid->n_grids)){

            grid->ids[particula] = nova_celula;
            
            //printf("Finalizou %d\n",grid->celulas[nova_celula][0]);
            grid->celulas[nova_celula][0] += 1;
            grid->celulas[nova_celula] = realloc(grid->celulas[nova_celula],(grid->celulas[nova_celula][0]+1)*sizeof(int));
            grid->celulas[nova_celula][grid->celulas[nova_celula][0]] = particula;
            if(celula_anterior != -1){

                int last = grid->celulas[celula_anterior][0];
                int i;
                for (i = 1; i <= last; i++)if(grid->celulas[celula_anterior][i] == particula) break;
                int aux =  grid->celulas[celula_anterior][last];
                grid->celulas[celula_anterior][last] = particula;
                grid->celulas[celula_anterior][i] = aux;
                grid->celulas[celula_anterior] = realloc(grid->celulas[celula_anterior],(grid->celulas[celula_anterior][0])*sizeof(int));
                grid->celulas[celula_anterior][0] -= 1;
            }
        }
    }
}