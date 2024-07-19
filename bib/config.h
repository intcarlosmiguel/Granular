#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "define.h"

void init_values(int N,int n_retas,double atrito_retas,double atrito_particulas,struct reta* retas,struct particula* particulas){
    int i;
    for ( i = 0; i < n_retas; i++){
        retas[i].A = 0.01;
        retas[i].atrito = atrito_retas;
        retas[i].gamma = 10;
    }
    for ( i = 0; i < N; i++){

        particulas[i].A = 0.01;
        particulas[i].aceleracao_angular = 0.;
        particulas[i].angular = 0.;
        particulas[i].atrito = atrito_particulas;
        particulas[i].gamma = 10.;
        particulas[i].raio =0.0075;
        particulas[i].Young = 1e9;
        
        particulas[i].massa = 4./3.*PI*7860.*pow(particulas[i].raio,3.);
        particulas[i].Inertia = 2./5.*particulas[i].massa*pow(particulas[i].raio,2.);
        particulas[i].Force.x = 0;
        particulas[i].Force.y = 0;
        particulas[i].Force_Rot = 0;
        particulas[i].angular = 0;
        particulas[i].aceleracao_angular = 0;

        particulas[i].velocidade.x = 0;
        particulas[i].velocidade.y = 0;

        particulas[i].aceleracao.x = 0;
        particulas[i].aceleracao.y = 0;
        
    }
}