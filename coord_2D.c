//
//  coord_2D.c
//  HexaCellLocal
//
//  Created by Felipe Viberti on 05/12/19.
//  Copyright © 2019 Felipe Viberti. All rights reserved.
//

#include "coord_2D.h"
#define H 0.001
double coord_fisica_2D(double s, double t, double *v) {
    double soma = 0.0;
    double *N = (double *) malloc(4*sizeof(double));
    
    N[0] = (1.0/4.0)*(1-s)*(1-t);
    N[1] = (1.0/4.0)*(1-s)*(1+t);
    N[2] = (1.0/4.0)*(1+s)*(1+t);
    N[3] = (1.0/4.0)*(1+s)*(1-t);
    
    for (int i=0; i<4; i++)
        soma += N[i]*v[i];
    return soma;
}

double derivada_parcial_s_2(double p, double s, double t, double *v) {
    
    double raiz_com_h = raiz_2D(p, s+H, t, v);
    double raiz_sem_h = raiz_2D(p, s, t, v);
    return (raiz_com_h - raiz_sem_h) / H;
}

double derivada_parcial_t_2(double p, double s, double t, double *v) {
    double raiz_com_h = raiz_2D(p, s, t+H, v);
    double raiz_sem_h = raiz_2D(p, s, t, v);
    return (raiz_com_h - raiz_sem_h) / H;
}

double raiz_2D(double p,double s,double t, double *v)
{
    return (p - coord_fisica_2D(s, t, v));
}

int coord_parametrica_2D(double x, double y, double *vx, double *vy, double *s, double *t, double tol)
{
    double *vf = (double *)malloc(2*sizeof(double));  // vetor com as funções u e v definidas no enunciado
    double *vh = (double *)malloc(2*sizeof(double));  // vetor passo a ser somado com a solução
    double *vs = (double *)malloc(2*sizeof(double));  // vetor solução
    double **J = mat_cria(2, 2);  // matriz Jacobiana
    
    int i = 0;
    
    if (vf == NULL || vh == NULL || vs == NULL) {
        puts("Erro de alocacao na memoria.");
        exit(1);
    }
    
    vs[0] = vs[1] = 0;
    
    while (fabs(x - coord_fisica_2D(vs[0], vs[1], vx)) >= tol || fabs(y - coord_fisica_2D(vs[0], vs[1], vy)) >= tol) {
        vf[0] = -raiz_2D(x, vs[0], vs[1], vx);  // função u(s,t)
        vf[1] = -raiz_2D(y, vs[0], vs[1], vy);  // função v(s,t)
        
        J[0][0] = derivada_parcial_s_2(x, vs[0], vs[1], vx);  // derivada parcial de u em s
        J[0][1] = derivada_parcial_t_2(x, vs[0], vs[1], vx);  // derivada parcial de u em t
        J[1][0] = derivada_parcial_s_2(y, vs[0], vs[1], vy);  // derivada parcial de v em s
        J[1][1] = derivada_parcial_t_2(y, vs[0], vs[1], vy);  // derivada parcial de v em t
        
        sist_linear(2, J, vf, vh);
        
        vs[0] += vh[0];
        vs[1] += vh[1];
        
        i++;
    }
    
    *s = vs[0];
    *t = vs[1];
    
    free(vf);
    free(vh);
    free(vs);
    mat_libera(2, J);
    
    return i;
}
