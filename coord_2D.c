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
    
    double raiz_com_h = (p - coord_fisica_2D(s + H, t, v));
    double raiz_sem_h = (p - coord_fisica_2D(s, t, v));
    return (raiz_com_h - raiz_sem_h) / H;
}

double derivada_parcial_t_2(double p, double s, double t, double *v) {
    double raiz_com_h = (p - coord_fisica_2D(s, t + H, v));
    double raiz_sem_h = (p - coord_fisica_2D(s, t, v));
    return (raiz_com_h - raiz_sem_h) / H;
}

int coord_parametrica_2D(double x, double y, double *vx, double *vy, double *s, double *t, double tol) {
    double *vf = (double *)malloc(2*sizeof(double));
    double *vh = (double *)malloc(2*sizeof(double));
    double *vs = (double *)malloc(2*sizeof(double));
    double **J = mat_cria(2, 2);
    
    int i = 0;
    
    if (vf == NULL || vh == NULL || vs == NULL) {
        puts("Erro malloc");
        exit(1);
    }
    
    vs[0] = 0;
    vs[1] = 0;
    
    double coord_fis_x = coord_fisica_2D(vs[0], vs[1], vx);
    double coord_fis_y = coord_fisica_2D(vs[0], vs[1], vy);
    
    while (fabs(x - coord_fis_x) >= tol || fabs(y - coord_fis_y) >= tol) {
        vf[0] = -(x - coord_fisica_2D(vs[0], vs[1], vx));
        vf[1] = -(y - coord_fisica_2D(vs[0], vs[1], vy));
        
        J[0][0] = derivada_parcial_s_2(x, vs[0], vs[1], vx);
        J[0][1] = derivada_parcial_t_2(x, vs[0], vs[1], vx);
        J[1][0] = derivada_parcial_s_2(y, vs[0], vs[1], vy);
        J[1][1] = derivada_parcial_t_2(y, vs[0], vs[1], vy);
        
        sist_linear(2, J, vf, vh);
        
        vs[0] += vh[0];
        vs[1] += vh[1];
        
        coord_fis_x = coord_fisica_2D(vs[0], vs[1], vx);
        coord_fis_y = coord_fisica_2D(vs[0], vs[1], vy);
        
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
