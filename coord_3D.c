//
//  coord_3D.c
//  HexaCellLocal
//
//  Created by Felipe Viberti on 05/12/19.
//  Copyright © 2019 Felipe Viberti. All rights reserved.
//

#include "coord_3D.h"
#define H 0.001
double coord_fisica_3D(double s, double t, double r, double *v) {
    double soma = 0.0;
    double *N = (double *)malloc(8*sizeof(double));
    
    N[0] = (1.0/8.0)*(1-s)*(1-t)*(1-r);
    N[1] = (1.0/8.0)*(1-s)*(1+t)*(1-r);
    N[2] = (1.0/8.0)*(1+s)*(1+t)*(1-r);
    N[3] = (1.0/8.0)*(1+s)*(1-t)*(1-r);
    N[4] = (1.0/8.0)*(1-s)*(1-t)*(1+r);
    N[5] = (1.0/8.0)*(1-s)*(1+t)*(1+r);
    N[6] = (1.0/8.0)*(1+s)*(1+t)*(1+r);
    N[7] = (1.0/8.0)*(1+s)*(1-t)*(1+r);
    
    for (int i=0; i<8; i++)
        soma += N[i]*v[i];
    
    free(N);
    
    return soma;
}

double derivada_parcial_s_3(double p, double s, double t, double r, double *v) {
    double raiz_com_h = (p - coord_fisica_3D(s + H, t, r, v));
    double raiz_sem_h = (p - coord_fisica_3D(s, t, r, v));
    return (raiz_com_h - raiz_sem_h) / H;
}

double derivada_parcial_t_3(double p, double s, double t, double r, double *v) {
    double raiz_com_h = (p - coord_fisica_3D(s, t + H, r, v));
    double raiz_sem_h = (p - coord_fisica_3D(s, t, r, v));
    return (raiz_com_h - raiz_sem_h) / H;
}

double derivada_parcial_r_3(double p, double s, double t, double r, double *v) {
    double raiz_com_h = (p - coord_fisica_3D(s, t, r + H, v));
    double raiz_sem_h = (p - coord_fisica_3D(s, t, r, v));
    return (raiz_com_h - raiz_sem_h) / H;
}

int coord_parametrica_3D(double x, double y, double z, double *vx, double *vy, double *vz, double *s, double *t, double *r, double tol) {
    double *vf = (double *)malloc(3*sizeof(double));
    double *vh = (double *)malloc(3*sizeof(double));
    double *vs = (double *)malloc(3*sizeof(double));
    double **J = mat_cria(3, 3);
    
    int i = 0;
    
    if (vf == NULL || vh == NULL || vs == NULL) {
        puts("Erro malloc");
        exit(1);
    }
    
    vs[0] = vs[1] = vs[2] = 0;
    
    double coord_fisica_x = coord_fisica_3D(vs[0], vs[1], vs[2], vx);
    double coord_fisica_y = coord_fisica_3D(vs[0], vs[1], vs[2], vy);
    double coord_fisica_z = coord_fisica_3D(vs[0], vs[1], vs[2], vz);
    
    while (fabs(x - coord_fisica_x) >= tol || fabs(y - coord_fisica_y) >= tol || fabs(z - coord_fisica_z ) >= tol) {
        vf[0] = -(x - coord_fisica_3D(vs[0], vs[1], vs[2], vx));
        vf[1] =  -(y - coord_fisica_3D(vs[0], vs[1], vs[2], vy));
        vf[2] = -(z - coord_fisica_3D(vs[0], vs[1], vs[2], vz));
        
        J[0][0] = derivada_parcial_s_3(x,vs[0], vs[1], vs[2], vx);
        J[0][1] = derivada_parcial_t_3(x, vs[0], vs[1], vs[2], vx);
        J[0][2] = derivada_parcial_r_3(x, vs[0], vs[1], vs[2], vx);
        J[1][0] = derivada_parcial_s_3(y, vs[0], vs[1], vs[2], vy);
        J[1][1] = derivada_parcial_t_3(y, vs[0], vs[1], vs[2], vy);
        J[1][2] = derivada_parcial_r_3(y, vs[0], vs[1], vs[2], vy);
        J[2][0] = derivada_parcial_s_3(z, vs[0], vs[1], vs[2], vz);
        J[2][1] = derivada_parcial_t_3(z,vs[0], vs[1], vs[2], vz);
        J[2][2] = derivada_parcial_r_3(z,vs[0], vs[1], vs[2], vz);
        
        sist_linear(3, J, vf, vh);
        
        vs[0] += vh[0];
        vs[1] += vh[1];
        vs[2] += vh[2];
        
        coord_fisica_x = coord_fisica_3D(vs[0], vs[1], vs[2], vx);
        coord_fisica_y = coord_fisica_3D(vs[0], vs[1], vs[2], vy);
        coord_fisica_z = coord_fisica_3D(vs[0], vs[1], vs[2], vz);
        i++;
    }
    
    *s = vs[0];
    *t = vs[1];
    *r = vs[2];
    
    free(vf);
    free(vh);
    free(vs);
    mat_libera(3, J);
    
    return i;
}
