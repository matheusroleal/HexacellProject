//
//  coord_3D.c
//  HexaCellLocal
//
//  Created by Felipe Viberti and Matheus Leal on 05/12/19.
//  Copyright © 2019 Felipe Viberti. All rights reserved.
//

#include "coord_3D.h"
#define H 0.001
double coord_fisica_3D(double s, double t, double r, double *v) {
    double soma = 0.0;
    double *N = (double *)malloc(8*sizeof(double));  // vetor com as funções de norma

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

double raiz_3D(double p, double s, double t, double r, double *v) {
    return (p - coord_fisica_3D(s, t, r, v));
}

double derivada_parcial_s_3(double p, double s, double t, double r, double *v) {
    double raiz_com_h = raiz_3D(p, s+H, t, r, v);
    double raiz_sem_h = raiz_3D(p,s,t,r,v);
    return (raiz_com_h - raiz_sem_h) / H;
}

double derivada_parcial_t_3(double p, double s, double t, double r, double *v) {
    double raiz_com_h = raiz_3D(p, s, t+H, r, v);
    double raiz_sem_h = raiz_3D(p,s,t,r,v);
    return (raiz_com_h - raiz_sem_h) / H;
}

double derivada_parcial_r_3(double p, double s, double t, double r, double *v) {
    double raiz_com_h = raiz_3D(p, s, t, r+H, v);
    double raiz_sem_h = raiz_3D(p,s,t,r,v);
    return (raiz_com_h - raiz_sem_h) / H;
}

int coord_parametrica_3D(double x, double y, double z, double *vx, double *vy, double *vz, double *s, double *t, double *r, double tol) {
    double *vf = (double *)malloc(3*sizeof(double));  // vetor com as funções u, v e w
    double *vh = (double *)malloc(3*sizeof(double));  // vetor passo a ser somado com a solução
    double *vs = (double *)malloc(3*sizeof(double));  // vetor solução
    double **J = matcria(3, 3);  // matriz Jacobiana

    int i = 0;

    if (vf == NULL || vh == NULL || vs == NULL) {
        puts("Erro de alocacao na memoria.");
        exit(1);
    }

    vs[0] = vs[1] = vs[2] = 0;

    while (fabs(x - coord_fisica_3D(vs[0], vs[1], vs[2], vx)) >= tol || fabs(y - coord_fisica_3D(vs[0], vs[1], vs[2], vy)) >= tol || fabs(z - coord_fisica_3D(vs[0], vs[1], vs[2], vz)) >= tol) {
        vf[0] = -raiz_3D(x, vs[0], vs[1], vs[2], vx);  // função u(s,t,r)
        vf[1] = -raiz_3D(y, vs[0], vs[1], vs[2], vy);  // função v(s,t,r)
        vf[2] = -raiz_3D(z, vs[0], vs[1], vs[2], vz);  // função w(s,t,r)

        J[0][0] = derivada_parcial_s_3(x,vs[0], vs[1], vs[2], vx);  // derivada parcial de u em s
        J[0][1] = derivada_parcial_t_3(x, vs[0], vs[1], vs[2], vx);  // derivada parcial de u em t
        J[0][2] = derivada_parcial_r_3(x, vs[0], vs[1], vs[2], vx);  // derivada parcial de u em r
        J[1][0] = derivada_parcial_s_3(y, vs[0], vs[1], vs[2], vy);  // derivada parcial de v em s
        J[1][1] = derivada_parcial_t_3(y, vs[0], vs[1], vs[2], vy);  // derivada parcial de v em t
        J[1][2] = derivada_parcial_r_3(y, vs[0], vs[1], vs[2], vy);  // derivada parcial de v em r
        J[2][0] = derivada_parcial_s_3(z, vs[0], vs[1], vs[2], vz);  // derivada parcial de w em s
        J[2][1] = derivada_parcial_t_3(z,vs[0], vs[1], vs[2], vz);  // derivada parcial de w em t
        J[2][2] = derivada_parcial_r_3(z,vs[0], vs[1], vs[2], vz);  // derivada parcial de w em r

        sist_linear(3, J, vf, vh);

        vs[0] += vh[0];
        vs[1] += vh[1];
        vs[2] += vh[2];

        i++;
    }

    *s = vs[0];
    *t = vs[1];
    *r = vs[2];

    free(vf);
    free(vh);
    free(vs);
    matlibera(3, J);

    return i;
}
