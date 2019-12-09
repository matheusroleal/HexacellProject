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

double calcula_derivada_parcial_3D(double p, double s, double t, double r, double *v,int coordenada) {
    double raiz_com_h,raiz_sem_h;
    if(coordenada == 0) {
        raiz_com_h = (p - coord_fisica_3D(s + H, t, r, v));
    }
    else if (coordenada == 1) {
        raiz_com_h = (p - coord_fisica_3D(s, t + H, r, v));
    }
    else if (coordenada == 2){
        raiz_com_h = (p - coord_fisica_3D(s, t, r + H, v));
    }
    raiz_sem_h = (p - coord_fisica_3D(s, t, r, v));
    return (raiz_com_h - raiz_sem_h) / H;
}



int coord_parametrica_3D(double x, double y, double z, double *vx, double *vy, double *vz, double *s, double *t, double *r, double tol) {
    double *vf = (double *)malloc(3*sizeof(double));
    double *vh = (double *)malloc(3*sizeof(double));
    double *vs = (double *)malloc(3*sizeof(double));
    double **J = matcria(3, 3);

    int i = 0;

    if (vf == NULL || vh == NULL || vs == NULL) {
        puts("Erro malloc");
        exit(1);
    }

    vs[0] = 0;
    vs[1] = 0;
    vs[2] = 0;

    double coord_fisica_x = coord_fisica_3D(vs[0], vs[1], vs[2], vx);
    double coord_fisica_y = coord_fisica_3D(vs[0], vs[1], vs[2], vy);
    double coord_fisica_z = coord_fisica_3D(vs[0], vs[1], vs[2], vz);

    while (fabs(x - coord_fisica_x) >= tol || fabs(y - coord_fisica_y) >= tol || fabs(z - coord_fisica_z ) >= tol) {
        vf[0] = -(x - coord_fisica_3D(vs[0], vs[1], vs[2], vx));
        vf[1] =  -(y - coord_fisica_3D(vs[0], vs[1], vs[2], vy));
        vf[2] = -(z - coord_fisica_3D(vs[0], vs[1], vs[2], vz));

        J[0][0] = calcula_derivada_parcial_3D(x,vs[0], vs[1], vs[2], vx,0);
        J[0][1] = calcula_derivada_parcial_3D(x, vs[0], vs[1], vs[2], vx,1);
        J[0][2] = calcula_derivada_parcial_3D(x, vs[0], vs[1], vs[2], vx,2);
        J[1][0] = calcula_derivada_parcial_3D(y, vs[0], vs[1], vs[2], vy,0);
        J[1][1] = calcula_derivada_parcial_3D(y, vs[0], vs[1], vs[2], vy,1);
        J[1][2] = calcula_derivada_parcial_3D(y, vs[0], vs[1], vs[2], vy,2);
        J[2][0] = calcula_derivada_parcial_3D(z, vs[0], vs[1], vs[2], vz,0);
        J[2][1] = calcula_derivada_parcial_3D(z,vs[0], vs[1], vs[2], vz,1);
        J[2][2] = calcula_derivada_parcial_3D(z,vs[0], vs[1], vs[2], vz,2);

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
    matlibera(3, J);

    return i;
}
