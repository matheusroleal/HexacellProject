//
//  coord_3D.h
//  HexaCellLocal
//
//  Created by Felipe Viberti on 05/12/19.
//  Copyright © 2019 Felipe Viberti. All rights reserved.
//

#ifndef coord_3D_h
#define coord_3D_h

#include <stdio.h>
#include <stdlib.h>

#include "sistlinear.h"
#include "matriz.h"
double coord_fisica_3D(double s, double t, double r, double *v);  // função que calcula a coordenada física de um ponto (s,t,r) dado um vetor de vértices v

double raiz_3D(double p, double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v);  // função genérica definida como u(s,t,r), v(s,t,r) e w(s,t,r)

double derivada_parcial_s_3(double p, double (*raiz) (double p, double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v), double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v);  // derivada parcial em s para três variáveis

double derivada_parcial_t_3(double p, double (*raiz) (double p, double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v), double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v);  // derivada parcial em t para três variáveis

double derivada_parcial_r_3(double p, double (*raiz) (double p, double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v), double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v);  // derivada parcial em r para três variáveis


int coord_parametrica_3D(double x, double y, double z, double *vx, double *vy, double *vz, double *s, double *t, double *r, double tol);  // função que calcula (s,t,r) dados o par de pontos (x,y,z), os vetores dos vértices e a tolerância, retornando o número de iterações
#endif /* coord_3D_h */
