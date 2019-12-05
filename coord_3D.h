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

double raiz_3D(double p, double s, double t, double r, double *v);

double derivada_parcial_s_3(double p, double s, double t, double r, double *v);// derivada parcial em s para três variáveis

double derivada_parcial_t_3(double p, double s, double t, double r, double *v);

double derivada_parcial_r_3(double p, double s, double t, double r, double *v);


int coord_parametrica_3D(double x, double y, double z, double *vx, double *vy, double *vz, double *s, double *t, double *r, double tol);  // função que calcula (s,t,r) dados o par de pontos (x,y,z), os vetores dos vértices e a tolerância, retornando o número de iterações
#endif /* coord_3D_h */
