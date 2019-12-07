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
double coord_fisica_3D(double s, double t, double r, double *v);

double calcula_derivada_parcial_3D(double p, double s, double t, double r, double *v,int coordenada);

int coord_parametrica_3D(double x, double y, double z, double *vx, double *vy, double *vz, double *s, double *t, double *r, double tol);
#endif /* coord_3D_h */
