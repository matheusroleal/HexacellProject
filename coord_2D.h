//
//  coord_2D.h
//  HexaCellLocal
//
//  Created by Felipe Viberti on 05/12/19.
//  Copyright © 2019 Felipe Viberti. All rights reserved.
//

#ifndef coord_2D_h
#define coord_2D_h

#include <stdio.h>
#include <stdlib.h>
#include "sistlinear.h"
#include "matriz.h"

double coord_fisica_2D(double s, double t, double *v);

double derivada_parcial_s_2(double p, double s, double t, double *v);

double derivada_parcial_t_2(double p, double s, double t, double *v);

int coord_parametrica_2D(double x, double y, double *vx, double *vy, double *s, double *t, double tol);

#endif /* coord_2D_h */
