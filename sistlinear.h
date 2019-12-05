//
//  sistlinear.h
//  HexaCellLocal
//
//  Created by Felipe Viberti on 05/12/19.
//  Copyright © 2019 Felipe Viberti. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int* fatoracao (int n, double** a);
double* substituicao (int n, double** a, int* p, double* b);
void sist_linear(int n, double **A, double *b, double *x);
