//
//  matriz.c
//  HexaCellLocal
//
//  Created by Felipe Viberti on 05/12/19.
//  Copyright © 2019 Felipe Viberti. All rights reserved.
//

#include "matriz.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double** mat_cria (int m, int n) {
  double** A = (double**) malloc(m*sizeof(double*));
  for (int i=0; i<m; ++i)
    A[i] = (double*) malloc(n*sizeof(double));
  return A;
}

void mat_libera (int m, double** A) {
  for (int i=0; i<m; ++i)
    free(A[i]);
  free(A);
}

void mat_transposta (int m, int n, double** A, double** T) {
  for (int i=0; i<m; ++i)
    for (int j=0; j<n; ++j)
      T[j][i] = A[i][j];
}
  
void mat_multv (int m, int n, double** A, double* v, double* w) {
  for (int i=0; i<m; ++i) {
    w[i] = 0.0;
    for (int j=0; j<n; ++j)
      w[i] += A[i][j] * v[j];
  }
}
  
void mat_multm (int m, int n, int q, double** A, double** B, double** C) {
  for (int i=0; i<m; ++i) {
    for (int k=0; k<q; ++k) {
      C[i][k] = 0.0;
      for (int j=0; j<n; ++j) 
        C[i][k] += A[i][j] * B[j][k];
    }
  }
}
  
int mat_iguais (int m, int n, double** A, double** B, double tol) {
  for (int i=0; i<m; ++i)
    for (int j=0; j<n; ++j)
      if (fabs(A[i][j] - B[i][j]) > tol)
        return 0;
  return 1;
}
  
void mat_imprime (int m, int n, double** A, char* format) {
  for (int i=0; i<m; ++i) {
    for (int j=0; j<n; ++j) {
      printf(format,A[i][j]);
      if (j==n-1)
        printf("\n");
      else
        printf(" ");
    }
  }
}
  
