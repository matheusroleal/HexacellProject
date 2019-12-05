#include "sistlinear.h"
#include <stdlib.h>
#include <math.h>
 
int* fatoracao (int n, double** a)
{
  int* p = (int*) malloc(n*sizeof(int));
  for (int i=0; i<n; ++i)
    p[i] = i;
  for (int j=0; j<n-1; ++j) {
    // find pivô
    int m = j;
    for (int i=j+1; i<n; ++i)
      if (fabs(a[i][j]) > fabs(a[m][j]))
        m = i;
    // swap lines: j <-> m
    for (int k=0; k<n; ++k) {
      double t = a[j][k];
      a[j][k] = a[m][k];
      a[m][k] = t;
    }
    // register permutation
    int t = p[j];
    p[j] = p[m];
    p[m] = t;
    // elimination
    for (int i=j+1; i<n; ++i) {
      double f = a[i][j]/a[j][j];
      for (int k=j+1; k<n; ++k)
        a[i][k] -= f*a[j][k];
      a[i][j] = f;
    }
  }
  return p;
}

double* substituicao (int n, double** a, int* p, double* b)
{
  double* x = (double*) malloc(n*sizeof(double));
  // forward substitution
  for (int i=0; i<n; ++i) {
    double s = 0;
    for (int j=0; j<i; ++j) 
      s += a[i][j]*x[j];
    x[i] = b[p[i]] - s;
  }
  // backward substitution
  for (int i=n-1; i>=0; --i) {
    double s = 0;
    for (int j=i+1; j<n; ++j) 
      s += a[i][j]*x[j];
    x[i] = (x[i] - s) / a[i][i];
  }
  return x;
}


