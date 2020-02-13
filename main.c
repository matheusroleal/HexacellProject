//
//  main.c
//  HexaCellLocal
//
//  Created by Felipe Viberti on 05/12/19.
//  Copyright © 2019 Felipe Viberti. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matriz.h"
#include "sistlinear.h"
#include "coord_2D.h"
#include "coord_3D.h"
#define MAX 1
#define MIN 0
#define H 0.001  // h da derivada numérica

int main(void)
{
	int i;  // número de iterações
	int p = 18;  // número de casas decimais para precisão

	double vx_2[] = {MIN, MIN, MAX, MAX};  // vértices de x no caso 2D
	double vy_2[] = {MIN, MAX, MAX, MIN};  // vértices de y no caso 2D

	double vx_3[] = {MIN, MIN, MAX, MAX, MIN, MIN, MAX, MAX};  // vértices de x no caso 3D
	double vy_3[] = {MIN, MAX, MAX, MIN, MIN, MAX, MAX, MIN};  // vértices de y no caso 3D
	double vz_3[] = {MIN, MIN, MIN, MIN, MAX, MAX, MAX, MAX};  // vértices de z no caso 3D

	double s, t, r;
	double x, y, z;

	double tol = 0.5*pow(10,-p);  // tolerância com precisão de p casas decimais

	s = 1;
	t = 0;

	printf("\nCaso 2D:\n");
	printf("Usando tolerância de %d casas decimais:\n", p);
	printf("\nUsando valores de Xi = {%f, %f, %f, %f}\n", vx_2[0], vx_2[1], vx_2[2], vx_2[3]);
	printf("Usando valores de Yi = {%f, %f, %f, %f}\n", vy_2[0], vy_2[1], vy_2[2], vy_2[3]);
	printf("Usando s = %.12f e t = %.12f\n\n", s, t);

	x = coord_fisica_2D(s, t, vx_2);
	y = coord_fisica_2D(s, t, vy_2);

	i = coord_parametrica_2D(x, y, vx_2, vy_2, &s, &t, tol);

	printf("Resultado:\n");
	printf("x: %.12f\ny: %.12f\n", x, y);
	printf("s: %.22f\nt: %.22f\n", s, t);
	printf("iterações: %d\n", i);
    printf("------------\n\n");

	s = 1;
	t = 0;
	r = 0.5;

	printf("\nCaso 3D:\n");
	printf("Usando tolerância de %d casas decimais:\n", p);
	printf("\nUsando valores de Xi = {%f, %f, %f, %f, %f, %f, %f, %f}\n", vx_3[0], vx_3[1], vx_3[2], vx_3[3], vx_3[4], vx_3[5], vx_3[6], vx_3[7]);
	printf("Usando valores de Yi = {%f, %f, %f, %f, %f, %f, %f, %f}\n", vy_3[0], vy_3[1], vy_3[2], vy_3[3], vy_3[4], vy_3[5], vy_3[6], vy_3[7]);
	printf("Usando valores de Zi = {%f, %f, %f, %f, %f, %f, %f, %f}\n", vz_3[0], vz_3[1], vz_3[2], vz_3[3], vz_3[4], vz_3[5], vz_3[6], vz_3[7]);
	printf("Usando s = %.12f, t = %.12f e r = %.12f\n\n", s, t, r);

	x = coord_fisica_3D(s, t, r, vx_3);
	y = coord_fisica_3D(s, t, r, vy_3);
	z = coord_fisica_3D(s, t, r, vz_3);

	i = coord_parametrica_3D(x, y, z, vx_3, vy_3, vz_3, &s, &t, &r, tol);

	printf("Resultado:\n");
	printf("x: %.12f\ny: %.12f\nz: %.12f\n", x, y, z);
	printf("s: %.22f\nt: %.22f\nr: %.22f\n", s, t, r);
	printf("iterações: %d\n", i);

	return 0;
}
