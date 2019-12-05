#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matriz.h"
#include "sistlinear.h"

#define MAX 6
#define MIN 0
#define H 0.001  // h da derivada numérica

double coordenada_fisica_2D(double s, double t, double *v);  // função que calcula a coordenada física de um ponto (s,t) dado um vetor de vértices v

double coordenada_fisica_3D(double s, double t, double r, double *v);  // função que calcula a coordenada física de um ponto (s,t,r) dado um vetor de vértices v

double raiz_2D(double p, double (*coord_fis_2D) (double s, double t, double *v), double s, double t, double *v);  // função genérica definida como u(s,t) e v(s,t) expressas no enunciado

double raiz_3D(double p, double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v);  // função genérica definida como u(s,t,r), v(s,t,r) e w(s,t,r)

double dps_2(double p, double (*raiz) (double p, double (*coord_fis_2D) (double s, double t, double *v), double s, double t, double *v), double (*coord_fis_2D) (double s, double t, double *v), double s, double t, double *v);  // derivada parcial em s para duas variáveis

double dpt_2(double p, double (*raiz) (double p, double (*coord_fis_2D) (double s, double t, double *v), double s, double t, double *v), double (*coord_fis_2D) (double s, double t, double *v), double s, double t, double *v);  // derivada parcial em t para duas variáveis

double dps_3(double p, double (*raiz) (double p, double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v), double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v);  // derivada parcial em s para três variáveis

double dpt_3(double p, double (*raiz) (double p, double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v), double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v);  // derivada parcial em t para três variáveis

double dpr_3(double p, double (*raiz) (double p, double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v), double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v);  // derivada parcial em r para três variáveis

void sist_linear(int n, double **A, double *b, double *x);  // função que resolve um sistema linear armazenando a solução no vetor x

int coordenada_parametrica_2D(double x, double y, double *vx, double *vy, double *s, double *t, double tol);  // função que calcula (s,t) dados o par de pontos (x,y), os vetores dos vértices e a tolerância, retornando o número de iterações

int coordenada_parametrica_3D(double x, double y, double z, double *vx, double *vy, double *vz, double *s, double *t, double *r, double tol);  // função que calcula (s,t,r) dados o par de pontos (x,y,z), os vetores dos vértices e a tolerância, retornando o número de iterações

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
	t = -0.5;
	
	printf("\nCaso 2D:\n");
	printf("Usando tolerância de %d casas decimais:\n", p);
	printf("\nUsando valores de Xi = {%f, %f, %f, %f}\n", vx_2[0], vx_2[1], vx_2[2], vx_2[3]);
	printf("Usando valores de Yi = {%f, %f, %f, %f}\n", vy_2[0], vy_2[1], vy_2[2], vy_2[3]);
	printf("Usando s = %.12f e t = %.12f\n\n", s, t);
	
	x = coordenada_fisica_2D(s, t, vx_2);
	y = coordenada_fisica_2D(s, t, vy_2);

	i = coordenada_parametrica_2D(x, y, vx_2, vy_2, &s, &t, tol);
	
	printf("Resultado:\n");
	printf("x: %.12f\ny: %.12f\n", x, y);
	printf("s: %.22f\nt: %.22f\n", s, t);
	printf("iterações: %d\n", i);
    printf("------------\n\n");

	s = 1;
	t = -0.5;
	r = 0.5;
	
	printf("\nCaso 3D:\n");
	printf("Usando tolerância de %d casas decimais:\n", p);
	printf("\nUsando valores de Xi = {%f, %f, %f, %f, %f, %f, %f, %f}\n", vx_3[0], vx_3[1], vx_3[2], vx_3[3], vx_3[4], vx_3[5], vx_3[6], vx_3[7]);
	printf("Usando valores de Yi = {%f, %f, %f, %f, %f, %f, %f, %f}\n", vy_3[0], vy_3[1], vy_3[2], vy_3[3], vy_3[4], vy_3[5], vy_3[6], vy_3[7]);
	printf("Usando valores de Zi = {%f, %f, %f, %f, %f, %f, %f, %f}\n", vz_3[0], vz_3[1], vz_3[2], vz_3[3], vz_3[4], vz_3[5], vz_3[6], vz_3[7]);
	printf("Usando s = %.12f, t = %.12f e r = %.12f\n\n", s, t, r);
	
	x = coordenada_fisica_3D(s, t, r, vx_3);
	y = coordenada_fisica_3D(s, t, r, vy_3);
	z = coordenada_fisica_3D(s, t, r, vz_3);

	i = coordenada_parametrica_3D(x, y, z, vx_3, vy_3, vz_3, &s, &t, &r, tol);
	
	printf("Resultado:\n");
	printf("x: %.12f\ny: %.12f\nz: %.12f\n", x, y, z);
	printf("s: %.22f\nt: %.22f\nr: %.22f\n", s, t, r);
	printf("iterações: %d\n", i);

	return 0;
}

double coordenada_fisica_2D(double s, double t, double *v)
{
	double soma = 0.0;
	double c = 1.0/4.0;
	double *N = (double *)malloc(4*sizeof(double));  // vetor com as funções de norma

	if (N == NULL) {
		puts("Erro de alocacao na memoria.");
		exit(1);
	}

	N[0] = c*(1-s)*(1-t);
	N[1] = c*(1-s)*(1+t);
	N[2] = c*(1+s)*(1+t);
	N[3] = c*(1+s)*(1-t);

	for (int i=0; i<4; i++)
		soma += N[i]*v[i];

	free(N);

	return soma;
}

double coordenada_fisica_3D(double s, double t, double r, double *v)
{
	double soma = 0.0;
	double c = 1.0/8.0;
	double *N = (double *)malloc(8*sizeof(double));  // vetor com as funções de norma

	if (N == NULL) {
		puts("Erro de alocacao na memoria.");
		exit(1);
	}

	N[0] = c*(1-s)*(1-t)*(1-r);
	N[1] = c*(1-s)*(1+t)*(1-r);
	N[2] = c*(1+s)*(1+t)*(1-r);
	N[3] = c*(1+s)*(1-t)*(1-r);
	N[4] = c*(1-s)*(1-t)*(1+r);
	N[5] = c*(1-s)*(1+t)*(1+r);
	N[6] = c*(1+s)*(1+t)*(1+r);
	N[7] = c*(1+s)*(1-t)*(1+r);

	for (int i=0; i<8; i++)
		soma += N[i]*v[i];

	free(N);

	return soma;
}

double raiz_2D(double p, double (*coord_fis_2D) (double s, double t, double *v), double s, double t, double *v)
{
	return (p - coord_fis_2D(s, t, v));
}

double raiz_3D(double p, double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v)
{
	return (p - coord_fis_3D(s, t, r, v));
}

double dps_2(double p, double (*raiz) (double p, double (*coord_fis_2D) (double s, double t, double *v), double s, double t, double *v), double (*coord_fis_2D) (double s, double t, double *v), double s, double t, double *v)
{
	return (raiz(p, coord_fis_2D, s+H, t, v)-raiz(p, coord_fis_2D, s, t, v))/H;
}

double dpt_2(double p, double (*raiz) (double p, double (*coord_fis_2D) (double s, double t, double *v), double s, double t, double *v), double (*coord_fis_2D) (double s, double t, double *v), double s, double t, double *v)
{
	return (raiz(p, coord_fis_2D, s, t+H, v)-raiz(p, coord_fis_2D, s, t, v))/H;
}

double dps_3(double p, double (*raiz) (double p, double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v), double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v)
{
	return (raiz(p, coord_fis_3D, s+H, t, r, v)-raiz(p, coord_fis_3D, s, t, r, v))/H;
}

double dpt_3(double p, double (*raiz) (double p, double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v), double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v)
{
	return (raiz(p, coord_fis_3D, s, t+H, r, v)-raiz(p, coord_fis_3D, s, t, r, v))/H;
}

double dpr_3(double p, double (*raiz) (double p, double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v), double (*coord_fis_3D) (double s, double t, double r, double *v), double s, double t, double r, double *v)
{
	return (raiz(p, coord_fis_3D, s, t, r+H, v)-raiz(p, coord_fis_3D, s, t, r, v))/H;
}

void sist_linear(int n, double **A, double *b, double *x)
{
	int *p = fatoracao(n, A);

	free(x);

	x = substituicao(n, A, p, b);

	free(p);
}

int coordenada_parametrica_2D(double x, double y, double *vx, double *vy, double *s, double *t, double tol)
{
	double *vf = (double *)malloc(2*sizeof(double));  // vetor com as funções u e v definidas no enunciado
	double *vh = (double *)malloc(2*sizeof(double));  // vetor passo a ser somado com a solução
	double *vs = (double *)malloc(2*sizeof(double));  // vetor solução
	double **J = mat_cria(2, 2);  // matriz Jacobiana

	int i = 0;

	if (vf == NULL || vh == NULL || vs == NULL) {
		puts("Erro de alocacao na memoria.");
		exit(1);
	}

	vs[0] = vs[1] = 0;

	while (fabs(x - coordenada_fisica_2D(vs[0], vs[1], vx)) >= tol || fabs(y - coordenada_fisica_2D(vs[0], vs[1], vy)) >= tol) {
		vf[0] = -raiz_2D(x, coordenada_fisica_2D, vs[0], vs[1], vx);  // função u(s,t)
		vf[1] = -raiz_2D(y, coordenada_fisica_2D, vs[0], vs[1], vy);  // função v(s,t)

		J[0][0] = dps_2(x, raiz_2D, coordenada_fisica_2D, vs[0], vs[1], vx);  // derivada parcial de u em s
		J[0][1] = dpt_2(x, raiz_2D, coordenada_fisica_2D, vs[0], vs[1], vx);  // derivada parcial de u em t
		J[1][0] = dps_2(y, raiz_2D, coordenada_fisica_2D, vs[0], vs[1], vy);  // derivada parcial de v em s
		J[1][1] = dpt_2(y, raiz_2D, coordenada_fisica_2D, vs[0], vs[1], vy);  // derivada parcial de v em t

		sist_linear(2, J, vf, vh);

		vs[0] += vh[0];
		vs[1] += vh[1];
		
		i++;
	}

	*s = vs[0];
	*t = vs[1];

	free(vf);
	free(vh);
	free(vs);
	mat_libera(2, J);

	return i;
}

int coordenada_parametrica_3D(double x, double y, double z, double *vx, double *vy, double *vz, double *s, double *t, double *r, double tol)
{
	double *vf = (double *)malloc(3*sizeof(double));  // vetor com as funções u, v e w
	double *vh = (double *)malloc(3*sizeof(double));  // vetor passo a ser somado com a solução
	double *vs = (double *)malloc(3*sizeof(double));  // vetor solução
	double **J = mat_cria(3, 3);  // matriz Jacobiana

	int i = 0;

	if (vf == NULL || vh == NULL || vs == NULL) {
		puts("Erro de alocacao na memoria.");
		exit(1);
	}

	vs[0] = vs[1] = vs[2] = 0;

	while (fabs(x - coordenada_fisica_3D(vs[0], vs[1], vs[2], vx)) >= tol || fabs(y - coordenada_fisica_3D(vs[0], vs[1], vs[2], vy)) >= tol || fabs(z - coordenada_fisica_3D(vs[0], vs[1], vs[2], vz)) >= tol) {
		vf[0] = -raiz_3D(x, coordenada_fisica_3D, vs[0], vs[1], vs[2], vx);  // função u(s,t,r)
		vf[1] = -raiz_3D(y, coordenada_fisica_3D, vs[0], vs[1], vs[2], vy);  // função v(s,t,r)
		vf[2] = -raiz_3D(z, coordenada_fisica_3D, vs[0], vs[1], vs[2], vz);  // função w(s,t,r)

		J[0][0] = dps_3(x, raiz_3D, coordenada_fisica_3D, vs[0], vs[1], vs[2], vx);  // derivada parcial de u em s
		J[0][1] = dpt_3(x, raiz_3D, coordenada_fisica_3D, vs[0], vs[1], vs[2], vx);  // derivada parcial de u em t
		J[0][2] = dpr_3(x, raiz_3D, coordenada_fisica_3D, vs[0], vs[1], vs[2], vx);  // derivada parcial de u em r
		J[1][0] = dps_3(y, raiz_3D, coordenada_fisica_3D, vs[0], vs[1], vs[2], vy);  // derivada parcial de v em s
		J[1][1] = dpt_3(y, raiz_3D, coordenada_fisica_3D, vs[0], vs[1], vs[2], vy);  // derivada parcial de v em t
		J[1][2] = dpr_3(y, raiz_3D, coordenada_fisica_3D, vs[0], vs[1], vs[2], vy);  // derivada parcial de v em r
		J[2][0] = dps_3(z, raiz_3D, coordenada_fisica_3D, vs[0], vs[1], vs[2], vz);  // derivada parcial de w em s
		J[2][1] = dpt_3(z, raiz_3D, coordenada_fisica_3D, vs[0], vs[1], vs[2], vz);  // derivada parcial de w em t
		J[2][2] = dpr_3(z, raiz_3D, coordenada_fisica_3D, vs[0], vs[1], vs[2], vz);  // derivada parcial de w em r

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
	mat_libera(3, J);

	return i;
}