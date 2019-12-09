//
//  matriz.h
//  HexaCellLocal
//
//  Created by Felipe Viberti and Matheus Leal on 05/12/19.
//  Copyright © 2019 Felipe Viberti. All rights reserved.
//

double** matcria (int m, int n);

void matlibera (int m, double** A);

void transposta (int m, int n, double** A, double** T);

void multv (int m, int n, double** A, double* v, double* w);

void multm (int m, int n, int q, double** A, double** B, double** C);
