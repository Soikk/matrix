#pragma once
#ifndef MATRIX_H
#define MATRIX_H

#define __MINGW_FEATURES__ 1

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


typedef struct matrix{
	int rows;
	int cols;
	long double **data;
} matrix;

matrix *newMatrix(int rows, int cols);

void freeMatrix(matrix **m);

matrix *copyMatrix(matrix *m);

matrix *identityMatrix(int order);

matrix *fillMatrix(matrix *m, long double n);

matrix *addMatrix(matrix *m, long double n);

matrix *subtractMatrix(matrix *m, long double n);

matrix *multiplyMatrix(matrix *m, long double n);

matrix *addMatrices(matrix *m1, matrix *m2);

matrix *subtractMatrices(matrix *m1, matrix *m2);

matrix *multiplyMatrices(matrix *m1, matrix *m2);

matrix *subMatrix(matrix *m, int row, int col);

long double determinant(matrix *m);

matrix *cofactor(matrix *m);

matrix *transpose(matrix *m);

matrix *adjugate(matrix *m);

matrix *inverse(matrix *m);

matrix *raiseMatrix(matrix *m, int n);

#endif
