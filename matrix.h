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

void copyMatrix(matrix *dest, matrix *src);

matrix *cloneMatrix(matrix *m);

matrix *identityMatrix(int order);

void fillMatrix(matrix *m, long double n);

void addMatrix(matrix *m, long double n);

void subtractMatrix(matrix *m, long double n);

void multiplyMatrix(matrix *m, long double n);

matrix *addMatrices(matrix *a, matrix *b);

matrix *subtractMatrices(matrix *a, matrix *b);

matrix *multiplyMatrices(matrix *a, matrix *b);

matrix *HadamardProduct(matrix *a, matrix *b);

matrix *subMatrix(matrix *m, int row, int col);

long double determinant(matrix *m);

matrix *cofactor(matrix *m);

matrix *transpose(matrix *m);

matrix *dotProduct(matrix *a, matrix *b);

matrix *adjugate(matrix *m);

void invert(matrix *m);

void raiseMatrix(matrix *m, int n);

#endif
