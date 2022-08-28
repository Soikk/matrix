#include "matrix.h"


matrix *newMatrix(int rows, int cols){
	matrix *m = malloc(sizeof(matrix));
	m->rows = rows;
	m->cols = cols;
	m->data = malloc(rows*sizeof(long double*));
	for(int i = 0; i < rows; ++i){
		m->data[i] = malloc(cols*sizeof(long double));
	}
	return m;
}

void freeMatrix(matrix **m){
	for(int i = 0; i < (*m)->rows; ++i){
		free((*m)->data[i]);
		(*m)->data[i] = NULL;
	}
	free((*m)->data);
	(*m)->data = NULL;
	free(*m);
	*m = NULL;
}

matrix *copyMatrix(matrix *m){
	matrix *c = newMatrix(m->rows, m->cols);
	for(int i = 0; i < m->rows; ++i){
		for(int j = 0; j < m->cols; ++j){
			c->data[i][j] = m->data[i][j];
		}
	}
	return c;
}

matrix *identityMatrix(int order){
	matrix *im = newMatrix(order, order);
	for(int i = 0; i < im->rows; ++i){
		for(int j = 0; j < im->cols; ++j){
			im->data[i][j] = (i == j);
		}
	}
	return im;
}

matrix *fillMatrix(matrix *m, long double n){
	matrix *r = newMatrix(m->rows, m->cols);
	for(int i = 0; i < r->rows; ++i){
		for(int j = 0; j < r->cols; ++j){
			r->data[i][j] = n;
		}
	}
	return r;
}

matrix *addMatrix(matrix *m, long double n){
	for(int i = 0; i < m->rows; ++i){
		for(int j = 0; j < m->cols; ++j){
			m->data[i][j] += n;
		}
	}
}

matrix *subtractMatrix(matrix *m, long double n){
	return addMatrix(m, -n);
}

matrix *multiplyMatrix(matrix *m, long double n){
	matrix *r = copyMatrix(m);
	for(int i = 0; i < r->rows; ++i){
		for(int j = 0; j < r->cols; ++j){
			r->data[i][j] *= n;
		}
	}
	return r;
}

static inline bool sameDimensions(matrix *m1, matrix *m2){
	return (m1->rows == m2->rows) && (m1->cols == m2->cols);
}

matrix *addMatrices(matrix *m1, matrix *m2){
	if(!sameDimensions(m1, m2)){
		fprintf(stderr, "Wrong dimensions (%dx%d != %dx%d)\n", m1->rows, m1->cols, m2->rows, m2->cols);
		return NULL;
	}
	matrix *r = newMatrix(m1->rows, m1->cols);
	for(int i = 0; i < m1->rows; ++i){
		for(int j = 0; j < m1->cols; ++j){
			r->data[i][j] = m1->data[i][j] + m2->data[i][j];
		}
	}
	return r;
}

matrix *subtractMatrices(matrix *m1, matrix *m2){
	m2 = multiplyMatrix(m2, (long double)-1);
	return addMatrices(m1, m2);
}

matrix *multiplyMatrices(matrix *m1, matrix *m2){
	if(m1->cols != m2->rows){
		fprintf(stderr, "Wrong dimensions (%dx%d != %dx%d)\n", m1->rows, m1->cols, m2->rows, m2->cols);
		return NULL;
	}
	matrix *r = newMatrix(m1->rows, m2->cols);
	for(int i = 0; i < r->rows; ++i){
		for(int j = 0; j < r->cols; ++j){
			long double sum = 0;
			for(int n = 0; n < m1->cols; ++n){
				sum += m1->data[i][n] * m2->data[n][j];
			}
			r->data[i][j] = sum;
		}
	}
	return r; 
}

matrix *HadamardProduct(matrix *m1, matrix *m2){
	if(!sameDimensions(m1, m2)){
		fprintf(stderr, "Wrong dimensions (%dx%d != %dx%d)\n", m1->rows, m1->cols, m2->rows, m2->cols);
		return NULL;
	}
	matrix *r = newMatrix(m1->rows, m1->cols);
	for(int i = 0; i < r->rows; ++i){
		for(int j = 0; j < r->cols; ++j){
			r->data[i][j] = m1->data[i][j] * m2->data[i][j];
		}
	}
	return r;
}

static inline bool isSquare(matrix *m){
	return m->rows == m->cols;
}

matrix *subMatrix(matrix *m, int row, int col){
	matrix *r = newMatrix(m->rows-1, m->cols-1);
	int ri = 0, rj = 0;
	for(int i = 0; i < m->rows; ++i){
		if(i == row){ continue; }
		rj = 0;
		for(int j = 0; j < m->cols; ++j){
			if(j == col){ continue; }
			r->data[ri][rj] = m->data[i][j];
			++rj;
		}
		++ri;
	}
	return r;
}

long double determinant(matrix *m){
	if(!isSquare(m)){
		fprintf(stderr, "Matrix is not square (%dx%d)\n", m->rows, m->cols);
		exit(EXIT_FAILURE);
	}
	long double d = 0;
	if(m->rows == 1){
		return m->data[0][0];
	}else{
		for(int i = 0; i < m->rows; ++i){
			matrix *s = subMatrix(m, 0, i);
			long double v = determinant(s)*m->data[0][i];
			d += !(i%2) ? v : -v;
		}
		return d;
	}
}

matrix *cofactor(matrix *m){
	if(!isSquare(m)){
		fprintf(stderr, "Matrix is not square (%dx%d)\n", m->rows, m->cols);
		return NULL;
	}
	matrix *r = newMatrix(m->rows, m->cols);
	for(int i = 0; i < r->rows; ++i){
		for(int j = 0; j < r->cols; ++j){
			long double v = determinant(subMatrix(m, i, j));
			r->data[i][j] = !((i+j)%2) ? v : -v;
		}
	}
	return r;
}

matrix *transpose(matrix *m){
	matrix *r = newMatrix(m->cols, m->rows);
	for(int i = 0; i < r->rows; ++i){
		for(int j = 0; j < r->cols; ++j){
			r->data[i][j] = m->data[j][i];
		}
	}
	return r;
}

matrix *dotProduct(matrix *m1, matrix *m2){
	return multiplyMatrices(m1, transpose(m2));
}

matrix *adjugate(matrix *m){
	return transpose(cofactor(m));
}

matrix *inverse(matrix *m){
	long double d = determinant(m);
	if(d == 0){
		fprintf(stderr, "Determinant is 0, the matrix is not invertible\n");
		return NULL;
	}
	return multiplyMatrix(adjugate(m), 1/d);
}

matrix *raiseMatrix(matrix *m, int n){
	if(!isSquare(m)){
		fprintf(stderr, "Matrix is not square (%dx%d)\n", m->rows, m->cols);
		return NULL;
	}
	if(n < 0){
		m = inverse(m);
		n = -n;
	}
	matrix *r = identityMatrix(m->rows);
	for(int i = 0; i < n; ++i){
		r = multiplyMatrices(r, m);
	}
	return r;
}


/*int main(){

	int input = 784, hidden = 300, output = 10;
	
	matrix *input_matrix = newMatrix(input, 1); input_matrix = fillMatrix(input_matrix, 1);
	matrix *hidden_weights = newMatrix(hidden, input); hidden_weights = fillMatrix(hidden_weights, 1);
	matrix *output_weights = newMatrix(output, hidden); output_weights = fillMatrix(output_weights, 1);
	
	matrix *hidden_inputs = multiplyMatrices(hidden_weights, input_matrix);
	matrix *final_inputs = multiplyMatrices(output_weights, hidden_inputs);
	
	printf("done\n");

	return 0;
}*/
