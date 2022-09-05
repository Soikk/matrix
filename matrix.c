#include "matrix.h"

static int global = 0;
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

void saveMatrix(matrix *m, FILE *fp){
	char header = 'M';
	fwrite(&header, sizeof(char), 1, fp);
	fwrite(&m->rows, sizeof(int), 1, fp);
	fwrite(&m->cols, sizeof(int), 1, fp);
	for(int i = 0; i < m->rows; ++i){
		fwrite(m->data[i], sizeof(long double), m->cols, fp);
	}
	char end = 'E';
	fwrite(&end, sizeof(char), 1, fp);
}

matrix *loadMatrix(FILE *fp){
	char header;
	fread(&header, sizeof(char), 1, fp);
	if(header != 'M'){
		fprintf(stderr, "Header is '%c' not 'M'\n", header);
		exit(EXIT_FAILURE);
	}
	int rows, cols;
	fread(&rows, sizeof(int), 1, fp);
	fread(&cols, sizeof(int), 1, fp);
	matrix *m = malloc(sizeof(matrix));
	m->rows = rows;
	m->cols = cols;
	m->data = malloc(rows*sizeof(long double*));
	for(int i = 0; i < m->rows; ++i){
		m->data[i] = malloc(cols*sizeof(long double));
		fread(m->data[i], sizeof(long double), cols, fp);
	}
	char end;
	fread(&end, sizeof(char), 1, fp);
	if(end != 'E'){
		fprintf(stderr, "End is '%c' not 'E'\n", end);
		exit(EXIT_FAILURE);
	}
	return m;
}

static inline bool sameDimensions(matrix *a, matrix *b){
	return (a->rows == b->rows) && (a->cols == b->cols);
}

static inline bool isSquare(matrix *m){
	return m->rows == m->cols;
}

void copyMatrix(matrix *dest, matrix *src){
	if(!sameDimensions(dest, src)){
		fprintf(stderr, "Wrong dimensions (%dx%d != %dx%d)\n", dest->rows, dest->cols, src->rows, src->cols);
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i < src->rows; ++i){
		for(int j = 0; j < src->cols; ++j){
			dest->data[i][j] = src->data[i][j];
		}
	}
}

matrix *cloneMatrix(matrix *m){
	matrix *r = newMatrix(m->rows, m->cols);
	copyMatrix(r, m);
	return r;
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

void fillMatrix(matrix *m, long double n){
	for(int i = 0; i < m->rows; ++i){
		for(int j = 0; j < m->cols; ++j){
			m->data[i][j] = n;
		}
	}
}

void addMatrix(matrix *m, long double n){
	for(int i = 0; i < m->rows; ++i){
		for(int j = 0; j < m->cols; ++j){
			m->data[i][j] += n;
		}
	}
}

void subtractMatrix(matrix *m, long double n){
	addMatrix(m, -n);
}

void multiplyMatrix(matrix *m, long double n){
	for(int i = 0; i < m->rows; ++i){
		for(int j = 0; j < m->cols; ++j){
			m->data[i][j] *= n;
		}
	}
}

matrix *addMatrices(matrix *a, matrix *b){
	if(!sameDimensions(a, b)){
		fprintf(stderr, "Wrong dimensions (%dx%d != %dx%d)\n", a->rows, a->cols, b->rows, b->cols);
		exit(EXIT_FAILURE);
	}
	matrix *r = newMatrix(a->rows, a->cols);
	for(int i = 0; i < r->rows; ++i){
		for(int j = 0; j < r->cols; ++j){
			r->data[i][j] = a->data[i][j] + b->data[i][j];
		}
	}
	return r;
}

matrix *subtractMatrices(matrix *a, matrix *b){
	matrix *t = cloneMatrix(b);
	multiplyMatrix(t, -1);
	matrix *r = addMatrices(a, t);
	freeMatrix(&t);
	return r;
}

matrix *multiplyMatrices(matrix *a, matrix *b){
	if(a->cols != b->rows){
		fprintf(stderr, "Wrong dimensions (%dx%d != %dx%d)\n", a->rows, a->cols, b->rows, b->cols);
		exit(EXIT_FAILURE);
	}
	matrix *r = newMatrix(a->rows, b->cols);
	for(int i = 0; i < r->rows; ++i){
		for(int j = 0; j < r->cols; ++j){
			long double sum = 0;
			for(int n = 0; n < a->cols; ++n){
				sum += a->data[i][n] * b->data[n][j];
			}
			r->data[i][j] = sum;
		}
	}
	return r; 
}

matrix *HadamardProduct(matrix *a, matrix *b){
	if(!sameDimensions(a, b)){
		fprintf(stderr, "Wrong dimensions (%dx%d != %dx%d)\n", a->rows, a->cols, b->rows, b->cols);
		exit(EXIT_FAILURE);
	}
	matrix *r = newMatrix(a->rows, a->cols);
	for(int i = 0; i < r->rows; ++i){
		for(int j = 0; j < r->cols; ++j){
			r->data[i][j] = a->data[i][j] * b->data[i][j];
		}
	}
	return r;
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
			freeMatrix(&s);
			d += !(i%2) ? v : -v;
		}
		return d;
	}
}

void cofactor(matrix *m){
	if(!isSquare(m)){
		fprintf(stderr, "Matrix is not square (%dx%d)\n", m->rows, m->cols);
		exit(EXIT_FAILURE);
	}
	matrix *r = newMatrix(m->rows, m->cols);
	for(int i = 0; i < r->rows; ++i){
		for(int j = 0; j < r->cols; ++j){
			matrix *s = subMatrix(m, i, j);
			long double ds = determinant(s);
			freeMatrix(&s);
			r->data[i][j] = !((i+j)%2) ? ds : -ds;
		}
	}
	copyMatrix(m, r);
	freeMatrix(&r);
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

matrix *dotProduct(matrix *a, matrix *b){
	matrix *tb = transpose(b);
	matrix *r = multiplyMatrices(a, tb);
	freeMatrix(&tb);
	return r;
}

void adjugate(matrix *m){
	cofactor(m);
	matrix *t = transpose(m);
	copyMatrix(m, t);
	freeMatrix(&t);
}

void invert(matrix *m){
	long double d = determinant(m);
	if(d == 0){
		fprintf(stderr, "Determinant is 0, the matrix is not invertible\n");
		exit(EXIT_FAILURE);
	}
	adjugate(m);
	multiplyMatrix(m, 1/d);
}

void raiseMatrix(matrix *m, int n){
	if(!isSquare(m)){
		fprintf(stderr, "Matrix is not square (%dx%d)\n", m->rows, m->cols);
		exit(EXIT_FAILURE);
	}
	if(n < 0){
		invert(m);
		n = -n;
	}
	matrix *r = identityMatrix(m->rows);
	for(int i = 0; i < n; ++i){
		matrix *t = multiplyMatrices(r, m);
		copyMatrix(r, t);
		freeMatrix(&t);
	}
	copyMatrix(m, r);
	freeMatrix(&r);
}


/*int main(){

	int input = 784, hidden = 300, output = 10;
	
	matrix *input_matrix = newMatrix(input, 1); input_matrix = fillMatrix(input_matrix, 1);
	matrix *hidden_weights = newMatrix(hidden, input); hidden_weights = fillMatrix(hidden_weights, 1);
	matrix *output_weights = newMatrix(output, hidden); output_weights = fillMatrix(output_weights, 1);
	
	matrix *hidden_inputs = multiplyMatrices(hidden_weights, input_matrix);
	matrix *final_inputs = multiplyMatrices(output_weights, hidden_inputs);
	
	printf("done\n");
	
	matrix in = newMatrix(1, 2); in.data[0][0] = 0; in.data[0][1] = 1;
	matrix layer = newMatrix(1, 2); layer.data[0][0] = 0.188; layer.data[0][1] = 0.812;
	matrix layerm = newMatrix(1, 2); layerm.data[0][0] = 0.812; layerm.data[0][1] = 0.188;
	
	matrix diff = subtractMatrices(in, layer);
	matrix der = HadamardProduct(layer, layerm);
	matrix res = HadamardProduct(diff, der);
	
	for(int i = 0; i < res.rows; ++i){
		for(int j = 0; j < res.cols; ++j){
			printf("%.3Lf\t", res.data[i][j]);
		}
		printf("\n");
	}

	return 0;
}*/
