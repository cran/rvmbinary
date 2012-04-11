#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "matrix.h"


void print_matrix( char* desc, int m, int n, double* a, int lda ) {
	int i, j;
	Rprintf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) Rprintf( " %6.2f", a[i+j*lda] );
		Rprintf( "\n" );
	}
}


void vprod(const matrix &a,const matrix &b, matrix &y, int trans){
	
	double alpha=1.0,beta=0.0;
	int inc=1;
	int rows=a.rows;
	int cols=a.cols;
	
	if (trans==0){
		y.reset(rows,1);
		F77_NAME(dgemv)("N",&rows,&cols,&alpha,a.data,&rows,b.data,&inc,&beta,y.data,&inc);
	}
	else {
		y.reset(cols,1);
		F77_NAME(dgemv)("T",&rows,&cols,&alpha,a.data,&rows,b.data,&inc,&beta,y.data,&inc);
	}
	
	
}

void mprod(const matrix &A,const matrix &B, matrix &y, int trans, double alpha){
	
	double beta=0.0;
	int nrows=A.rows;
	int ncols=A.cols;
	int bcols=B.cols;
	int brows=B.rows;
	
	int lda=A.rows;
	int ldb=B.rows;
	
	if (trans==0) {
		y.reset(nrows,bcols);
		F77_NAME(dgemm)("N", "N", &nrows, &bcols, &ncols, &alpha, A.data, &lda, B.data, &ldb, &beta, y.data, &nrows);
		
	}
	else if(trans==1){
		y.reset(ncols,bcols);
		F77_NAME(dgemm)("T", "N", &ncols, &bcols, &brows, &alpha, A.data, &lda, B.data, &ldb, &beta, y.data, &ncols);
	}
	else if(trans==2){
		y.reset(nrows,brows);
		F77_NAME(dgemm)("N", "T", &nrows, &brows, &ncols, &alpha, A.data, &lda, B.data, &ldb, &beta, y.data, &nrows);
	}
	else if(trans==3){
		y.reset(ncols,brows);
		F77_NAME(dgemm)("T","T", &ncols, &brows, &nrows, &alpha, A.data, &lda, B.data, &ldb, &beta, y.data, &ncols);
	}
	
}

void linalg(const matrix &a,const matrix &b,matrix &c,int trans){
	
	
	matrix Atemp(a.rows,a.cols,a.data);
	matrix Btemp(b.rows,b.cols,b.data);
	
	
	int m=a.rows;
	int n=a.cols;
	int nrhs=b.cols;
	int lda=a.rows;
	int ldb=b.rows;
	int lwork;
	double wkopt;
	int info;
	int i;
	double *work;
	lwork = -1;
	if (trans==0){
		F77_NAME(dgels)( "N", &m, &n, &nrhs, Atemp.data, &lda, Btemp.data, &ldb, &wkopt, &lwork,&info );
		lwork = (int)wkopt;
		work = (double*)malloc( lwork*sizeof(double) );
		F77_NAME(dgels)("N", &m, &n, &nrhs,  Atemp.data, &lda,  Btemp.data, &ldb,  work, &lwork, &info);
		c.reset(a.cols,1);
		for( i = 0; i < a.cols; i++ ) {
			c.data[i]=Btemp.data[i];
		}
	}
	else {
		F77_NAME(dgels)( "T", &n, &m, &nrhs, Atemp.data, &lda, Btemp.data, &ldb, &wkopt, &lwork,&info );
		lwork = (int)wkopt;
		work = (double*)malloc( lwork*sizeof(double) );
		F77_NAME(dgels)("T", &n, &m, &nrhs,  Atemp.data, &lda,  Btemp.data, &ldb,  work, &lwork, &info);
		c.reset(a.rows, 1);
		for (i=0; i<a.rows; i++) {
			c.data[i]=Btemp.data[i];
		}
	}

	free(work);
	if( info > 0 ) {
		printf( "The diagonal element %i of the triangular factor ", info );
		printf( "of A is zero, so that A does not have full rank;\n" );
		printf( "the least squares solution could not be computed.\n" );
		exit( 1 );
	}
	
	
	
}


int chol(const matrix &A,matrix &B){
	
	int N=A.rows;
	char T='U';
	int info;
	
	B.reset(A.rows,A.cols);
	
	for (int i=0; i<A.rows*A.cols; i++) {
			B.data[i]=A.data[i];
	}
	
	
	F77_NAME(dpotrf)(&T, &N, B.data, &N, &info);
	
	for (int i=1; i<N; i++) {
		for (int k=0; k<i; k++) {
			B.data[k*B.rows+i]=0.0;
		}
	}
	
	return(info);
}


int inverse(const matrix &A, matrix&B){
	
	int M=A.rows;
	int N=A.cols;
	
	int *IPIV;
	IPIV = (int*)malloc( A.rows*sizeof(int) );
	
	B.reset(A.rows,A.cols);
	
	for (int i=0; i<A.rows; i++) {
		for (int k=0; k<A.cols; k++) {
			B.data[(i*A.cols)+k]=A.data[(i*A.cols)+k];
		}
	}
	
	int info;
	double *work;
	int lwork=A.rows*A.cols;
	work = (double*)malloc( lwork*sizeof(double) );

	
	dgetrf_(&M, &N, B.data, &M, IPIV, &info);
	dgetri_(&N, B.data, &M, IPIV, work, &lwork, &info);
	
	free(IPIV);
	free(work);
	
	return(info);
}



