/*
 *  matrix.h
 *  RVM-Speed
 *
 *  Created by Robert Lowe on 27/08/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MATRIXOPS_H
#define MATRIXOPS_H

void print_matrix( char* desc, int m, int n, double* a, int lda );
void vprod(const matrix &a,const matrix &b, matrix &y, int trans);
void mprod(const matrix &A,const matrix &B, matrix &y, int trans, double alpha);
void linalg(const matrix &a, const matrix &b,matrix &c, int trans);
int inverse(const matrix &A, matrix&B);
int chol(const matrix &A,matrix &B);

#endif
