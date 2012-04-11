/*
 *  fullstatistics.h
 *  RVM-Speed
 *
 *  Created by Robert Lowe on 27/08/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "matrix.h"
#include <vector>
#include <math.h>
#include "matrixops.h"
#include <limits>

int fullstatistics(int likelihood,const matrix &PHI,const matrix &BASIS,const matrix &BASIS2,matrix &beta,matrix& Sigma,matrix& Mu
					, matrix &Alpha,double &logML,const matrix &Targets,const std::vector<int> &Used,matrix& Factor,
					matrix &S_out, matrix &Q_in,matrix &S_in, matrix &Q_out,matrix &betaBASIS_PHI,matrix &BASIS_PHI,matrix &BASIS_Targets,matrix &Gamma);

int PosteriorMode(matrix &U,const matrix &BASIS,matrix &beta,const matrix &Targets, const matrix &Alpha,matrix &Mu,int itsMax,int likelihood, double &dataLikely);
double DataError(int likelihood,const matrix &BASIS_Mu,const matrix &Targets, matrix &y);
void Sigmoid(const matrix &A,matrix&B);
