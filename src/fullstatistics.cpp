/*
 *  fullstatistics.cpp
 *  RVM-Speed
 *
 *  Created by Robert Lowe on 27/08/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "fullstatistics.h"
#include <R.h>

using namespace std;

void fullstatistics(int likelihood,const matrix &PHI,const matrix &BASIS,const matrix &BASIS2,matrix &beta,matrix& Sigma,matrix& Mu
				   , matrix &Alpha,double &logML,const matrix &Targets,const std::vector<int> &Used,matrix& Factor
				   ,matrix &S_out, matrix &Q_in,matrix &S_in, matrix &Q_out,matrix &betaBASIS_PHI,matrix &Gamma){

	int MAX_POSTMODE_ITS=25;
	//Check these are the right way round
	int N=BASIS.rows;
	int M_FULL=BASIS.cols;
	int M=PHI.cols;
	matrix y(Targets.rows,1),e(Targets.rows,1);
	
	matrix U;
	double dataLikely;
	matrix Ui;

	if (likelihood==0){
		//Gaussian Case
	}
	else{

		PosteriorMode(U,PHI,beta,Targets,Alpha,Mu,MAX_POSTMODE_ITS,likelihood,dataLikely);

		if (U.rows==1 and U.cols==1){
			Ui.reset(1, 1);	
			Ui.data[0]=1.0/U.data[0];
		}
		else{
			Ui.reset(U.rows, U.cols);
			inverse(U, Ui);
		}
		mprod(Ui,Ui,Sigma,2,1.0);
		if (likelihood==1){
			matrix temp;
			mprod(PHI, Mu, temp, 0, 1.0);
			Sigmoid(temp, y);
		}
		for (int i=0; i<y.rows; i++) {
			e.data[i]=Targets.data[i]-y.data[i];
		}
		
		double logdetHOver2=0.0;
		for (int i=0; i<U.rows; i++){
			logdetHOver2+=log(U.data[i*U.cols+i]);
		}
		
		double sumlogalpha=0.0;
		
		
		for (int i=0; i<Alpha.rows; i++){
			sumlogalpha+=log(Alpha.data[i]);
		}
		
		matrix temp(Mu.rows,Mu.cols);
		for (int i=0; i<temp.rows; i++){
			for(int k=0; k<temp.cols; k++){
				temp.data[i+k*temp.rows]=pow(Mu.data[i+k*temp.rows],2);
			}
		}
		
		matrix temp2;
		mprod(temp, Alpha, temp2, 1, 1.0);
		logML=dataLikely-temp2.data[0]/2.0+sumlogalpha/2.0-logdetHOver2;
		
		matrix DiagC(Ui.rows,1);
		for (int i=0; i<Ui.rows; i++) {
			DiagC.data[i]=0.0;
			for(int k=0; k<Ui.cols; k++){
				DiagC.data[i]+=(Ui.data[i+k*Ui.rows]*Ui.data[i+k*Ui.rows]);
			}
		}
		Gamma.reset(Ui.rows,1);
		for(int i=0; i<Gamma.rows; i++){
			Gamma.data[i]=1-Alpha.data[i]*DiagC.data[i];
		}		

		//NON GAUSSIAN CASE!!!!
		matrix temp3(beta.rows,M);		
		for(int i=0; i<beta.rows; i++){
			for (int k=0; k<M; k++){
				temp3.data[i+temp3.rows*k]=PHI.data[i+temp3.rows*k]*beta.data[i];
			}
		}
		//cout << "3" <<endl;
		
		mprod(BASIS,temp3,betaBASIS_PHI,1,1.0);
		
		//cout << "4" <<endl;
		
		//BASIS2 removed outside of loop 1/09/10


		matrix bb_phi_Ui;
		mprod(betaBASIS_PHI,Ui,bb_phi_Ui,0,1.0);
		
		matrix sumbb_phi_Ui(bb_phi_Ui.rows,1);

		
		for (int i=0; i<bb_phi_Ui.rows; i++) {
			sumbb_phi_Ui.data[i]=0.0;
			for(int k=0; k<bb_phi_Ui.cols; k++){
				sumbb_phi_Ui.data[i]+=(bb_phi_Ui.data[i+bb_phi_Ui.rows*k]*bb_phi_Ui.data[i+bb_phi_Ui.rows*k]);
			}
		}

		
		mprod(beta, BASIS2, S_in, 1, 1.0);
		

		S_in.rows=S_in.cols;
		S_in.cols=1;
		
		
		S_out.reset(S_in.rows,S_in.cols);
		
		for (int i=0; i<S_in.rows; i++) {
				S_in.data[i]-=sumbb_phi_Ui.data[i];
				S_out.data[i]=S_in.data[i];
			}

		mprod(BASIS,e,Q_in,1,1.0);
		
		Q_out.reset(Q_in.rows,Q_in.cols);
		
		for (int i=0; i<Q_in.rows; i++) {
				Q_out.data[i]=Q_in.data[i];
		}

		for (int i=0; i<Used.size(); i++){
			S_out.data[Used[i]]= (Alpha.data[i]*S_in.data[Used[i]]/(Alpha.data[i]-S_in.data[Used[i]]));
			Q_out.data[Used[i]]= (Alpha.data[i]*Q_in.data[Used[i]]/(Alpha.data[i]-S_in.data[Used[i]]));
		}
		
		Factor.reset(S_out.rows,1);
		for(int i=0; i<Factor.rows; i++){
			Factor.data[i]=(Q_out.data[i]*Q_out.data[i])-S_out.data[i];
		}
		 
	}

}
	
void PosteriorMode(matrix &U,const matrix &BASIS,matrix &beta, const matrix &Targets,const matrix &Alpha,matrix &Mu,int itsMax,int likelihood, double &dataLikely){
	
	double GRADIENT_MIN=0.000001;
	double STEP_MIN=1/(pow(2.0,8.0));
	
	int N=BASIS.rows;
	int M=BASIS.cols;
	
	
	matrix BASIS_Mu;
	mprod(BASIS,Mu,BASIS_Mu,0,1.0);
	matrix y;
	
	double dataError=DataError(likelihood,BASIS_Mu,Targets,y);
	
	double regulariser=0.0;
	for (int i=0; i<Mu.rows; i++){
		regulariser+=(Alpha.data[i]*(Mu.data[i]*Mu.data[i]))/2.0;
	}
	double newTotalError=dataError+regulariser;
	matrix errorLog(itsMax,1);

	for (int iteration=0; iteration<itsMax; iteration++){
		errorLog.data[iteration]=newTotalError;
		//Rprintf("PM cycle: %2d\t error: %.6f\n",iteration, errorLog.data[iteration]);
		matrix tempAMu(Alpha.rows,Alpha.cols);
		
		matrix e(y.rows,1);
		for (int i=0; i<y.rows; i++) {
			e.data[i]=Targets.data[i]-y.data[i];
			beta.data[i]=y.data[i]*(1-y.data[i]);
		}
		
		for (int i=0; i<Alpha.rows; i++){
			tempAMu.data[i]=Alpha.data[i]*Mu.data[i];
		}
		
		matrix g;
		mprod(BASIS, e, g, 1, 1.0);
		
		for (int i=0; i<g.rows; i++) {
			g.data[i]=g.data[i]-tempAMu.data[i];
		}
		
		matrix BASIS_B(beta.rows,M);		
		for(int i=0; i<beta.rows; i++){
			for (int k=0; k<M; k++){
				BASIS_B.data[i+beta.rows*k]=BASIS.data[i+beta.rows*k]*beta.data[i];
			}
		}
		
		matrix H;
		mprod(BASIS_B, BASIS, H, 1, 1.0);
		
		for (int i=0; i<Alpha.rows; i++) {
			H.data[i+i*H.rows]+=Alpha.data[i];
		}
		
		int error=0;
		error=chol(H, U);
		
		if (error!=0){
			Rprintf("***Warning ** Ill conditioned Hessian (%f)\n",error);
			exit(1);
		}
		
		bool gtest=1;
		double maxg=0.0;
		for(int i=0; i<g.rows; i++){
			if(fabs(g.data[i])>maxg)
				maxg=fabs(g.data[i]);
			if (fabs(g.data[i])>GRADIENT_MIN) {
				gtest=0;
				break;
			}
		}
		
		if (gtest==1){
			//Is therean easier way to select part of matrix
			matrix temperrorLog(iteration,1);
			for (int i=0; i<iteration; i++) {
				temperrorLog.data[i]=errorLog.data[i];
			}
			//errorLog=temperrorLog;
			//ERROR LOG INCOMPLETE
			//Rprintf("PM convergence (<%g) after %d iterations, |g| = %g\n",GRADIENT_MIN,iteration,maxg);
			break;
		}
		matrix DeltaMu(g.rows,g.cols);
		
		if (U.rows==1 and U.cols==1 and g.rows==1 and g.cols==1)
			DeltaMu.data[0]=(g.data[0]/U.data[0])/U.data[0];
		else{
			matrix temp;
			linalg(U, g, temp, 1);
			linalg(U, temp, DeltaMu, 0);

		}
		
		


		double step=1;
		while (step>STEP_MIN) {
			matrix Mu_new(Mu.rows,1);
			for(int i=0; i<Mu.rows; i++){
					Mu_new.data[i]=Mu.data[i] + step*(DeltaMu.data[i]);
			}
			
			mprod(BASIS,Mu_new,BASIS_Mu,0,1.0);
			
			dataError=DataError(likelihood,BASIS_Mu,Targets,y);
			regulariser=0.0;
			for (int i=0; i<Mu_new.rows; i++){
				regulariser+=(Alpha.data[i]*(Mu_new.data[i]*Mu_new.data[i]))/2.0;
			}
			
			newTotalError=dataError+regulariser;
			if (newTotalError>=errorLog.data[iteration]) {
				step=step/2.0;
			}
			else {
				for(int i=0; i<Mu_new.rows; i++){
					Mu.data[i]=Mu_new.data[i];
				}
				step=0;
			}
		}
		if (step!=0){
			//Could print out for higher verbosity 
			Rprintf("ARE WE HERE \n");
			break;
		}
		
	}
	dataLikely=-dataError;
	
}

void Sigmoid(const matrix &A,matrix &B){
	
	B.reset(A.rows,A.cols);
	for (int i=0; i<A.rows; i++){
		for(int k=0; k<A.cols; k++){
			B.data[i+k*A.rows]=1.0/(1.0+exp(-A.data[i+k*A.rows]));
		}
	}

}

double DataError(int likelihood,const matrix &BASIS_Mu,const matrix &Targets, matrix &y){
	
	double e=0.0;
	
	if (likelihood==1){
		Sigmoid(BASIS_Mu,y);
		for (int i=0; i<Targets.rows; i++){
			
			if ((y.data[i]==1.0 and Targets.data[i]==0.0) or (y.data[i]==0.0 and Targets.data[i]==1.0)){
				e=std::numeric_limits<double>::infinity();
				break;
			}
			if(y.data[i]!=0.0){
				e+= -(Targets.data[i]*log(y.data[i]));
			}
			if(y.data[i]!=1.0){
				e+=-(1-Targets.data[i])*log(1-y.data[i]);
			}
			
		}
		
	}
	return(e);
}

