#include "matrix.h"
#include "matrixops.h"
#include "fullstatistics.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <R.h>
#include <vector>
#include <string>

extern "C" {
	
	
	void initial (double *a, double *b,int *arows,int *acols,int *PARAMATERrev, double *PARAMATERval, int * maxits) {
		
		int ItNum=maxits[0];
		double MinDeltaLogAlpha=1e-3,MinDeltaLogBeta = 1e-6,AlignmentMax=1-1e-3;
		

		
		bool PriorityAddition=0,PriorityDeletion=1,BasisAlignmentTest=1;
		//Reporting on the iterations
		int monitor_its=10;
		int rows=arows[0];
		int cols=acols[0];
		
		//Assume BASIS has number cols same as number of DATA points
		matrix BASIS(rows,cols,a);
		matrix Targets(cols,1,b);
		std::vector<int> Used;

		//***************************************************
		//***************** INITIALISE **********************
		//***************************************************
		
		Rprintf("Initialise\n");
		double GAUSSIAN_SNR_INIT	= 0.1;
		double	INIT_ALPHA_MAX	= 1e3;
		double INIT_ALPHA_MIN	= 1e-3;
		//Normalise each Basis vector
		matrix Scales(1,BASIS.cols);
		
		for (int k=0; k<BASIS.rows; k++){
			Scales.data[k]=0;
			for (int i=0; i<BASIS.cols; i++){
				Scales.data[k]+=pow(BASIS.data[(k*BASIS.rows)+i],2);
			}
			Scales.data[k]=sqrt(Scales.data[k]);
			if (Scales.data[k]==0){
				Scales.data[k]=1;
			}
		}
		
		for (int i=0; i<BASIS.rows; i++){
			for (int k=0; k<BASIS.cols; k++){
				BASIS.data[(k*BASIS.rows)+i]=BASIS.data[(k*BASIS.rows)+i]/Scales.data[k];
			}
		}
		
		matrix LogOut(arows[0],1);
		matrix TargetsPseudoLinear(arows[0],1);
		
		for (int i=0; i<arows[0]; i++){
			TargetsPseudoLinear.data[i]=(2*Targets.data[i]-1);
			LogOut.data[i]	= (TargetsPseudoLinear.data[i]*0.9+1)/2.0;
			LogOut.data[i]=log(LogOut.data[i]/(1-LogOut.data[i]));
		}
		
		matrix proj;
		vprod(BASIS,TargetsPseudoLinear,proj,1);

		double max=0.0;
		int maxindex=0;
		for(int i=0; i<proj.rows; i++){
			if (fabs(proj.data[i])>max) {
				max=fabs(proj.data[i]);
				maxindex=i;
			}
		}
		
		Rprintf("Initialising with maximally aligned basis vector (%d)\n",maxindex+1);
		matrix PHI,Mu,Alpha;
		PHI.AddColumn(BASIS, maxindex);
		Used.push_back(maxindex);
		linalg(PHI, LogOut,Mu,0);

		Alpha.data=new double[1];
		Alpha.rows=1;
		Alpha.cols=1;
		
		Rprintf("Initial Alpha = ");
		if (Mu.data[0]==0)
			Alpha.data[0]=1;
		else{
			double value=1/(Mu.data[0]*Mu.data[0]);
			if (value<INIT_ALPHA_MIN) 
				Alpha.data[0]=INIT_ALPHA_MIN;
			else if(value>INIT_ALPHA_MAX)
				Alpha.data[0]=INIT_ALPHA_MAX;
			else
				Alpha.data[0]=value;
		}
		Rprintf("%f\n", Alpha.data[0]);
		
		//***************************************************
		//**************** END OF INITIALISE ****************
		//***************************************************
		matrix BASIS2(BASIS.rows,BASIS.cols);
		for (int i=0; i<BASIS.rows; i++) {
			for (int k=0; k<BASIS.cols; k++) {
				BASIS2.data[i+BASIS.rows*k]=(BASIS.data[i+BASIS.rows*k]*BASIS.data[i+BASIS.rows*k]);
			}
		}
		matrix BASIS_PHI,BASIS_B_PHI;
		matrix BASIS_Targets;
		mprod(BASIS,Targets,BASIS_Targets,1,1.0);
		matrix S_out,Q_in,S_in,Q_out,Gamma;
		matrix beta(rows,1),SIGMA,Factor;
		double logML;
		
		fullstatistics(1, PHI, BASIS,BASIS2,beta,SIGMA,Mu,Alpha,logML,Targets,Used,Factor,S_out,Q_in,S_in,Q_out,BASIS_B_PHI,Gamma);
		
		int N=BASIS.rows;
		int M_full=BASIS.cols;
		int M=PHI.rows;
		
		int addCount=0;
		int deleteCount=0;
		int updateCount=0;
		
		//Control not present for betaupdatestart and BetaUpdateFrequency
		int maxLogSize=ItNum+10+(int)(ItNum/5);
		
		matrix logMarginalLog(maxLogSize,1);
		int count=0;
		
		std::vector<double> Aligned_out,Aligned_in;
		int alignDeferCount=0;
		
		//Action Codes
		const int ACTION_REESTIMATE=0;
		const int ACTION_ADD=1;
		const int ACTION_DELETE=-1;
		
		const int ACTION_TERMINATE=10;
		const int ACTION_NOISE_ONLY=11;
		const int ACTION_ALIGNMENT_SKIP=12;
		
		int selectedAction;
		/*
		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		 %%
		 %% MAIN LOOP
		 %%
		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		 */
		int iteration_counter=0;
		bool LAST_ITERATION=1;
		
		while (LAST_ITERATION) {
			
			iteration_counter+=1;
			//NEVER Update  iteration??
			bool UpdateIteration=0;
			/*
			 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			 
			 %% DECISION PHASE
			 %%
			 %% Assess all potential actions
			 %%
			 */
			std::vector<double> GoodFactor(M_full);
			matrix DeltaML(M_full,1);
			matrix Action(M_full,1);
			
			for(int i=0; i<M_full; i++){
				DeltaML.data[i]=0;
				if (Factor.data[i]>1e-12) {
					GoodFactor[i]=1;
				}
				Action.data[i]*=ACTION_REESTIMATE;
			}
			
			
			matrix UsedFactor(Used.size(),1);
			std::vector<int> index_c,index_c_neg;
			std::vector<int> iu,iu_neg;
			
			for (int i=0; i<UsedFactor.rows; i++) {
				UsedFactor.data[i]=Factor.data[Used[i]];
				GoodFactor[Used[i]]=0;
				if (Factor.data[Used[i]]>1e-12){
					index_c.push_back(Used[i]);
					iu.push_back(i);
				}
				else {
					index_c_neg.push_back(Used[i]);
					iu_neg.push_back(i);
				}
				
			}
			if (BasisAlignmentTest){
				for (int i=0; i<Aligned_out.size(); i++) {
					GoodFactor[Aligned_out[i]]=0;
				}
			}
			
			matrix NewAlpha(index_c.size(),1);
			matrix Delta(index_c.size(),1);
			
			
			for(int i=0; i<index_c.size(); i++){
				NewAlpha.data[i]=(S_out.data[index_c[i]]*S_out.data[index_c[i]])/Factor.data[index_c[i]];
				Delta.data[i]=(1.0/NewAlpha.data[i])-(1.0/Alpha.data[iu[i]]);
				DeltaML.data[index_c[i]]= (Delta.data[i]*(Q_in.data[index_c[i]]*Q_in.data[index_c[i]])/(Delta.data[i]*S_in.data[index_c[i]]+1)-log(1+S_in.data[index_c[i]]*Delta.data[i]))/2.0;
			}
			
			
			//FREE BASIS OPTION NOT AVAILABLE
			bool anytoDelete=0;
			if(index_c_neg.size()!=0 and M>1){
				for(int i=0; i<index_c_neg.size(); i++){
					DeltaML.data[index_c_neg[i]]= -(Q_out.data[index_c_neg[i]]*Q_out.data[index_c_neg[i]]/(S_out.data[index_c_neg[i]]+Alpha.data[iu_neg[i]])
													-log(1+S_out.data[index_c_neg[i]]/Alpha.data[iu_neg[i]]))/2.0;
					Action.data[index_c_neg[i]]=ACTION_DELETE;
					anytoDelete=1;
				}
			}
			
			
			//ANYTHING TO ADD
			index_c.clear();
			for(int i=0; i<GoodFactor.size(); i++){
				if(GoodFactor[i]==1){
					index_c.push_back(i);
				}
			}
			bool anytoADD=0;
			if(index_c.size()!=0){
				matrix quot(index_c.size(),1);
				for(int i=0; i<index_c.size(); i++){
					quot.data[i]=Q_in.data[index_c[i]]*Q_in.data[index_c[i]]/S_in.data[index_c[i]];
					DeltaML.data[index_c[i]]=(quot.data[i]-1-log(quot.data[i]))/2.0;
					Action.data[index_c[i]]=ACTION_ADD;
					anytoADD=1;
				}
			}
			
			
			//NOT TESTED? 
			if ((anytoADD && PriorityAddition) || (anytoDelete && PriorityDeletion)){
				//We won't perform re-estimation this iteration, which we achieve by
				//zero-ing out the delta
				for(int i=0; i<Action.rows; i++){
					if (Action.data[i]==ACTION_REESTIMATE)
						DeltaML.data[i]	= 0;
					//Furthermore, we should enforce ADD if preferred and DELETE is not
					// - and vice-versa
					if (anytoADD && PriorityAddition && PriorityDeletion){
						if (Action.data[i]==ACTION_DELETE)
							DeltaML.data[i]	= 0;
						
					}
					if (anytoDelete && PriorityDeletion && !PriorityAddition){
						if (Action.data[i]==ACTION_ADD)
							DeltaML.data[i]	= 0;	
						
					}
				}
			}
			
			//Choose the one with largest likelihood
			double deltaLogMallrginal=0.0;
			int nu=0;
			selectedAction=0;
			bool anyWorthwhileAction;
			for (int i=0; i<DeltaML.rows; i++) {
				if (DeltaML.data[i]>deltaLogMallrginal){
					//Updated 23/03/10 so that we do not delete when only one relevance vector
					if(Action.data[i]==-1 && Mu.rows==1){
						std::cout << "TRYING TO DELETE WHEN MU=1" << std::endl;
					}
					else{
						deltaLogMallrginal=DeltaML.data[i];
						nu=i;
						selectedAction=Action.data[i];
					}
				}
			}
			
			int j=0;
			anyWorthwhileAction=deltaLogMallrginal>0;
			if (selectedAction==ACTION_REESTIMATE || selectedAction==ACTION_DELETE){
				for	(int i=0; i<Used.size(); i++){
					if (Used[i]==nu)
						j=i;
				}
			}
			
			//DIFFERENCE TO MATLAB WITH NU being selected./RVM-Speed -b 0.99 -k Binary -d kbd420
			//MATLAB 72	0.7476216971706305	0.7476216971706305
			//C++ 
			/*
			 if (iteration_counter>220) {
			 printf("%d\t%.16f\t%.16f\n",nu,deltaLogMallrginal ,DeltaML.data[71]);
			 }
			 */
			
			matrix Phi;
			std::string act;
			Phi.AddColumn(BASIS, nu);
			
			double newAlpha=S_out.data[nu]*S_out.data[nu]/Factor.data[nu];
			
			if (!anyWorthwhileAction || (selectedAction==ACTION_REESTIMATE && abs(log(newAlpha)-log(Alpha.data[j]))<MinDeltaLogAlpha && !anytoDelete)){
				selectedAction=ACTION_TERMINATE;
				act="potential termination";
			}
			
			if (BasisAlignmentTest){
				if (selectedAction==ACTION_ADD){
					matrix p;
					mprod(Phi,PHI,p,1,1.0);
					if (p.rows>1){
						Rprintf("BELIEVED ERROR: check p, should be a single row\n");
						exit(1);
					}
					std::vector<int> findAligned;
					for (int i=0; i<p.cols; i++){
						if(p.data[i]>AlignmentMax){
							findAligned.push_back(i);
						}
					}
					int numAligned=findAligned.size();
					if (numAligned>0){
						selectedAction=ACTION_ALIGNMENT_SKIP;
						act="alignment-deferred addition";
						alignDeferCount+=1;
						for(int i=0; i<numAligned; i++){
							Aligned_out.push_back(nu);
							Aligned_in.push_back(Used[findAligned[i]]);
						}
					}
				}
				if (selectedAction==ACTION_DELETE){
					std::vector<int> findAligned;
					for(int i=0; i<Aligned_in.size(); i++){
						if(Aligned_in[i]==nu){
							findAligned.push_back(i);
						}
					}
					int numAligned=findAligned.size();
					//COuld include DIAGNOSTICS here
					if(numAligned>0){
						for (int i=(numAligned-1); i>-1; i--) {
							Aligned_in.erase(Aligned_in.begin()+findAligned[i]);
							Aligned_out.erase(Aligned_out.begin()+findAligned[i]);
						}
					}
				}
				
			}
			
			
			/*		
			 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			 
			 %% ACTION PHASE
			 %%
			 %% Implement above decision
			 %%
			 */
			bool UPDATE_REQUIRED=0;
			matrix SIGMANEW;
			switch (selectedAction) {
				case ACTION_REESTIMATE:{
					
					double oldAlpha=Alpha.data[j];
					Alpha.data[j]=newAlpha;
					matrix s_j;
					
					s_j.AddColumn(SIGMA, j);
					double deltaInv=1.0/(newAlpha-oldAlpha);
					double kappa=1.0/(SIGMA.data[j+j*SIGMA.rows]+deltaInv);
					
					matrix tmp(s_j.rows,s_j.cols);
					matrix deltaMu(tmp.rows,tmp.cols);
					
					for (int i=0; i<s_j.rows; i++) {
						for (int k=0; k<s_j.cols; k++) {
							tmp.data[i+s_j.rows*k]=s_j.data[i+s_j.rows*k]*kappa;
							deltaMu.data[i+s_j.rows*k]=tmp.data[i+s_j.rows*k]*-Mu.data[j];
							Mu.data[i+s_j.rows*k]+=deltaMu.data[i+s_j.rows*k];
						}
					}
					matrix SIGMANEW;
					mprod(tmp, s_j, SIGMANEW, 2, 1.0);
					
					for (int i=0; i<SIGMA.rows; i++) {
						for (int k=0; k<SIGMA.cols; k++) {
							SIGMANEW.data[i+SIGMANEW.rows*k]=SIGMA.data[i+SIGMA.rows*k]-SIGMANEW.data[i+SIGMANEW.rows*k];
						}
					}
					
					if (UpdateIteration){
						matrix tempbbpsj;
						mprod(BASIS_B_PHI,s_j,tempbbpsj,0,1.0);
						for (int i=0; i<tempbbpsj.rows; i++){
							for(int k=0; k<tempbbpsj.cols; k++){
								tempbbpsj.data[i+tempbbpsj.rows*k]=(tempbbpsj.data[i+tempbbpsj.rows*k]*tempbbpsj.data[i+tempbbpsj.rows*k]);
							}
						}
						for (int i=0; i<S_in.rows; i++) {
							S_in.data[i]=S_in.data[i]+kappa*(tempbbpsj.data[i]);
							
						}
						//printoutMatrix(S_in);
						matrix tmpq;
						mprod(BASIS_B_PHI, deltaMu, tmpq, 0, 1.0);
						for (int i=0; i<Q_in.rows; i++) {
							Q_in.data[i]=Q_in.data[i]-tmpq.data[i];
							
						}
					}
					updateCount+=1;
					act="re-estimation";
					UPDATE_REQUIRED=1;
				}
					break;
				case ACTION_ADD:{
					matrix B_Phi(beta.rows,beta.cols);
					for(int i=0; i<B_Phi.rows; i++){
						for (int k=0; k<B_Phi.cols; k++) {
							B_Phi.data[i+B_Phi.rows*k]=Phi.data[i+B_Phi.rows*k]*beta.data[i+B_Phi.rows*k];
						}
					}
					
					matrix BASIS_B_phi;
					
					mprod(BASIS,B_Phi,BASIS_B_phi,1,1.0);
					
					matrix tmp0;
					matrix tmp;
					mprod(B_Phi, PHI, tmp0, 1, 1.0);
					mprod(tmp0, SIGMA, tmp, 0, 1.0);
					tmp.rows=tmp.cols;
					tmp.cols=1;
					
					//cout << "tmp "<<endl;
					double *tmpalphadata=new double[Alpha.rows*Alpha.cols];
					
					
					Alpha.resize(Alpha.rows+1, Alpha.cols);
					Alpha.data[Alpha.rows-1]=newAlpha;
					PHI.AddColumn(Phi, 0);
					
					double s_ii=1/(newAlpha+S_in.data[nu]);
					matrix s_i(tmp.rows,1);
					for (int i=0; i<tmp.rows; i++) {
						s_i.data[i]=-s_ii*tmp.data[i];
					}
					matrix TAU;
					mprod(s_i,tmp,TAU,2,-1.0);
					SIGMANEW.resize(s_i.rows+1, s_i.rows+1);
					
					for(int i=0; i<SIGMANEW.rows; i++){
						for(int k=0; k<SIGMANEW.cols; k++){
							if(i<SIGMA.rows and k<SIGMA.cols){
								SIGMANEW.data[i+SIGMANEW.rows*k]=SIGMA.data[i+SIGMA.rows*k]+TAU.data[i+TAU.rows*k];
							}
							else if (i==SIGMA.rows and k<s_i.rows){
								SIGMANEW.data[i+SIGMANEW.rows*k]=s_i.data[k];
							}
							else if(k==SIGMA.rows and i<s_i.rows){
								SIGMANEW.data[i+SIGMANEW.rows*k]=s_i.data[i];
							}
							else if(i==SIGMA.rows and k==SIGMA.rows){
								SIGMANEW.data[i+SIGMANEW.rows*k]=s_ii;
							}
							else {
								Rprintf("ERROR IN ACTION ADD STATEMENT\n");
								exit(1);
							}
						}
					}
					double mu_i=s_ii*Q_in.data[nu];
					matrix deltaMu(tmp.rows+1,1);
					for (int i=0; i<deltaMu.rows-1; i++) {
						deltaMu.data[i]=-mu_i*tmp.data[i];
					}
					deltaMu.data[deltaMu.rows-1]=mu_i;
					//cout << "MU is being updated" <<endl;
					//printoutMatrix(deltaMu);
					Mu.resize(Mu.rows+1, Mu.cols);
					Mu.data[Mu.rows-1]=0;
					for (int i=0; i<Mu.rows; i++) {
						Mu.data[i]+=deltaMu.data[i];
						
					}
					//printoutMatrix(Mu);
					if(UpdateIteration){
						matrix mctmp;
						mprod(BASIS_B_PHI,tmp,mctmp,0,1.0);
						matrix mCi(mctmp.rows,mctmp.cols);
						matrix mCi2=mCi;
						
						for (int i=0; i<mCi.rows; i++) {
							for (int k=0; k<mCi.cols; k++) {
								mCi.data[i*mCi.rows+k]=BASIS_B_phi.data[i*mCi.rows+k]-mCi.data[i*mCi.rows+k];
								mCi2.data[i*mCi.rows+k]=mCi.data[i*mCi.rows+k]*mCi.data[i*mCi.rows+k];
							}
						}
						for (int i=0; i<S_in.rows; i++) {
							S_in.data[i]=S_in.data[i]-s_ii*mCi2.data[i];
							
						}
						//printoutMatrix(S_in);
						for (int i=0; i<Q_in.rows; i++) {
							Q_in.data[i]=Q_in.data[i]-mu_i*mCi.data[i];
							
						}
						
					}
					Used.push_back(nu);
					addCount+=1;
					act="addition";
					UPDATE_REQUIRED=true;
				}
					break;
				case ACTION_DELETE:{
					PHI.RemoveColumn(j);
					Alpha.RemoveRow(j);
					double s_jj=SIGMA.data[j+SIGMA.rows*j];
					matrix s_j;
					s_j.AddColumn(SIGMA, j);
					
					matrix tmp(s_j.rows,s_j.cols);
					for(int i=0; i<s_j.rows; i++){
						tmp.data[i]=s_j.data[i]/s_jj;
					}
					matrix SIGMANEW;
					
					mprod(tmp, s_j, SIGMANEW, 2, 1.0);
					
					for(int i=0; i<SIGMANEW.rows; i++){
						for (int k=0; k<SIGMANEW.cols; k++) {
							SIGMANEW.data[i+SIGMANEW.rows*k]=SIGMA.data[i+SIGMA.rows*k]-SIGMANEW.data[i+SIGMANEW.rows*k];
						}
					}
					SIGMANEW.RemoveRow(j);
					SIGMANEW.RemoveColumn(j);
					
					matrix deltaMu(tmp.rows,1);
					for (int i=0; i<tmp.rows; i++) {
						deltaMu.data[i]=-Mu.data[j]*tmp.data[i];
						Mu.data[i]+=deltaMu.data[i];
						
					}
					double mu_j=Mu.data[j];
					Mu.RemoveRow(j);
					if (UpdateIteration){
						matrix jPm;
						mprod(BASIS_B_PHI,s_j,jPm,0,1.0);
						matrix jPm2(jPm.rows,jPm.cols);
						for (int i=0; i<jPm.rows; i++) {
							for (int k=0; k<jPm.cols; k++) {
								jPm2.data[i+jPm.rows*k]=jPm.data[i+jPm.rows*k]*jPm.data[i+jPm.rows*k];
							}
						}
						for (int i=0; i<S_in.rows; i++) {
							S_in.data[i]=S_in.data[i]+jPm2.data[i]/s_jj;
							
						}
						//printoutMatrix(S_in);
						for (int i=0; i<Q_in.rows; i++) {
							Q_in.data[i]=Q_in.data[i]+jPm.data[i]*(mu_j/s_jj);
							
						}
						
					}
					Used.erase(Used.begin()+(j));
					deleteCount+=1;
					act="deletion";
					UPDATE_REQUIRED=true;
				}
					break;
					
				default:
					break;
			}
			M=Used.size();
			
			//cout << "ACTION: " << act << " of " << nu << " ("<<deltaLogMallrginal<<")"<<endl;
			
			
			/*
			 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			 
			 %% UPDATE STATISTICS
			 
			 % If we've performed a meaningful action,
			 % update the relevant variables
			 % 
			 */
			if (UPDATE_REQUIRED){
				if(UpdateIteration){
					Rprintf("UPDATE ITERATION NOT IMPLEMENTED \n");
					/*
					 S_out=S_in;
					 Q_out=Q_in;
					 matrix tmp(Used.size(),1);
					 for (int i=0; i<Used.size(); i++) {
					 tmp.data[i]=Alpha.data[i]/(Alpha.data[i]-S_in.data[Used[i]]);
					 S_out.data[Used[i]]=tmp.data[i]*S_in.data[Used[i]];
					 Q_out.data[Used[i]]=tmp.data[i]*Q_in.data[Used[i]];
					 }
					 for (int i=0; i<Q_out.size1(); i++) {
					 Factor.data[i]=(Q_out.data[i]*Q_out.data[i])-S_out.data[i];
					 }
					 SIGMA=SIGMANEW;
					 for (int i=0; i<Alpha.size1(); i++) {
					 Gamma.data[i]=1-Alpha.data[i]*SIGMA.data[i*SIGMA.cols+i];
					 }
					 matrix temp3(beta.size1(),M);		
					 for(int i=0; i<beta.size1(); i++){
					 for (int k=0; k<M; k++){
					 temp3.data[i*temp3.cols+k]=PHI.data[i*PHI+k]*beta.data[i];
					 }
					 }
					 BASIS_B_PHI=trans(matrixprod(temp3,BASIS,1));
					 */
				}
				else{
					double newLogML=0.0;
					fullstatistics(1, PHI, BASIS,BASIS2,beta,SIGMA,Mu,Alpha,newLogML,Targets,Used,Factor,S_out,Q_in,S_in,Q_out,BASIS_B_PHI,Gamma);
					deltaLogMallrginal=newLogML-logML;
				}
				if(UpdateIteration && deltaLogMallrginal<0){
					Rprintf("** Alert **  DECREASE IN LIKELIHOOD !!\n");
				}
				logML=logML+deltaLogMallrginal;
				count+=1;
				logMarginalLog.data[count]=logML;
			}
			/*
			 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			 
			 % END OF CYCLE PROCESSING
			 % 
			 % Check if termination still specified, and output diagnostics
			 % 
			 */
			double sumgamma=0.0;
			for (int i=0; i<Gamma.rows; i++) {
				sumgamma+=Gamma.data[i];
			}
			
			double sqrtbeta=0.0;
			for (int i=0; i<beta.rows; i++) {
				sqrtbeta+=sqrt(1.0/beta.data[i]);
			}		
			if(selectedAction==ACTION_TERMINATE){
				Rprintf("** Stopping at iteration %d (Max_delta_ml=%f) **",iteration_counter,deltaLogMallrginal);
				Rprintf("'%4d> L = %.6f\t Gamma = %.2f (M = %d)\n",iteration_counter, logML/N,sumgamma, M);
				break;
			}
			
			//NO TIME LIMIT AS OF YET!!!!!!
			
			
			if ((iteration_counter%monitor_its==0 || iteration_counter==1)){
				Rprintf("%5d> L = %.6f\t Gamma = %.2f (M = %d)\n",iteration_counter, logML/N, sumgamma, M);
			}
			
			if(iteration_counter==ItNum){
				break;
			}
			
		}
		/*
		 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		 
		 %%
		 %% POST-PROCESSING
		 %%
		 
		 %
		 % Warn if we exited the main loop without terminating automatically
		 % 
		 */
		if (selectedAction!=ACTION_TERMINATE){
			Rprintf("** Iteration limit: algorithm did not converge\n");
		}
		int total	= addCount + deleteCount + updateCount;
		if(BasisAlignmentTest)
			total	= total+alignDeferCount;
		
		Rprintf("Action Summary\n");
		Rprintf("==============\n");
		Rprintf("Added\t\t\t%6d (%.0f%%)\n",addCount, 100*addCount/(total*1.0));
		Rprintf("Deleted\t\t%6d (%.0f%%)\n",deleteCount, 100*deleteCount/(total*1.0));
		Rprintf("Reestimated\t%6d (%.0f%%)\n",updateCount, 100*updateCount/(total*1.0));
		if (BasisAlignmentTest && alignDeferCount){
			Rprintf("--------------\n");
			Rprintf("Deferred\t%6d (%.0f%%)\n",alignDeferCount, 100*alignDeferCount/(total*1.0));
		}
		Rprintf("==============\n");
		Rprintf("Total of %d likelihood updates\n", count);
		
		//std::vector<int> PARAMATERrev;
		//matrix PARAMATERval;
		
		//int *PARAMATERrev;
		//double *PARAMATERval;
		//PARAMATERrev = (int*)R_alloc( Used.size(),sizeof(int) );
		//PARAMATERval = (double*)R_alloc(Used.size(),sizeof(double));
		
		//std::vector<size_t> index;
		for (unsigned i = 0; i < Used.size(); ++i)
			PARAMATERrev[i]=Used[i];
		
		//	index.push_back(i);
		//sort(index.begin(), index.end(), index_cmp<std::vector<int>&>(PARAMATERrev));
		//sort (PARAMATERrev.begin(), PARAMATERrev.end());
		
		//PARAMATERval.reset(index.size(),1);
		for (int i=0; i<Used.size(); i++) {
			PARAMATERval[i]=(Mu.data[i]/(Scales.data[Used[i]]));
		}
		
		
	}


}

