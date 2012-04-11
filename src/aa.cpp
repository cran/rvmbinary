

extern "C" {
	#include <math.h>
	
	void aakernel(double *a, double *b,double *c,int *arows, int *acols) {
		
		int rows=arows[0];
		int cols=acols[0];
		double basisWidth=b[0];
		
		for (int i=0; i<rows; i++) {
			for (int j=i; j<rows; j++) {
				double sum=0.0;
				for (int k=0; k<cols; k++) {
					sum+=fabs(a[i+k*rows]-a[j+k*rows]);
				}
				
				c[i+j*rows]=pow(basisWidth,(cols-sum))*pow(1-basisWidth,sum);
				if (j!=i)
					c[j+rows*i]=pow(basisWidth,(cols-sum))*pow(1-basisWidth,sum);
			}
		}
	}
	
	
	void aapredict(double *a,double *t, double *b,double *c,int *arows, int *acols, int *trows) {
		
		int rows=arows[0];
		int cols=acols[0];
		int rowstest=trows[0];
		double basisWidth=b[0];
		
		for (int i=0; i<rows; i++) {
			for (int j=0; j<rowstest; j++) {
				double sum=0.0;
				for (int k=0; k<cols; k++) {
					sum+=fabs(a[i+k*rows]-t[j+k*rowstest]);
				}
				
				c[i+j*rows]=pow(basisWidth,(cols-sum))*pow(1-basisWidth,sum);
			}
		}
	}
}
