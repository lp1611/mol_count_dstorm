#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "matrix.h"
#define nu_in prhs[0]
#define B0_in prhs[1]
#define B1_in prhs[2]
#define n_F_in prhs[3]
//#define Mwork prhs[4] //this is an N_F by n_states*(N_F+1) zeros matrix supplied by MATLAB
#define Mwork_probs plhs[0] 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    mxArray *Mwork,*Mwork2;

    double sum=0,*nu,*B0, *B1,*M,*M2,*probs;
	int n,n_states;
    int i,j,k,m,x;
    
    n = (int) mxGetScalar(n_F_in); /*number of frames*/
    nu = mxGetPr(nu_in);
    B0 = mxGetPr(B0_in);
    n_states = mxGetN(B0_in); 
    B1 = mxGetPr(B1_in);
    //M = mxGetPr(Mwork);
    Mwork = mxCreateDoubleMatrix(n+1, n_states, mxREAL); /*M matrix BEFORE*/
    M = mxGetPr(Mwork);
    Mwork2 = mxCreateDoubleMatrix(n+1, n_states, mxREAL); /*M matrix AFTER*/
    M2 = mxGetPr(Mwork2);
    
    Mwork_probs = mxCreateDoubleMatrix(1, n+1, mxREAL); /*end probabilities*/
    probs = mxGetPr(Mwork_probs);
    
    //This calculates the base case
	for (i=0;i<n_states;i++) {
        for(m=0;m<n_states;m++) { 
            M[i*(n+1)] += nu[m]*B0[m+n_states*i]; //this is correct 
            M[i*(n+1) + 1] += nu[m]*B1[m+n_states*i]; //this is correct 
        }
	}
    
    if (n > 1) { 
    
        for(j=1;j<n;j++) { //Works out M[:,j]

            for (i=0;i<n_states;i++) {
                for (m=0;m<n_states;m++) {
                    M2[i*(n+1)] += M[m*(n+1)]*B0[m+n_states*i]; 
                }
            }

            for (k=1;k<=j;k++) {
                for (i=0;i<n_states;i++) { 
                    for (m=0;m<n_states;m++) {
                        M2[i*(n+1)+k] += M[m*(n+1)+k]*B0[m+n_states*i] + M[m*(n+1)+(k-1)]*B1[m+n_states*i];
                    }
                }
            }

            for (i=0;i<n_states;i++) {
                for (m=0;m<n_states;m++) {
                    M2[i*(n+1)+(j+1)] += M[m*(n+1)+j]*B1[m+n_states*i];
                }
            }

            for (m=0; m<n_states; m++) { 
               for (i=0;i<j+2;i++)  {                
                    M[m*(n+1) +i] = M2[m*(n+1) + i]; 
                    M2[m*(n+1) +i] = 0.; 
                }
            }

        }
    
    }
    
    mxDestroyArray(Mwork2);
    
    for(j=0;j<=n;j++) { 
        for (m=0;m<n_states;m++) {
            probs[j] += M[m*(n+1)+j];
        }
    }

    mxDestroyArray(Mwork);
    return;
}

