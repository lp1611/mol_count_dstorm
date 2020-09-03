#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "matrix.h"
#define X_in prhs[0]
#define nu_in prhs[1]
#define BP0_in prhs[2]
#define BP1_in prhs[3]
#define log_lik plhs[0] 
double sum_log_vec(double*, int);
double sum_vec(double*, int); 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    mxArray *Mwork, *Mwork2,*Mwork4;
    mxArray *Mwork3;
    mxArray *prodwork;
    mxArray *cwork; 
    double sum=0.,*c,*lik,*X,*nu,*BP0,*BP1,*likwork,*alpha,*beta,*prods,log_likelihood = 0.;
	int n,emit;
    int i,j,k,m,x,nstates;
    
    emit = mxGetM(X_in);
    n = mxGetN(X_in);
    nstates = mxGetM(BP0_in); 
    X = mxGetPr(X_in);
	nu = mxGetPr(nu_in);
    BP0 = mxGetPr(BP0_in);
    BP1 = mxGetPr(BP1_in);
  
    Mwork = mxCreateDoubleMatrix(1, nstates, mxREAL);
    likwork = mxGetPr(Mwork);
    
    Mwork2 = mxCreateDoubleMatrix(n+1, nstates*emit, mxREAL);
    alpha = mxGetPr(Mwork2);
    
    Mwork4 = mxCreateDoubleMatrix(n+1, nstates*emit, mxREAL);
    beta = mxGetPr(Mwork4);
    
    Mwork3 = mxCreateDoubleMatrix(1, emit, mxREAL);
    lik = mxGetPr(Mwork3);
    
    prodwork = mxCreateDoubleMatrix(nstates, nstates, mxREAL);
    prods = mxGetPr(prodwork);
    
    cwork = mxCreateDoubleMatrix(1,n,mxREAL); 
    c = mxGetPr(cwork);
   
    
    for (j=0;j<emit;j++) {
            for (k=0;k<nstates;k++) {
                alpha[(nstates*j+k)*(n+1)] = nu[k];
            }
                x=X[j];
                for(k=0;k<nstates;k++) {
                    if (x==0) {
                        for (m=0;m<nstates;m++) {
                            likwork[k] += BP0[nstates*k+m]*nu[m];
                        }
                    } else {
                        for (m=0;m<nstates;m++) {
                            likwork[k] += BP1[nstates*k+m]*nu[m];
                        }
                    }
                        alpha[1+(nstates*j+k)*(n+1)] = likwork[k];
                }
                
                c[j] = sum_vec(likwork,nstates); 
                    for(k=0;k<nstates;k++) { 
                        likwork[k] = likwork[k]/c[j];
                        alpha[1+(nstates*j+k)*(n+1)] = likwork[k];
                    }
                
                
            for (i=1;i<n;i++) {
                x = X[j+i*emit];
                         for(m=0;m<nstates;m++) { 
                             if (x==0) {  
                                 for(k=0;k<nstates;k++) {
                                        alpha[i+1+(nstates*j+m)*(n+1)] += likwork[k]*BP0[k+nstates*m];
                                }
                             } else {
                                    for(k=0;k<nstates;k++) {
                                        alpha[i+1+(nstates*j+m)*(n+1)] += likwork[k]*BP1[k+nstates*m];
                                    }
                             }   
                         }

                    for(k=0;k<nstates;k++) { 
                        likwork[k] = alpha[i+1+(nstates*j+k)*(n+1)];
                    }
                c[i] = sum_vec(likwork,nstates); 
                       for(k=0;k<nstates;k++) { 
                            likwork[k] = likwork[k]/c[i];
                            alpha[i+1+(nstates*j+k)*(n+1)] = likwork[k]; 
                        }
            }
                log_likelihood += sum_log_vec(c,n);
      } 
    
    log_lik = mxCreateDoubleScalar(log_likelihood); 
    mxDestroyArray(cwork); 
    mxDestroyArray(Mwork);
    mxDestroyArray(Mwork2);
    mxDestroyArray(Mwork4);
    mxDestroyArray(prodwork);
    mxDestroyArray(Mwork3);
    return;
}

double sum_vec(double *x, int n) {
    int i; 
    double sum=0; 
    
    for (i=0;i<n;i++) { 
        sum+=x[i];
    }
    return(sum); 
}

double sum_log_vec(double *x, int n) {
    int i; 
    double sum=0.; 
    
        for (i=0;i<n;i++) { 
            sum+=log(x[i]);
        }
    return(sum); 
}

                
    
    

    
    