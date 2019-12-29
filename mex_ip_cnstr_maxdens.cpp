// Solves the inner problem in RatioDCA 
// for the general constrained maximum density subgraph problem
//
// (C)2012-19 Thomas Buehler, Syama Rangapuram, Simon Setzer and Matthias Hein

#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "float.h"

double* ProjKSimplex(double*, int, int);


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    // Test number of parameters.
    if (nrhs != 8 || nlhs != 2) {
        mexErrMsgTxt(
           "Usage: [Xbest,PrimalObjBest] = \n"
           "       mex_ip_cnstr_maxdens(W,c2,alpha,MAXITER,L,lambda,c1,debug)");
    }
  
    // get important parameters
    unsigned int rows = (unsigned int)mxGetM(prhs[0]);     // number of rows
    unsigned int cols = (unsigned int)mxGetN(prhs[0]);     // number of columns
    unsigned int len  = (unsigned int)mxGetM(prhs[1]);     // dimension
    unsigned int lenalpha = (unsigned int)mxGetM(prhs[2]); // alpha

    // check if input has right format
    if(!mxIsSparse(prhs[0])) { mexErrMsgTxt("Matrix is not sparse."); }
    if(rows!=cols) { mexErrMsgTxt("Sparse matrix is not square."); }
    if(rows!=len) { mexErrMsgTxt("Dimensions of W and c2 mismatch."); }
    
    // Create output array and compute values
    double* sr = mxGetPr(prhs[0]);     // get values
    mwIndex* irs = mxGetIr(prhs[0]);   // get row
    mwIndex* jcs = mxGetJc(prhs[0]);   // get columns
  
    // read additional input params
    double* c2 = mxGetPr(prhs[1]);
    double* alpha = mxGetPr(prhs[2]);             // main variable
    int MAXITER = (int) mxGetScalar(prhs[3]);     // max. number of iterations
    double L = mxGetScalar(prhs[4]);              // Lipschitz constant
    double lambda = (double)mxGetScalar(prhs[5]); // lambda
    double c1 = (double) mxGetScalar(prhs[6]);
    bool debug = (bool) mxGetScalar(prhs[7]);

    if(MAXITER<=0) { mexErrMsgTxt("MAXITER has to be positive."); }
    if(L<=0) { mexErrMsgTxt("Lipschitz constant has to be positive"); }
        
    // allocate memory for output
    plhs[0] = mxCreateDoubleMatrix(len,1,mxREAL);     // best output vector
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);       // best primal objective
  
    // create pointers to output
    double* Xbest = mxGetPr(plhs[0]);
    double* PrimalObjBest = mxGetPr(plhs[1]);

    // some helpers
    int counter = 0;
    unsigned int i, j, iter = 0;
    double Dcur, betacur, factor;
    double dummy, normD, Fval;
    double* dummyPointer;
    double lamByL = lambda/L;

    double* X = new double[len];
    double* D = new double[len];
    double* Dproj = new double[len];
    double* DprojBest = new double[len];
    double* beta = new double[lenalpha];
    double* betaold = new double[lenalpha];

    for(i=0; i<lenalpha; i++) { beta[i] = 0; }
    for(i=0; i<lenalpha; i++) { betaold[i] = 0; }
    for(i=0; i<len; i++) { DprojBest[i] = 0; }

    double temp, max_f;
    double primalVal = 0, primalValBest = 0, primalVal_orig = 0;

    unsigned int MAXITER_temp=50;
    Fval=-DBL_MAX;
   
    //main loop
    while(iter<= (unsigned int) MAXITER) {
        // exchange beta and betaold
        dummyPointer = beta; beta = betaold; betaold = dummyPointer;
 
        // compute X=Aalpha
        for(i=0; i<len; i++) { X[i] = 0; }
        counter = 0;
        for(j=0; j<cols; j++) {
            for(i=0; i<jcs[j+1]-jcs[j]; i++) {
                dummy = sr[counter]*alpha[counter];
                X[j] -= dummy;
                X[irs[counter]] += dummy;
                counter++;
            }
        }
 
        // compute D = (-c2- lambda Aalpha) / c1
        for(i=0; i<len; i++) {
            D[i] = (-c2[i]-lambda*X[i])/c1;
            if (D[i]<0) Dproj[i] = 0;
            else Dproj[i] = D[i];
        }
        Dproj = ProjKSimplex(Dproj, len, 1);

        for (i=0;i<len;i++) {
            Dproj[i] = (D[i]-Dproj[i])*c1;
            if (Dproj[i]<0) Dproj[i] = 0;
        }

        // compute primal and dual objective and check if converged
        if (iter==0 || iter==MAXITER_temp || iter== (unsigned int) MAXITER) {
            // compute dual objective
            normD = 0;
            for(i=0; i<len; i++) { normD += Dproj[i]*Dproj[i]; }
            Fval = -normD;

            // compute original inner objective
            primalVal_orig = 0;            
            counter = 0;
            for(j=0; j<cols; j++) {
                Dcur = Dproj[j];
                for(i=0; i<jcs[j+1]-jcs[j]; i++) {
                    primalVal_orig += fabs(sr[counter]*( Dproj[irs[counter]] - Dcur));
                    counter++;
                }
            }
            primalVal_orig = lambda * primalVal_orig;

            max_f = 0;
            for (i=0; i<len; i++) {
                if (Dproj[i]>max_f)
                    max_f = Dproj[i];
            }
            primalVal_orig += c1*max_f;

            for (i=0; i<len; i++) {
                primalVal_orig += Dproj[i]*c2[i];
            }

            // compute modified primal objective
            primalVal = primalVal_orig + 0.5*normD;

            if (normD >0) {
                primalVal_orig = primalVal_orig / sqrt(normD);
            } else {
                primalVal_orig = 0;
            }

            // check if we are better than previous best one
            if (primalVal_orig<primalValBest) {
                primalValBest = primalVal_orig;
                for(i=0; i<len; i++) { DprojBest[i] = Dproj[i]; }
            }

            if (debug) {
                mexPrintf("...... it=%i\tinnerobj=%1.6f\tprimalobj=%1.6f\tdualobj=%1.6f\tgap=%1.6f\n",
                          iter, primalVal_orig, primalVal, Fval, primalVal-Fval);
            }

            if (primalVal<0)
                break;
            else if(iter==MAXITER_temp) 
                MAXITER_temp = MAXITER_temp*2;
        }

        // update beta and alpha
        counter = 0;
        factor = (double) iter/(iter+3);
        for(j=0; j<cols; j++) {
            Dcur = Dproj[j];
            for(i=0; i<jcs[j+1]-jcs[j]; i++) {
                temp = sr[counter]*( Dproj[irs[counter]] - Dcur);
                betacur = alpha[counter] + lamByL*temp;
                if(betacur>1) betacur = 1;
                else if(betacur<-1) betacur = -1;
                beta[counter] = betacur;
                alpha[counter] = betacur + factor*(betacur-betaold[counter]);
                counter++;
            }
        }
        iter++;
    }

    for(i=0; i<len; i++) { Xbest[i] = DprojBest[i]; }
    PrimalObjBest[0] = primalValBest;
    delete[] X; delete[] betaold; delete[] beta; delete[] Dproj; delete[] DprojBest; delete[] D;
}



// projection on the k simplex
double* ProjKSimplex(double* f, int len, int K) {
    int MAXIT = 10000;
    int i;
    int* IX = new int[len];
    double sumxnew = 0;
    double sumxnew2;
    int counter = len;
    int counter2 = 0;
    for(i=0; i<len; i++) { IX[i] = i; sumxnew += f[i]; }
    int iteration = 0;

    while(iteration<MAXIT) {
        counter2 = 0;
        sumxnew2 = 0;
        for(i=0; i<counter; i++) {
            f[IX[i]] = f[IX[i]] - (sumxnew-K)/counter;
            if(f[IX[i]]>0) {
                IX[counter2] = IX[i];
                counter2++;
                sumxnew2 += f[IX[i]];
            } else f[IX[i]] = 0;
        }
        if(counter2==counter) break;
        counter = counter2;
        sumxnew = sumxnew2;
        iteration++;
    }
    if(iteration==MAXIT)
        mexPrintf("Bug in Projection: %i, Fval: %i\n", counter, counter2);
    delete[] IX;
    return f;
}

