// Solves the inner problem in RatioDCA 
// for the constrained ncut problem with subset and volume constraints
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
    if (nrhs != 9 || nlhs != 2) {
        mexErrMsgTxt(
           "Usage: [Xbest,PrimalObjBest] = \n"
           "       mex_ip_cnstr_ncut_subset(W,c2,alpha,MAXITER,EPS,MaxSumSquaredWeights,c1,MAXITER_temp,debug)");
    }
  
    // get important parameters
    mwSize rows = mxGetM(prhs[0]);      // number of rows of W
    mwSize cols = mxGetN(prhs[0]);      // number of columns of W (should be the same)
    mwSize len  = mxGetM(prhs[1]);      // the desired output
    mwSize lenalpha = mxGetM(prhs[2]);  // alpha

    // check if input has right format
    if(!mxIsSparse(prhs[0])) { mexErrMsgTxt("Matrix is not sparse."); }
    if(rows!=cols) { mexErrMsgTxt("Sparse matrix is not square."); }
    if(rows!=len) { mexErrMsgTxt("Dimensions of W and c2 mismatch."); }
    
    // Create output array and compute values
    double* sr = mxGetPr(prhs[0]);     // get values
    mwIndex* irs = mxGetIr(prhs[0]);   // get row
    mwIndex* jcs = mxGetJc(prhs[0]);   // get columns
  
    double* c2 = mxGetPr(prhs[1]);           // c2
    double* alpha = mxGetPr(prhs[2]);  	     // alpha
    mwIndex MAXITER = mxGetScalar(prhs[3]); 
    double EPS = mxGetScalar(prhs[4]);       //for criterion fabs(dualval)<EPS
    double MaxSumSquaredWeights = mxGetScalar(prhs[5]); 
    double c1 = mxGetScalar(prhs[6]); 
    mwIndex MAXITER_temp = mxGetScalar(prhs[7]); 
    bool debug = (bool) mxGetScalar(prhs[8]);

    if(MaxSumSquaredWeights<=0) {
        mexErrMsgTxt("Lipschitz constant has to be positive");
    }

    // allocate memory for output
    plhs[0] = mxCreateDoubleMatrix(len,1,mxREAL);     // best output vector
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);       // best primal objective
  
    // pointers for output
    double* Xbest = mxGetPr(plhs[0]);
    double* PrimalObjBest = mxGetPr(plhs[1]);

    mwIndex counter = 0, i, j, iter = 0;
    double Dcur, betacur, factor;
    double dummy, normD, Fval;
    double* dummyPointer;
    double L = 1/MaxSumSquaredWeights;

    double* X = new double[len];
    double* D = new double[len];
    double* Dproj = new double[len];
    double* DprojBest = new double[len];
    double* beta = new double[lenalpha];
    double* betaold = new double[lenalpha];

    for(i=0; i<lenalpha; i++) { beta[i] = 0; }
    for(i=0; i<lenalpha; i++) { betaold[i] = 0; }
    for(i=0; i<len; i++) { DprojBest[i] = 0; }

    double temp, max_f, primalVal_orig;
    double primalVal = 0;
    double primalValBest = 0;

    // main loop
    Fval = EPS+1;
    while(iter<=MAXITER) {
        // exchange beta and betaold
        dummyPointer = beta; beta = betaold; betaold = dummyPointer;
 
        // compute X=Aalpha
        for(i=0; i<len; i++) { X[i] = 0; }
        counter = 0;
        for(j=0; j<cols; j++) {
            for(i=0; (unsigned int) i<jcs[j+1]-jcs[j]; i++) {
                dummy = sr[counter]*alpha[counter];
                X[j] -= dummy;
                X[irs[counter]] += dummy;
                counter++;
            }
        }
 
        // compute D = (-Aalpha-c2) / c1
        for(i=0; i<len; i++) {
            D[i] = (-X[i]-c2[i])/c1;
            if (D[i]<0) Dproj[i] = 0;
            else Dproj[i] = D[i];
        }
        Dproj = ProjKSimplex(Dproj, len, 1);

        // this projection is for the gradient.
        for (i=0; i<len; i++) {
            Dproj[i] = (D[i]-Dproj[i])*c1;
            if (Dproj[i]<0) Dproj[i] = 0;
        }
   
        // compute primal and dual objective and check if converged
        if (iter==0 || iter==MAXITER_temp || iter==MAXITER) {
            // compute dual objective
            normD = 0;
            for(i=0; i<len; i++) { normD += Dproj[i]*Dproj[i]; }
            Fval = -normD;

            // compute original inner objective
            primalVal_orig = 0;
            counter = 0;
            for(j=0; j<cols; j++) {
                Dcur = Dproj[j];
                for(i=0; (unsigned int) i<jcs[j+1]-jcs[j]; i++) {
                    primalVal_orig += fabs(sr[counter]*( Dproj[irs[counter]] - Dcur));
                    counter++;
                }
            }

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
            for(i=0; (unsigned int) i<jcs[j+1]-jcs[j]; i++) {
                temp = sr[counter]*( Dproj[irs[counter]] - Dcur);
                betacur = alpha[counter] + L*temp;
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

