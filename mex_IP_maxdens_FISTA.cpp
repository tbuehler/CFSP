// Solves the inner problem in RatioDCA 
// for the general constrained maximum density subgraph problem
//

#include <math.h>
#include "mex.h"
#include "matrix.h"
#include <time.h>
#include "float.h"

double* ProjKSimplex(double*, int, int);


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    // Test number of parameters.
    if (nrhs != 8 || nlhs != 2) {
        mexErrMsgTxt(
           "Usage: [Xbest,PrimalObjBest] = \n"
           "       mex_IP_maxdens_FISTA(W,c2,alpha,MAXITER,L,lambda,c1,debug)");
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
    double* Dproj2 = new double[len];
    double* DprojBest = new double[len];
    double* beta = new double[lenalpha];
    double* betaold = new double[lenalpha];

    for(i=0; i<lenalpha; i++) { beta[i] = 0; }
    for(i=0; i<lenalpha; i++) { betaold[i] = 0; }
    for(i=0; i<len; i++) { DprojBest[i] = 0; }

    double temp, primalVal2;
    double primalVal = 0, primalValBest = 0, primalVal_orig = 0;
    double max_f;

    clock_t start, ende;
    float time1 = 0, time2 = 0, time3 = 0, time4 = 0;
    
    unsigned int MAXITER_temp=50;
    Fval=-DBL_MAX;
   
    //main loop
    while(iter< (unsigned int) MAXITER) {
        iter++;

        // exchange beta and betaold
        dummyPointer = beta; beta = betaold; betaold = dummyPointer;
       
        // compute X=Aalpha
        for(i=0; i<len; i++) { X[i]=0; }
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
            D[i] = (-c2[i]-lambda*0.5*X[i])/c1;
            if (D[i]<0) Dproj[i]=0;
            else Dproj[i]=D[i];
        }

        // project onto simplex
        Dproj = ProjKSimplex(Dproj, len, 1);
        
        for (i=0;i<len;i++) {
            Dproj[i] = (D[i]-Dproj[i])*c1;
            if (Dproj[i]<0) Dproj[i] = 0;
        }
        // z computed

        // update beta and alpha
        counter = 0;
        //tnew = (1 + sqrt(1+4*told*told))/2;
        //factor = (told-1)/tnew;
        factor = (double) iter/(iter+3);
        for(j=0; j<cols; j++) {
            Dcur = Dproj[j];
            for(i=0; i<jcs[j+1]-jcs[j]; i++) {
                temp = sr[counter]*( Dproj[irs[counter]] - Dcur);
                // update of beta
                betacur = alpha[counter] + lamByL*temp;
                // projection onto l_inf-cube
                if(betacur>1) betacur = 1;
                else if(betacur<-1) betacur = -1;
                beta[counter] = betacur;
                // update of alpha
                alpha[counter] = betacur + factor*(betacur-betaold[counter]);
                counter++;
            }
        }
        
        // compute primal and dual objective
        if (iter==1 || iter==MAXITER_temp || iter== (unsigned int) MAXITER) {
            start = clock();
       
            // compute X=Aalpha
            for(i=0; i<len; i++) { X[i] = 0; }
            counter = 0;
            for(j=0; j<cols; j++) {
                for(i=0; i<jcs[j+1]-jcs[j]; i++) {
                    dummy = sr[counter]*beta[counter];
                    X[j] -= dummy;
                    X[irs[counter]] += dummy;
                    counter++;
                }
            }
            
            // compute D = (-c2- lambda Aalpha) / c1
            for(i=0; i<len; i++) {
                D[i] = (-c2[i]-lambda*0.5*X[i])/c1;
                if (D[i]<0) Dproj2[i]=0;
                else Dproj2[i]=D[i];
            }
            
            // project on simplex
            Dproj2 = ProjKSimplex(Dproj2, len, 1);
            
            for (i=0;i<len;i++) {
                Dproj2[i] = (D[i]-Dproj2[i])*c1;
                if (Dproj2[i]<0) Dproj2[i] = 0;
            }
    
            // compute dual objective
            normD = 0;
            for(i=0; i<len; i++) { normD += Dproj2[i]*Dproj2[i]; }
            Fval = -normD;
                        
            // compute primal objective (original)
            counter = 0;
            primalVal_orig = 0;
            // 0.5 sum w_ij |f_i - f_j|
            for(j=0; j<cols; j++) {
                Dcur = Dproj2[j];
                for(i=0; i<jcs[j+1]-jcs[j]; i++) {
                    primalVal_orig += fabs(sr[counter]*( Dproj2[irs[counter]] - Dcur));
                    counter++;
                }
            }
            primalVal_orig = lambda * 0.5 * primalVal_orig;

            // c1 max_f
            max_f = 0;
            for (i=0; i<len; i++) {
                if (Dproj2[i]>max_f)
                    max_f = Dproj2[i];
            }
            primalVal_orig += c1*max_f;

            // <c2,f>
            for (i=0;i<len; i++) {
                primalVal_orig += Dproj2[i]*c2[i];
            }
            primalVal = primalVal_orig + 0.5*normD;


            if (normD >0) {
                primalVal_orig = primalVal_orig/ sqrt(normD);
            } else {
                primalVal_orig = 0;
            }
            
            // check if we are better than previous best one
            if (primalVal_orig<primalValBest) {
                primalValBest = primalVal_orig;
                for(i=0; i<len; i++) {
                    DprojBest[i] = Dproj2[i];
                }
            }
                        
            ende = clock();
            time4+=(float) (ende-start);
            
            // output values
            if (debug) {
                mexPrintf("...... it=%i\tinnerobj=%1.6f\tprimalobj=%1.6f\tdualobj=%1.6f\tgap=%1.6f\n",
                          iter, primalVal_orig, primalVal, Fval, primalVal-Fval);
            }

            // check if converged
            if (primalVal<0)
                break;
            else if(iter==MAXITER_temp) 
                MAXITER_temp = MAXITER_temp*2; 
        }
    }
    
    // output final results
    if (debug) {
        mexPrintf("...... time1: %.2f   time2: %.2f  time3: %.2f   time4: %.2f\n",
                  time1 / (float)CLOCKS_PER_SEC, time2 / (float)CLOCKS_PER_SEC, time3 / (float)CLOCKS_PER_SEC, time4 / (float)CLOCKS_PER_SEC);
        mexPrintf("...... FINAL it=%i\tinnerobj=%1.6f\tprimalobj=%1.6f\tdualobj=%1.6f\tgap=%1.6f\n",
                  iter, primalVal_orig, primalVal, Fval, primalVal-Fval);
    }
    
    // write output
    for(i=0; i<len; i++) { Xbest[i] = DprojBest[i]; }
    PrimalObjBest[0] = primalValBest;
    delete[] X; delete[] betaold; delete[] beta; delete[] Dproj; delete[] Dproj2; delete[] DprojBest; delete[] D;
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

