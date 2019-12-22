// Solves the inner problem in RatioDCA 
// for the constrained ncut problem with volume constraints
//
// (C)2012-19 Thomas Buehler, Syama Rangapuram, Simon Setzer and Matthias Hein

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
           "       mex_ip_cnstr_ncut(W,c2,alpha,MAXITER,EPS,MaxSumSquaredWeights,MAXITER_temp,debug)");
    }
  
    // get important parameters
    mwSize rows = mxGetM(prhs[0]);      // number of rows of W
    mwSize cols = mxGetN(prhs[0]);      // number of columns of W (should be the same)
    mwSize len  = mxGetM(prhs[1]);      // the desired output
    mwSize lenalpha = mxGetM(prhs[2]);  // alpha
  
    // check if input has right format
    if(!mxIsSparse(prhs[0])) { mexErrMsgTxt("Matrix is not sparse"); }
    if(rows!=cols) { mexErrMsgTxt("Sparse matrix is not square"); }
    if(rows!=len) { mexErrMsgTxt("Dimensions of W and c2 mismatch."); }
    
    // Create output array and compute values
    double* sr = mxGetPr(prhs[0]);     // get values
    mwIndex* irs = mxGetIr(prhs[0]);   // get row
    mwIndex* jcs = mxGetJc(prhs[0]);   // get columns
  
    double* c2 = mxGetPr(prhs[1]);		     
    double* alpha = mxGetPr(prhs[2]);
    mwIndex MAXITER =  mxGetScalar(prhs[3]); 
    double EPS = mxGetScalar(prhs[4]);       //for criterion fabs(dualval)<EPS
    double MaxSumSquaredWeights = mxGetScalar(prhs[5]); 
    mwIndex MAXITER_temp = mxGetScalar(prhs[6]); 
    bool debug = (bool) mxGetScalar(prhs[7]);

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
    double* Dproj2 = new double[len];
    double* DprojBest = new double[len];
    double* beta = new double[lenalpha];
    double* betaold = new double[lenalpha];

    for(i=0; i<lenalpha; i++) { beta[i] = 0; }
    for(i=0; i<lenalpha; i++) { betaold[i] = 0; }
    for(i=0; i<len; i++) DprojBest[i] = 0;

    double temp, primalVal2, primalVal_orig;
    double primalVal = 0;
    double primalValBest = 0;
  
    // double C = 0; 
    // for(i=0;i<len;i++) { C += c2[i]*c2[i]; }
    // C = 0.5*C;

    clock_t start, ende;
    float time1 = 0, time2 = 0, time3 = 0, time4 = 0;

    // main loop
    Fval = EPS+1;
    while(iter<MAXITER) {
        iter++;
        start=clock();
      
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
        
        // compute D = -Aalpha-c2
        for(i=0; i<len; i++) { 
            D[i] = (-X[i]-c2[i]);
        
            //Dproj[i]=(Y[i]-lambda*X[i])/c1; 
            //if (Dproj[i]<0) {Dproj[i]=0;}
            // Dproj[i] = D[i];
            if (D[i]<0) Dproj[i] = 0;
            else Dproj[i] = D[i];
        }
    
        ende = clock();
        time1 += (float) (ende-start);
        start=clock();

        // update beta and alpha 
        counter = 0;
        factor = (double) iter/(iter+3); 
        for(j=0; j<cols; j++) {   
            Dcur = Dproj[j];
            for(i=0;  i<jcs[j+1]-jcs[j]; i++) {             
                // update of beta
                temp = sr[counter]*( Dproj[irs[counter]] - Dcur);
                betacur = alpha[counter] + L*temp;
                // projection onto l_inf-cube 
                if(betacur>1) betacur = 1;
	            else if(betacur<-1) betacur = -1;
                beta[counter] = betacur;
                // update of alpha
                alpha[counter] = betacur + factor*(betacur-betaold[counter]);
                counter++;
            }
        }

        ende = clock();
        time2 += (float) (ende-start);
    
        // compute primal and dual objective and check if converged
        if (iter==1 || iter==MAXITER_temp || iter==MAXITER) {
            start = clock();
       
            // compute X=Aalpha
            for(i=0; i<len; i++) { X[i] = 0; }
            counter = 0;
            for(j=0; j<cols; j++) {   
                for(i=0;  i<jcs[j+1]-jcs[j]; i++) {  
                    dummy = sr[counter]*beta[counter];
                    X[j] -= dummy;
                    X[irs[counter]] += dummy;
                    counter++;
                }
            }
          
            // compute D = -Aalpha-c2
            for(i=0; i<len; i++) { 
                D[i] = (-X[i]-c2[i]);///eta;
                // Dproj2[i] = D[i];
                if (D[i]<0) Dproj2[i] = 0;
                else Dproj2[i] = D[i];
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
        
            // <c_2, f> 
            for (i=0; i<len; i++) {
                primalVal_orig += Dproj2[i]*c2[i];
            }

            primalVal = primalVal_orig+ 0.5*normD;
        
            if (normD >0) {
                primalVal_orig = primalVal_orig / sqrt(normD);
            }
            else {
                primalVal_orig = 0;
            }

            ende = clock();
            time2 += (float) (ende-start);
      
        
            // check if we are better than previous best one 
            if (primalVal_orig<primalValBest) {
                primalValBest = primalVal_orig;
                for(i=0; i<len; i++) { DprojBest[i] = Dproj2[i]; }
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
    }
  
    if (debug) {
        mexPrintf("...... time1: %.2f   time2: %.2f  time3: %.2f  time4: %.2f\n",time1 / (float)CLOCKS_PER_SEC, time2 / (float)CLOCKS_PER_SEC,time3 / (float)CLOCKS_PER_SEC, time4 / (float)CLOCKS_PER_SEC);
    }

    for(i=0; i<len; i++) { Xbest[i] = DprojBest[i]; }
    PrimalObjBest[0] = primalValBest;
    delete[] X; delete[] betaold; delete[] beta; delete[] Dproj; delete[] Dproj2; delete[] DprojBest; delete[] D;
}



// projection on the k simplex
double* ProjKSimplex(double* f, int len, int K) {
    int i;
    //double* xnew = new double[len];
    int* IX = new int[len];
    double sumxnew = 0;
    double sumxnew2;
    int counter = len;
    int counter2 = 0;
    for(i=0; i<len; i++) { IX[i] = i; sumxnew += f[i]; }
    int iteration = 0;
    while(iteration<10000) {
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
    if(iteration==10000)
        mexPrintf("Bug in Projection: %i, Fval: %i\n", counter, counter2);
    delete[] IX;
    return f;
    /*x = f;  IX=1:num;  
    while(cond)
      xnew = x;
      xnew(IX) = x(IX) - (sum(x)-k)/length(IX);
    if(sum(xnew>=0)==num)
      Projf=xnew; break;
    else
      %ixx = find(xnew<0);
      %IX = setdiff(IX,ixx);
      IX = find(xnew>0);
      x = max(xnew,0);
    end
    end*/
}

