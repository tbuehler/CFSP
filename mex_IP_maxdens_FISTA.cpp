// Solves the inner problem in RatioDCA 
// for the general constrained maximum density subgraph problem
//

#include <math.h>
#include "mex.h"
#include "matrix.h"
#include <time.h>
#include "float.h"

double* ProjKSimplex(double*, int,int);


void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    // Test number of parameters.
    if (nrhs != 8 || nlhs != 8) {
        mexWarnMsgTxt("Usage: [X,alpha,PrimalObj,DualObj,iter,Xbest,PrimalObjBest,DualObjBest]=mex_IP_maxdens_FISTA(W,Y,alpha,MAXITER,L,lambda,eta,debug)");
        return;
    }
    
    // get important parameters
    unsigned int rows = (unsigned int)mxGetM(prhs[0]); // number of rows
    unsigned int cols = (unsigned int)mxGetN(prhs[0]); // number of columns
    unsigned int len  = (unsigned int)mxGetM(prhs[1]); // dimension
    unsigned int lenalpha = (unsigned int)mxGetM(prhs[2]); // alpha
    
    // check if input as right format
    if(!mxIsSparse(prhs[0])) {
        mexWarnMsgTxt("Matrix is not sparse.");
        return;
    }
    
    if(rows!=cols){
        mexWarnMsgTxt("Sparse matrix is not square.");
        return;
    }
    
    if(rows!=len){
        mexWarnMsgTxt("Dimensions of W and Y mismatch.");
        return;
    }
    
    // Create output array and compute values
    double* sr = mxGetPr(prhs[0]);     // get values
    mwIndex* irs = mxGetIr(prhs[0]);   // get row
    mwIndex* jcs = mxGetJc(prhs[0]);   // get columns
    
    double* Y = mxGetPr(prhs[1]);
    double* alpha = mxGetPr(prhs[2]);  // main variable
    
    // get maximum number of iterations
    int MAXITER = (int) mxGetScalar(prhs[3]);
    if(MAXITER<=0){
        mexWarnMsgTxt("MAXITER has to be positive.");
        return;
    }
    
    // get Lipschitz constant
    double L = mxGetScalar(prhs[4]);
    if(L<=0){
        mexWarnMsgTxt("Lipschitz constant has to be positive");
        return;
    }
    
    double lambda = (double)mxGetScalar(prhs[5]); // lambda
    double eta = (double) mxGetScalar(prhs[6]);
    bool debug = (bool) mxGetScalar(prhs[7]);
    
    
    // allocate memory for output
    plhs[0] = mxCreateDoubleMatrix(len,1,mxREAL);     // output vector
    plhs[1] = mxCreateDoubleMatrix(lenalpha,1,mxREAL);// dual variable
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);       // primal objective
    plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);       // dual objective
    plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL);       // final iteration
    plhs[5] = mxCreateDoubleMatrix(len,1,mxREAL);     // best output vector
    plhs[6] = mxCreateDoubleMatrix(1,1,mxREAL);       // best objective
    plhs[7] = mxCreateDoubleMatrix(1,1,mxREAL);       // best dual objective
    
    // create pointers to output
    double* X = mxGetPr(plhs[0]);
    double* alphaout = mxGetPr(plhs[1]);
    double* PrimalObj = mxGetPr(plhs[2]);
    double* DualObj = mxGetPr(plhs[3]);
    double* FinalIter = mxGetPr(plhs[4]);
    double* Xbest = mxGetPr(plhs[5]);
    double* PrimalObjBest = mxGetPr(plhs[6]);
    double* DualObjBest = mxGetPr(plhs[7]);
    
    // some helpers
    int counter=0;
    unsigned int i,j,iter=0;
    double Dcur,betacur,factor;
    double dummy,normD,Fval;
    double* dummyPointer;
    double lamByL=lambda/L;
    
    double* D =new double[len];
    double* Dproj =new double[len];
    double* Dproj2 =new double[len];
    double* DprojBest =new double[len];
    double* beta    = new double[lenalpha];
    double* betaold = new double[lenalpha];
    for(i=0; i<lenalpha; i++) { beta[i]=0; }
    for(i=0; i<lenalpha; i++) { betaold[i]=0; }
    double temp,primalVal2,dualValBest;
    double primalVal=0,primalValBest=0,primalVal_orig=0;double max_f;
    bool computeObjective;
    
    // step sizes in fista
   // double tnew=1, told=1;
    
 //   double C=0;
 //   for(i=0;i<len;i++){C+=Y[i]*Y[i];}
 //   C=0.5*C;
      
    for(i=0; i<len; i++) DprojBest[i]=0;
    
    clock_t start,ende;
    float time1=0,time2=0, time3=0,time4=0;
    
    unsigned int MAXITER_temp=50;
    
    Fval=-DBL_MAX;
    dualValBest=Fval;
    
   
    //main loop
    while(iter< (unsigned int) MAXITER)
    {
        // defines strategy when to check if converged
        if (iter== (unsigned int) MAXITER-1 || iter==MAXITER_temp)     {
            computeObjective=true;
        } else {
            computeObjective =false;
        }
        //if (iter % 10 ==0){
        //    computeObjective=true;
        //}
        
        // exchange beta and betaold
        dummyPointer=beta; beta=betaold; betaold=dummyPointer;
        
        // exchange tnew and told
        //told=tnew;
    
        
        // initialize X with zeros
        for(i=0; i<len; i++) { X[i]=0; }
        
        // A alpha
        dummy=0; counter=0;
        for(j=0; j<cols; j++)
        {
            for(i=0; i<jcs[j+1]-jcs[j]; i++)
            {
                dummy = sr[counter]*alpha[counter];
                X[j] -= dummy;
                X[irs[counter]] += dummy;
                counter++;
            }
        }
        
        //D = Y- lambda Aalpha
        for(i=0; i<len; i++) {
            D[i]=(Y[i]-lambda*0.5*X[i])/eta;
            
            if (D[i]<0) Dproj[i]=0;
            else Dproj[i]=D[i];
        }
        
		//for(i=0;i<10;i++) mexPrintf("%1.8f ",Dproj[irs[i]]);
	 	//mexPrintf("\n");
        
        // project onto simplex
        Dproj = ProjKSimplex(Dproj,len,1);
        
        for (i=0;i<len;i++){
            Dproj[i]=(D[i]-Dproj[i])*eta;
            if (Dproj[i]<0) Dproj[i]=0;
        }
        
        //for(i=0;i<10;i++) mexPrintf("%1.8f ",Dproj[irs[i]]);
        //mexPrintf("\n");
        
        
        
        //start = clock();
        //ende = clock();
        //time2+=(float) (ende-start);
        
        // update beta and alpha
        counter=0;
        //tnew = (1 + sqrt(1+4*told*told))/2;
        //factor = (told-1)/tnew;
        factor=(double) iter/(iter+3);
        for(j=0; j<cols; j++)
        {
            Dcur=Dproj[j];
            for(i=0; i<jcs[j+1]-jcs[j]; i++)
            {
                temp=sr[counter]*( Dproj[irs[counter]] - Dcur);
                // update of beta
                betacur=alpha[counter] + lamByL*temp;
                // projection onto l_inf-cube
                if(betacur>1) betacur=1;
                else if(betacur<-1) betacur=-1;
                beta[counter]=betacur;
                // update of alpha
                alpha[counter] = betacur + factor*(betacur-betaold[counter]);
                counter++;
            }
        }
        
        // every now and then, compute objective to check if we have made 
        // sufficient progress
        if (computeObjective)
        {
            start=clock();
            
            // compute dual val
            for(i=0; i<len; i++) { X[i]=0; }
            
            dummy=0; counter=0;
            for(j=0; j<cols; j++)
            {
                for(i=0; i<jcs[j+1]-jcs[j]; i++)
                {
                    dummy = sr[counter]*beta[counter];
                    X[j] -= dummy;
                    X[irs[counter]] += dummy;
                    counter++;
                }
            }
            
            
            for(i=0; i<len; i++) {
                D[i]=(Y[i]-lambda*0.5*X[i])/eta;
                
                if (D[i]<0) Dproj2[i]=0;
                else Dproj2[i]=D[i];
                
            }
            
            Dproj2 = ProjKSimplex(Dproj2,len,1);
            
            for (i=0;i<len;i++){
                Dproj2[i]=(D[i]-Dproj2[i])*eta;
                if (Dproj2[i]<0) Dproj2[i]=0;
            }
            
            
            normD = 0;
            for(i=0; i<len; i++) {
                normD+= Dproj2[i]*Dproj2[i];
            }
            Fval = -normD;//  + C;
            
            //mexPrintf("Iteration: %i, Fval =%1.8f ",iter,Fval);
            
            
            primalVal_orig=0;
            for (i=0;i<len; i++)
            {
                primalVal_orig-=Dproj2[i]*Y[i];
            }
            
            //  mexPrintf("Iteration: %i, primalval_orig =%1.8f ",iter,primalVal_orig);
            
            start=clock();
            max_f=0;
            for (i=0;i<len;i++)
            {
                if (Dproj2[i]>max_f)
                    max_f=Dproj2[i];
            }
            // mexPrintf(" eta*max_f =%1.8f ",eta*max_f);
            
            
            primalVal=primalVal_orig+eta*max_f;
            ende=clock();
            time2+=(float) (ende-start);
            
            counter=0;
            primalVal2=0;
            for(j=0; j<cols; j++)
            {
                Dcur=Dproj2[j];
                for(i=0; i<jcs[j+1]-jcs[j]; i++)
                {
                    primalVal2+=fabs(sr[counter]*( Dproj2[irs[counter]] - Dcur));
                    counter++;
                }
            }
            
            //     mexPrintf(" primalVal2 =%1.8f \n",primalVal2);
            
            
     
            
            primalVal_orig=primalVal_orig+eta*max_f+lambda*0.5*primalVal2;
            primalVal=primalVal_orig+ 0.5*normD;// +C;

 			//mexPrintf(" primalVal_orig=%1.8f \n",primalVal_orig);
            

            if (normD >0) {
                primalVal_orig=primalVal_orig/ sqrt(normD);
            }
            else {
                primalVal_orig=0;
            }
            

		    //mexPrintf(" primalVal_orig=%1.8f \n",primalVal_orig);
            
            
            // check if we are better than previous best one
            if (primalVal_orig<primalValBest)
            {
                primalValBest=primalVal_orig;
                dualValBest=Fval;
                
                for(i=0; i<len; i++) {
                    DprojBest[i]=Dproj2[i];
                }
            }
            
            
            ende = clock();
            time4+=(float) (ende-start);
            
          // if (iter==0) break;
            
            //if (iter==MAXITER_temp)
            //{
            // check if converged
            if (primalVal<0)
                break;
            else {
                MAXITER_temp=MAXITER_temp*2;
                //told=1;tnew=1;
            }
            //}
            
            // output values
            if (debug) {
                if(iter % 10 ==0 || iter<10) {
                    mexPrintf("Iteration: %i, primalval_orig =%1.8f primalval =%1.8f dualval: %1.8f, primalVal-dualval=%1.8f\n",iter,primalVal_orig,primalVal,Fval,(primalVal-Fval)/fabs(Fval));
                }
            }
        }
        
        iter++;
    }
    
    // output final results
    if (debug) {
        mexPrintf("time1: %.2f   time2: %.2f  time3: %.2f   time4: %.2f\n",time1 / (float)CLOCKS_PER_SEC,time2 / (float)CLOCKS_PER_SEC,time3 / (float)CLOCKS_PER_SEC,time4 / (float)CLOCKS_PER_SEC);
        mexPrintf("FINAL Iteration: %i,  primalval_orig =%1.8f  primalval =%1.8f dualval: %1.8f, primalVal-dualval=%1.8f\n",iter,primalVal_orig,primalVal,Fval,(primalVal-Fval)/fabs(Fval));
    }
    
    // write output
    for(i=0; i<len; i++) { X[i]=Dproj2[i];}
    for(i=0; i<lenalpha; i++) { alphaout[i]=alpha[i];}
    PrimalObj[0]=primalVal;
    DualObj[0]=Fval;
    FinalIter[0]=iter;
    for(i=0; i<len; i++) { Xbest[i]=DprojBest[i];}
    PrimalObjBest[0]=primalValBest;
    DualObjBest[0]=dualValBest;
    delete betaold; delete beta; delete Dproj;delete Dproj2;delete DprojBest;delete D;
    
}



// projection on the k simplex
double* ProjKSimplex(double* f,int len, int K)
{
    int MAXIT=10000;
    int i;
    int* IX = new int[len];
    double sumxnew = 0; double sumxnew2;
    int counter=len; int counter2=0;
    for(i=0; i<len; i++)
    { IX[i]=i; sumxnew+=f[i];}
    int iteration=0;
    
    // main loop
    while(iteration<MAXIT)
    {
        counter2=0; sumxnew2=0;
        for(i=0;i<counter; i++)
        {
            f[IX[i]] = f[IX[i]] - (sumxnew-K)/counter;
            if(f[IX[i]]>0)
            { IX[counter2]=IX[i]; counter2++; sumxnew2+=f[IX[i]]; }
            else
                f[IX[i]]=0;
        }
        if(counter2==counter) break;
        counter=counter2; sumxnew=sumxnew2;
        iteration++;
    }
    // this should never happen
    if(iteration==MAXIT)
        mexPrintf("Bug in Projection: %i, Fval: %i\n",counter,counter2);
    delete IX;
    return f;
}

