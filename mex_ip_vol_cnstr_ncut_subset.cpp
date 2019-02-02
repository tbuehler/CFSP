#include <math.h>
#include "mex.h"
#include "matrix.h"
#include <time.h>
#include "float.h"

double* ProjKSimplex(double*, int,int);



void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]) {
  // Test number of parameters.
  if (nrhs != 9 || nlhs != 8) {
    mexWarnMsgTxt("Usage: [X,alpha,PrimalObj,DualObj,iter,Xbest,PrimalObjBest,DualObjBest]=mex_ip_vol_cnstr_ncut_subset(W,Y,alpha,MAXITER,EPS,MaxSumSquaredWeights,eta,MAXITER_temp,debug)");
    return;
  }
  
  // get important parameters
  mwSize rows = mxGetM(prhs[0]);      // number of rows of W
  mwSize cols = mxGetN(prhs[0]);      // number of columns of W (should be the same)
  mwSize len  = mxGetM(prhs[1]);      // the desired output
  mwSize lenalpha = mxGetM(prhs[2]);  // alpha
  
  // check if input has right format
  if(!mxIsSparse(prhs[0])) { 
    mexWarnMsgTxt("Matrix is not sparse");
    return;
  }
  
  if(rows!=cols){
    mexWarnMsgTxt("Sparse matrix is not square");
    return;
  }

  if(rows!=len){
    mexWarnMsgTxt("Length of the vector is not the same as the number of the rows of the sparse matrix");
    return;
  }
    
  // Create output array and compute values
  double* sr = mxGetPr(prhs[0]);     // get values
  mwIndex* irs = mxGetIr(prhs[0]);   // get row
  mwIndex* jcs = mxGetJc(prhs[0]);   // get columns
  
  double* Y = mxGetPr(prhs[1]);		// Y
  double* alpha = mxGetPr(prhs[2]);     // alpha

  mwIndex MAXITER =  mxGetScalar(prhs[3]); 
  double EPS = mxGetScalar(prhs[4]); //for criterion fabs(dualval)<EPS
  double MaxSumSquaredWeights = mxGetScalar(prhs[5]); 
  double eta = mxGetScalar(prhs[6]); 
  mwIndex MAXITER_temp =  mxGetScalar(prhs[7]); 
  bool debug=(bool) mxGetScalar(prhs[8]);

  if(MaxSumSquaredWeights<=0){
	  mexWarnMsgTxt("Lipschitz constant has to be positive");
    return;
  }

  // allocate memory for output
  plhs[0] = mxCreateDoubleMatrix(len,1,mxREAL);     // output vector 
  plhs[1] = mxCreateDoubleMatrix(lenalpha,1,mxREAL); // dual variable 
  plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);       // objective value
  plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);       // dual objective
  plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL);       // final iteration value
  plhs[5] = mxCreateDoubleMatrix(len,1,mxREAL); 	// best output vector
  plhs[6] = mxCreateDoubleMatrix(1,1,mxREAL);       // best primal objective
  plhs[7] = mxCreateDoubleMatrix(1,1,mxREAL);       // best dual objective
  
  // pointers for output
  double* X = mxGetPr(plhs[0]); 
  double* Z = mxGetPr(plhs[1]);
  double* PrimalObj = mxGetPr(plhs[2]);
  double* DualObj = mxGetPr(plhs[3]);
  double* FinalIter = mxGetPr(plhs[4]);
  double* Xbest = mxGetPr(plhs[5]); 
  double* PrimalObjBest = mxGetPr(plhs[6]);
  double* DualObjBest = mxGetPr(plhs[7]);

  mwIndex counter=0,i,j,iter=0;
  double Dcur,betacur,factor;
  double dummy,normD,Fval;
  double* dummyPointer;
  double L=1/MaxSumSquaredWeights;

  double* D =new double[len];
  double* Dproj =new double[len];
  double* Dproj2 =new double[len];
  double* DprojBest =new double[len];
  double* beta    = new double[lenalpha];
  double* betaold = new double[lenalpha];
  for(i=0; i<lenalpha; i++) { beta[i]=0; }
  for(i=0; i<lenalpha; i++) { betaold[i]=0; }
  double temp,primalVal2,primalVal_orig;
  double primalVal=0;
  double dualValBest= -DBL_MAX;
  
  double C=0;
  for(i=0;i<len;i++){C+=Y[i]*Y[i];}
  C=0.5*C;
  
  double primalValBest=0;double max_f;
  bool computeObjective;

  for(i=0; i<len; i++) DprojBest[i]=0;

  clock_t start,ende;
  float time2=0, time4=0;

//  int MAXITER_temp=50;
 
  Fval=EPS+1;
  while(iter<MAXITER){
    // defines strategy when to check if converged
	if (iter==MAXITER-1 || iter==MAXITER_temp )     {
        computeObjective=true;
    } else {   
		computeObjective =false;
    }
      
    // exchange beta and betaold
	dummyPointer=beta; beta=betaold; betaold=dummyPointer;
 
	// initialize X with zeros 
    for(i=0; i<len; i++) { X[i]=0; }

    // A alpha
	dummy=0; counter=0;
    for(j=0; j<cols; j++) 
    {   
	   for(i=0; (unsigned int) i<jcs[j+1]-jcs[j]; i++)
	   {  
         dummy = sr[counter]*alpha[counter];
	     X[j] -= dummy;
	     X[irs[counter]] += dummy;
         counter++;
	   }
	} 
 
    //D = Y - Aalpha
    for(i=0; i<len; i++) { 
        D[i]=(Y[i]-X[i])/eta;
        
		//Dproj[i]=(Y[i]-lambda*X[i])/eta; 
        //if (Dproj[i]<0) {Dproj[i]=0;}
        
        
        if (D[i]<0) Dproj[i]=0;
        else Dproj[i]=D[i];
        
    }
    
    //OptSimplex(Dproj,eta, C,len);
    Dproj = ProjKSimplex(Dproj,len,1);
    
    // This projection is for the gradient.
    for (i=0;i<len;i++){
        Dproj[i]=(D[i]-Dproj[i])*eta;
        if (Dproj[i]<0) Dproj[i]=0;
    }
   
      
    // primalval2 = lambda*sum(w_ij|f_i-f_j|)
    // update beta  and alpha 
    counter=0;
    
	//start = clock();
	//ende = clock();
	//time2+=(float) (ende-start);
    
    //primalVal2=0;
    factor=(double) iter/(iter+3);
    for(j=0; j<cols; j++) 
    {   
       Dcur=Dproj[j];
	   for(i=0; (unsigned int) i<jcs[j+1]-jcs[j]; i++)
	   {             
	      temp=sr[counter]*( Dproj[irs[counter]] - Dcur);
        //  if (computeObjective){
		//	primalVal2+=fabs(temp);
		//  }
          // update of beta
          betacur=alpha[counter] + L*temp;
          // projection onto l_inf-cube 
          if(betacur>1) betacur=1;
	      else if(betacur<-1) betacur=-1;
		  beta[counter]=betacur;
		  // update of alpha
          alpha[counter] = betacur + factor*(betacur-betaold[counter]);
          counter++;
	   }	  
    }
    

	
	if (computeObjective)
	{
		start=clock();
       
		// compute dual val
		for(i=0; i<len; i++) { X[i]=0; }

		dummy=0; counter=0;
		for(j=0; j<cols; j++) 
		{   
			for(i=0; (unsigned int) i<jcs[j+1]-jcs[j]; i++)
			{  
				dummy = sr[counter]*beta[counter];
		        X[j] -= dummy;
				X[irs[counter]] += dummy;
				counter++;
			}
		} 
	        
        
        for(i=0; i<len; i++) { 
            D[i]=(Y[i]-X[i])/eta;
            
        
            if (D[i]<0) Dproj2[i]=0;
            else Dproj2[i]=D[i];
        
        }
    
        Dproj2 = ProjKSimplex(Dproj2,len,1);
    
        for (i=0;i<len;i++){
            Dproj2[i]=(D[i]-Dproj2[i])*eta;
            if (Dproj2[i]<0) Dproj2[i]=0;
     }
        
        
        
        
        primalVal_orig=0;
		for (i=0;i<len; i++)
		{
            primalVal_orig-=Dproj2[i]*Y[i];
		}
        
       
        start=clock();
        max_f=0;
        for (i=0;i<len;i++)
        {
            if (Dproj2[i]>max_f)
            	max_f=Dproj2[i];
        }
         
        
        primalVal=primalVal_orig+eta*max_f;
        ende=clock();
        time2+=(float) (ende-start);
       
        counter=0;
        primalVal2=0;
        for(j=0; j<cols; j++) 
        {   
            Dcur=Dproj2[j];
            for(i=0; (unsigned int) i<jcs[j+1]-jcs[j]; i++)
            {  
                primalVal2+=fabs(sr[counter]*( Dproj2[irs[counter]] - Dcur));
	            counter++;
            }	  
        }
        
  
 
        normD = 0; 
        for(i=0; i<len; i++) { 
        	normD+= Dproj2[i]*Dproj2[i];
		}
		Fval = -0.5*normD;//  + C;
   
        primalVal_orig=primalVal_orig+eta*max_f+primalVal2;
	    primalVal=primalVal_orig+ 0.5*normD;// +C;
        
        Fval=0.5*Fval;
        primalVal_orig=0.5*primalVal_orig;
        primalVal=0.5*primalVal;
       
        
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
  
        if (debug) {
            if(iter % 10 ==0 || iter<10 || iter==MAXITER-1) { 
            mexPrintf("Iteration: %i, primalval_orig =%1.8f primalval =%1.8f dualval: %1.8f, primalVal-dualval=%1.8f\n",iter,primalVal_orig,primalVal,Fval,(primalVal-Fval)/fabs(Fval));
            }
        }

        
        if (primalVal<0)
                break;
        else
            MAXITER_temp=MAXITER_temp*2;
        
       
        
                
	}
   
	iter++;
 }

  
//  mexPrintf("time1: %.2f   time2: %.2f  time3: %.2f   time4: %.2f\n",time1 / (float)CLOCKS_PER_SEC,time2 / (float)CLOCKS_PER_SEC,time3 / (float)CLOCKS_PER_SEC,time4 / (float)CLOCKS_PER_SEC);
// mexPrintf("FINAL Iteration: %i,  primalval_orig =%1.8f  primalval =%1.8f dualval: %1.8f, primalVal-dualval=%1.8f\n",iter,primalVal_orig,primalVal,Fval,(primalVal-Fval)/fabs(Fval));

 for(i=0; i<len; i++) { X[i]=Dproj2[i];}
 for(i=0; i<lenalpha; i++) { Z[i]=alpha[i];}
 PrimalObj[0]=primalVal;
 DualObj[0]=Fval;
 FinalIter[0]=iter;
 for(i=0; i<len; i++) { Xbest[i]=DprojBest[i];}
  PrimalObjBest[0]=primalValBest;
 DualObjBest[0]=dualValBest;
 delete[] betaold; delete[] beta; delete[] Dproj;delete[] Dproj2;delete[] DprojBest;delete[] D;
	  
	
}








// projection on the k simplex
double* ProjKSimplex(double* f,int len, int K)
{
  int i;
  //double* xnew = new double[len];
  int* IX = new int[len];
  double sumxnew = 0; double sumxnew2;
  int counter=len; int counter2=0;
  for(i=0; i<len; i++)
   { IX[i]=i; sumxnew+=f[i];}
  int iteration=0;
  while(iteration<10000)
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
  if(iteration==10000)
	  mexPrintf("Bug in Projection: %i, Fval: %i\n",counter,counter2);
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

