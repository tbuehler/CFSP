function [dc,maxdens,lambda] = ... 
    cnstr_maxdens_subset_single_run(W,k1,k2,gdeg,start,subset,gamma)
% Performs one run with initialization given by start of the RatioDCA for
% the constrained maximum density problem.
%
% Usage: [dc,maxdens, lambda_new] ... 
%    = cnstr_maxdens_subset_single_run(W, k1,k2, gdeg,start, subset, gamma)
%
% Input: 
% W         The weight matrix.
% k1        Lower bound
% k2        Upper bound
% gdeg      Generalized degrees
% subset    Index of seed set
% gamma     Penalty parameter gamma.
%
% Output:
% dc        Indicator vector of the found community
% maxdens   The corresponding density
% lambda    The value of the objective
% 


    debug=false;

    % do some checks
    if k2>size(W,1)
        k2=size(W,1);
        fprintf('Setting k2 to %d.\n',k2);
    end
    assert(k1>=1,'Error: k1 has to be at least 1');
    assert(k1<=k2,'Error:k1 has to be less or equal to k2');

    %initialization
    start=abs(start);
    start=start/norm(start,2);
    f=start;

    inner_obj_all=0;

    num=length(f);
    deg=sum(W,2);
   
    % restrict graph
    ind_rest=setdiff((1:num)',subset);
    W_rest=W(ind_rest,ind_rest);
    [ix_rest,jx_rest,wval_rest]=find(W_rest);

    degJ=sum(W(ind_rest,subset),2);
    deg_tilde=deg(ind_rest)+degJ;

    assocJ=sum(sum(W(subset,subset)));
    gvolJ=sum(gdeg(subset));

    f=f(ind_rest);
    gdeg_rest=gdeg(ind_rest);
    
    assert(nnz(W_rest)==length(wval_rest));
    
    k1prime=k1-length(subset);
    k2prime=k2-length(subset);

    % if no gamma is given, we compute it using a subset of V with
    % cardinality k2
    if nargin<7
        ixx=ix_rest(1);
        if ixx-1+k2prime<=length(ind_rest)
            indgam=ixx:1:ixx-1 +k2prime;
        else
            indgam= num-k2prime+1:num;
        end
        indtemp=[subset; ind_rest(indgam)];
        gamma=sum(deg)*sum(gdeg(indtemp))/sum(sum(W(indtemp,indtemp)));

    end

    % compute the value of the objective
    [lambda,indvec,indvec1]=compute_lambda(f,gamma,k1,k2,deg(ind_rest),gvolJ,degJ,assocJ,wval_rest,ix_rest,jx_rest,gdeg_rest);

    
    lambda_all=lambda;
    
    % make some output
    if (debug)
        fprintf('gamma=%.5g \t lambda=%.5g \n',full(gamma), lambda);
    end
    converged=false;
    l=2;
    while (~converged)
        
         % solve inner problem
        [f_new, lambda_new, indvec_new, indvec_new1,primalobj,dualobj,iter1,converged] = ... 
            solveInnerProblem(f,indvec1,indvec,gamma,deg,W_rest,wval_rest,ix_rest,jx_rest,gdeg_rest,ind_rest,lambda,deg_tilde,assocJ,k1prime,k2prime,gvolJ,degJ,debug);
               
        % make some output
        if (debug)
            fprintf('lambda=%.5g  primalobj=%.5g dualobj=%.5g iter=%d \n', lambda_new,primalobj,dualobj,iter1);
        end     
             
      
        % update variables
        lambda=lambda_new;
        f=f_new;
        indvec=indvec_new;
        indvec1=indvec_new1;

        inner_obj_all(l)=primalobj;
        lambda_all(l)=lambda_new;

        l=l+1;
    end

    % make some output
    if (debug)
        fprintf('  \n');
    end

    % perform optimal thresholding
    dc_temp = opt_thresh_cnstr_maxdens_subset(W_rest, f_new, gdeg(ind_rest),max(1,k1prime), k2prime,gamma, assocJ,sum(gdeg(subset)),W(ind_rest,subset));
    dc=ones(size(W,1),1);
    dc(ind_rest) = dc_temp>0;
    maxdens = full(sum(sum(W(dc==1, dc==1)))/ sum(gdeg(dc>0)));


end



% solves the convex inner problem in RatioDCA
function [f_new, lambda_new, indvec_new, indvec_new1,primalobj,dualobj,iter1,converged,primalval2,dualval2] = ... 
    solveInnerProblem(f,indvec1,indvec,gamma,deg,W_rest,wval_rest,ix_rest,jx_rest,gdeg_rest,ind_rest,lambda,deg_tilde,assocJ,k1prime,k2prime,gvolJ,degJ,debug)

    % some parameters
    eps1=1E-4;
    MAXITER=1600;%800;

    % compute constants
    indmax=zeros(size(f,1),1);
    indmax(f==max(f))=1;
    c2=gdeg_rest -gamma*indvec1+gamma*indvec -lambda*deg_tilde;
    c1=k1prime*gamma+gvolJ-lambda*assocJ;
   
    % Lipschitz constant
    L = 2*lambda^2*max(sum(W_rest.^2));
    num=nnz(W_rest);
     
    % solve inner problem
    [X1,alpha1,PrimalObj1,DualObj1,iter1,f_new,primalobj,dualobj]= ... 
            mex_IP_maxdens_FISTA(W_rest,-c2,zeros(num,1),MAXITER,L,lambda,c1,debug);
   
    % recompute dual objective (for testing)    
    if (debug)
        dualval2=-norm(f_new)^2;
    else 
        dualval2=inf;
    end
    
    % renormalize
    f_new = f_new/norm(f_new);
    
    % recompute primal objective (for testing)
    if (debug)
        primalval2= lambda*0.5*sum(wval_rest.*abs(f_new(ix_rest)-f_new(jx_rest))) + c1*max(f_new) + f_new'*c2;
    else
        primalval2=inf;
    end

    % compute objective
    if primalobj==0
        f_new=f;
        lambda_new=lambda;
        indvec_new=indvec;
        indvec_new1=indvec1;
    else
        [lambda_new,indvec_new,indvec_new1]= ... 
            compute_lambda(f_new,gamma,k1prime,k2prime,deg(ind_rest),gvolJ,degJ,assocJ,wval_rest,ix_rest,jx_rest,gdeg_rest);
    end
    
    % check if converged
    diff=abs(lambda_new-lambda)/lambda;
    if (diff<eps1)
        converged=true;
    else
        converged=false;
    end
        
end


