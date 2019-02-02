function [clusters, ncut,feasible,lambda]= vol_cnstr_ncut_subset_direct(W,k,hdeg,start,subset,gamma)
% Solves the normalized cut problem with generalized volume constraint and
% subset constraint, where the subset constraint has been directly 
% incorporated into the objective (leading to a problem of lower dimension) 
% and the volume constraint has been incorporated into the objective as 
% penalty term.
%
% Usage: [clusters, ncut,feasible,lambda]= vol_cnstr_ncut_subset_direct(W,k,hdeg,start,subset,gamma)
%
% Input:
% W                 Weight matrix (full graph).
% k                 Upper bound
% hdeg              Generalized degrees used in constraint.
% start             Start vector
% subset            Indices of seed subset
% gamma             Penalty parameter for volume constraint.
%
% Output:
% clusters          Thresholded vector f yielding the best objective.
% ncut              Ncut value of resulting clustering.
% feasible          True if all constraints are fulfilled.
% lambda            Corresponding objective.

    %%  check inputs
    if k>sum(hdeg)
        k=sum(hdeg);
        fprintf('Setting k to %f\n',k);
    end
    assert(k>=sum(hdeg(subset)),'Error! Problem is unfeasible.');
    assert(gamma>=0,'Error! gamma cannot be negative.');

    %% initialization
    start=abs(start);
    start=start/norm(start,2);
    f=start;
 
    
    %% restrict graph
    num=length(f);
    ind_rest=setdiff((1:num)',subset);
    W_rest=W(ind_rest,ind_rest);
    [ix_rest, jx_rest, wval_rest] = find(W_rest);
    f=f(ind_rest);
    kprime=k-sum(hdeg(subset));
    degJ=full(sum(W(ind_rest,subset),2));
    totVol_rest = full(sum(sum(W(ind_rest, ind_rest))));
    deg=full(sum(W,2));
    

    %% evaluate objective
    [lambda,indvec]=functional_vol_cnstr_ncut_subset_direct(f,gamma,... 
        num-length(subset),k,deg(ind_rest),wval_rest,ix_rest,jx_rest, ... 
        hdeg(ind_rest), totVol_rest, degJ);
    inner_obj_all=0;
    lambda_all=lambda;


    %% main loop
    %fprintf('gamma=%.5g \t lambda=%.5g \n',full(gamma), lambda);
    converged=false;
    l=2;
    while (~converged)
 
        % solve inner problem
        [f_new, lambda_new,indvec_new,converged,primalobj]= ... 
            solveInnerProblem(f,lambda,subset,k,W_rest,wval_rest,... 
            ix_rest,jx_rest,num,deg,hdeg,ind_rest,gamma,indvec,degJ,totVol_rest);


        % update variables
        lambda=lambda_new;
        f=f_new;
        indvec=indvec_new;

        % store results
        inner_obj_all(l)=primalobj;
        lambda_all(l)=lambda_new;

        l=l+1;
    end
    %fprintf('  \n');

    
    %% Perform optimal thresholding
    [clusters_temp, ncut, lambda, feasible] = ...
        opt_thresh_vol_cnstr_ncut_subset_direct(f_new,W_rest, ...
        deg(ind_rest), hdeg(ind_rest),kprime, gamma,sum(deg(subset)),degJ);
    
    %% construct vector with respect to original graph
    clusters=zeros(num,1);
    clusters(subset)=1;
    clusters(ind_rest)=clusters_temp;

    %% check if partition (S=seed set, S^c) is optimal
    if ~isempty(subset)
        clusters_subset = zeros(num,1);
        clusters_subset(subset) = 1;
        ncut_subset = balanced_cut(W, sum(W,2), clusters_subset );
        feasible_subset = k>=sum(hdeg(subset)); 
        lambda_subset= functional_vol_cnstr_ncut_subset_direct( ... 
            zeros(num-length(subset),1),gamma,num-length(subset),k, ... 
            deg(ind_rest),wval_rest,ix_rest,jx_rest,hdeg(ind_rest), ... 
            totVol_rest, degJ);
        
        if feasible_subset && ncut > ncut_subset

            clusters = clusters_subset;
            ncut = ncut_subset;
            feasible=feasible_subset;
            lambda= lambda_subset;
        end
    end

end



%% solves the inner problem in RatioDCA
function [f_new, lambda_new,indvec_new,converged,primalobj]=  ... 
    solveInnerProblem(f,lambda,subset,k,W_rest,wval_rest,ix_rest,jx_rest, ... 
    num,deg,hdeg,ind_rest,gamma,indvec,degJ,totVol_rest)
        
    % set parameters
    MAXITER=1600;
    MAXITER_start=50;
    eps1=1E-4;
    debug=false;

	% compute subgradient of subset penalty
    indmax=zeros(size(f,1),1);
    indmax(f==max(f))=1;
    
    % compute subgradient of ncut balancing term
    Deg = spdiags(deg(ind_rest),0,size(f,1),size(f,1));
    Pf = f - (f'*deg(ind_rest)/sum(deg(ind_rest)));
    sg_bal = Deg*sign(Pf) - deg(ind_rest)*sum( Deg*sign( Pf))/sum(deg(ind_rest));
    
    % compute constants
    c2=gamma*indvec - degJ - (lambda/2)*totVol_rest*sg_bal ... 
        - sum(deg(subset))*totVol_rest*lambda*indmax ... 
        + lambda* sum(deg(subset)) *deg(ind_rest);
    vec=-c2;
    c = sum(degJ);
    
    % solve inner problem with FISTA
    [f_new1,alpha_new1,primalobj1,dualobj1,iter1,f_new,primalobj,dualobj]= ... 
        mex_ip_vol_cnstr_ncut_subset(W_rest,2*vec,zeros(length(wval_rest),1), ... 
        MAXITER,1E-8,4*max(sum(W_rest.^2)),2*c,MAXITER_start,debug);
    f_new = f_new/norm(f_new);
    
    % in this case take last iterate
    if primalobj==0
        f_new=f;
        lambda_new=lambda;
        indvec_new=indvec;
    else
        [lambda_new, indvec_new]=functional_vol_cnstr_ncut_subset_direct(f_new, ... 
            gamma,num-length(subset),k,deg(ind_rest),wval_rest,ix_rest, ...
            jx_rest,hdeg(ind_rest), totVol_rest, degJ);
    end
    
    % check if converged
    diff=abs(lambda_new-lambda)/lambda;
    if (diff<eps1)
        converged=true;
    else
        converged=false;
    end

end



    

