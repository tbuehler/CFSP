function [clusters, ncut, feasible,lambda]= vol_cnstr_ncut_subset_penalty(W, k, hdeg, start, subset, gamma1, gamma2)
% Solves the normalized cut problem with generalized volume constraint and
% subset constraint, where both subset constraint and volume constraint
% have been incorporated into the objective as penalty term.
%
% Usage: [clusters, ncut,feasible,lambda]= vol_cnstr_ncut_subset_penalty(W,k,hdeg,start,subset,gamma1,gamma2)
%
% Input:
% W                 Weight matrix (full graph).
% k                 Upper bound
% hdeg              Generalized degrees used in constraint.
% start             Start vector
% subset            Indices of seed subset
% gamma1             Penalty parameter for seed constraint.
% gamma2             Penalty parameter for volume constraint.
%
% Output:
% clusters          Thresholded vector f yielding the best objective.
% ncut              Ncut value of resulting clustering.
% feasible          True if all constraints are fulfilled.
% lambda            Corresponding objective.

    %% check inputs
    assert(sum(diag(W))==0, 'Error! Diagonal entries of W are non-zero.');
    assert(k <= sum(hdeg), 'Error! Upper bound exceeds total volume.');
    assert(k >= sum(hdeg(subset)),  ... 
        'Error! Upper bound is smaller than volume of the seed subset.');
    assert(gamma1>=0,'Error! gamma1 cannot be negative.');
    assert(gamma2>=0,'Error! gamma2 cannot be negative.');
    assert(~issparse(hdeg),'Error! hdeg should not be sparse.');
    
    %% initialization
    start=abs(start);
    start=start/norm(start,2);
    f=start;
    
    %% some helpers 
    num=length(f);
    deg=full(sum(W,2));   
    [ix,jx,wval]=find(W);
    totVol = full(sum(sum(W)));
    inner_obj_all=0;

    %% compute objective
    f(subset) = max(f);
    [lambda,indvec]=functional_vol_cnstr_ncut_subset_penalty(f, ... 
        gamma1,gamma2,num,k,subset,deg,wval,ix,jx,hdeg, totVol);
    lambda_all=lambda;


    %% main loop
    %fprintf('gamma1=%.5g \t gamma2=%.5g \t lambda=%.5g \n', ... 
    %    full(gamma1), full(gamma2), lambda);
    converged=false;
    l=2;
    while (~converged && max(abs(f)) ~= 0) % avoid the zero starting 
                                           % point (on the subgraph).

        % solve inner problem
        [f_new, lambda_new,indvec_new,converged,primalobj]=  ... 
            solveInnerProblem(f,lambda,subset,k,W,wval,ix,jx,num, ... 
            deg,hdeg,gamma1,gamma2,indvec,totVol);


        % update variables
        lambda=lambda_new;
        f=f_new;
        indvec=indvec_new;

        % keep track of results
        inner_obj_all(l)=primalobj;
        lambda_all(l)=lambda_new;

        l=l+1;
    end
    %fprintf('  \n');

    
    %% perform optimal thresholding
    if (max(abs(f))==0) % in this case we take the start vector
        f=start;
    end
    [clusters, ncut_thresh, lambda_thresh,feasible] = ... 
        opt_thresh_vol_cnstr_ncut_subset_penalty(f, W, ... 
        deg,deg,hdeg, k, subset, gamma1, gamma2 );

    ncut = balanced_cut(W, sum(W,2), clusters);
    assert ( abs(ncut/sum(sum(W)) - ncut_thresh) <= 1e-8 );
   
    feasible=(sum(hdeg(clusters==1))<=k) && (sum(clusters(subset))==length(subset));
    
    %% compute objective
    lambda=functional_vol_cnstr_ncut_subset_penalty(clusters,gamma1, ... 
        gamma2,num,k,subset,deg,wval,ix,jx,hdeg, totVol);                  
  
    if gamma1==0 && gamma2==0
        assert( abs(ncut/sum(sum(W)) - lambda) <= 1e-8 );
    else
        assert(  lambda - lambda_thresh >= -1e-8 ); 
    end
end



%% solves the inner problem in RatioDCA
function [f_new, lambda_new,indvec_new,converged,primalobj]= ... 
    solveInnerProblem(f,lambda,subset,k,W,wval,ix,jx,num,deg,hdeg, ... 
    gamma1,gamma2,indvec,totVol)
       
    % set parameters
    MAXITER=800;
    MAXITER_start=50;
    eps1=1E-4;      
    debug=false;

    % compute subgradient of ncut balancing term
    Deg = spdiags(deg,0,size(f,1),size(f,1));
    Pf = f - (f'*deg/sum(deg));
    sg_bal = Deg*sign(Pf) - deg*sum( Deg*sign( Pf))/sum(deg);

    % subgradient for subset constraint
    gamma_subset = zeros(num,1);
    gamma_subset(subset) = gamma1;
    
    % compute constants
    c2=gamma2*indvec - gamma_subset - (lambda/2)*totVol*sg_bal;
    vec=-c2;
    c = gamma1*length(subset);   

    % solve inner problem with FISTA
    if gamma1 ~= 0
        [f_new1,alpha_new1,primalobj1,dualobj1,iter1,f_new,primalobj,dualobj]= ... 
            mex_ip_vol_cnstr_ncut_subset(W,2*vec,zeros(length(wval),1), ... 
            MAXITER,1E-8,4*max(sum(W.^2)),2*c, MAXITER_start,debug);
    else
        %[f_new,alpha_new,dualobj_new,iter1]=mex_ip_vol_cnstr_ncut(W,lambda*full(vec)/2, ... 
        %    zeros(length(wval),1),MAXITER,1E-8,4*sqrt(max(deg)));
        %primalobj=sum(0.5*wval.*abs(f_new(ix)-f_new(jx))) - lambda*full(vec)'*f_new;
        [f_new1,alpha_new1,primalobj1,dualobj1,iter1,f_new,primalobj,dualobj]= ... 
            mex_ip_vol_cnstr_ncut(W,2*vec,zeros(length(wval),1), ... 
            MAXITER,1E-8,4*max(sum(W.^2)), MAXITER_start,debug);
        
    end
    f_new = f_new/norm(f_new);

    % make some output
   % if ~isnan(max(f_new))
   %     display(['pt1 = ', num2str(gamma1*max(f_new)*length(subset) - ... 
   %         gamma1*sum(f_new(subset))), ' pt2 = ', num2str(indvec'*f_new)]);
   % end

    % in this case take the last one
    if primalobj==0
        f_new=f;
        lambda_new=lambda;
        indvec_new=indvec;
    else
        [lambda_new, indvec_new]=functional_vol_cnstr_ncut_subset_penalty(f_new, ... 
            gamma1,gamma2,num,k,subset,deg,wval,ix,jx,hdeg, totVol);
    end

    % check if converged
    diff=abs(lambda_new-lambda)/lambda;
    if (diff<eps1 || lambda_new == 0 || lambda==0)
        converged=true;
    else
        converged=false;
    end

end
    

