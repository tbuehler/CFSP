function [clusters, ncut, feasible,lambda]= vol_cnstr_ncut_subset(W,k,hdeg,start,subset,gamma)

    %% check inputs
    assert(sum(diag(W))==0, 'Diagonal entries of W are non-zero. ');
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

    inner_obj_all=0;

    num=length(f);
    deg=full(sum(W,2));
    [ix,jx,wval]=find(W);
    
    %% restrict the graph
    ind_rest=setdiff((1:num)',subset);  % V\J
    W_rest=W(ind_rest,ind_rest);
    [ix_rest,jx_rest,wval_rest]=find(W_rest);
    degJ=full(sum(W(ind_rest,subset),2)); 
    f=f(ind_rest);


    %% compute objective
    totVol = full(sum(sum(W)));
    totVol_rest = full(sum(sum(W(ind_rest, :))));
    ftemp=zeros(num,1);
    ftemp(ind_rest)=f;
    ftemp(subset)=max(f);
    [lambda,indvec]=functional_vol_cnstr_ncut(ftemp,gamma,num,k,deg, ... 
        wval,ix,jx,hdeg, totVol);
    indvec=indvec(ind_rest);


    % compute again differentily for check
    kprime=k-sum(hdeg(subset));
    [lambda_direct,indvec_direct]=functional_vol_cnstr_ncut_subset_direct(f, ... 
        gamma,length(f),kprime,deg(ind_rest),wval_rest,ix_rest,jx_rest, ... 
        hdeg(ind_rest), totVol_rest,degJ);
    %assert(lambda_direct==lambda);
    assert(sum(indvec_direct~= indvec)==0);


    lambda_all=lambda;

    %% main loop
    fprintf('gamma=%.5g \t lambda=%.5g \n',full(gamma), lambda);
    converged=false;
    l=2;
    while (~converged && max(abs(f)) ~= 0) % avoid the zero starting point 


        % solve inner problem
        [f_new, lambda_new,indvec_new,converged]= solveInnerProblem(f, ... 
            lambda,subset,k,W_rest,wval_rest,ix_rest,jx_rest,wval,ix,jx, ... 
            num,deg,hdeg,ind_rest,gamma,indvec,degJ,totVol,totVol_rest);

        % update variables
        lambda=lambda_new;
        f=f_new;
        indvec=indvec_new;

        % keep track of results
        inner_obj_all(l)=primalobj;
        lambda_all(l)=lambda_new;

        l=l+1;
    end

    fprintf('  \n');

    %% perform optimal thresholding
    if max(abs(f)~=0)
        f_new_full = zeros(num,1);
        f_new_full(subset) = max(f_new)+1e-12;
        f_new_full(ind_rest) = f_new;
    else
        f_new_full=start;
    end
    
    [clusters, ncut_thresh, lambda_thresh,feasible] =  ... 
        opt_thresh_vol_cnstr_ncut_subset_penalty(f_new_full, W, ... 
        deg,deg,hdeg, k, subset, gamma ,gamma);
     
    ncut = balanced_cut(W, sum(W,2), clusters);
    assert ( abs(ncut - ncut_thresh) <= 1e-8 );

    %% check if partition (S=seed set, S^c) is optimal
    if ~isempty(subset)

        clusters_subset = zeros(num,1);
        clusters_subset(subset) = 1;
        ncut_subset = balanced_cut(W, sum(W,2), clusters_subset );
        feasible_subset = k>=sum(hdeg(subset)); 

        if feasible_subset && ncut > ncut_subset
            clusters = clusters_subset;
            ncut = ncut_subset;
            feasible=feasible_subset;
        end
    end
  

    %% compute objective
    [lambda,indvec]=functional_vol_cnstr_ncut(clusters,gamma,num,k,... 
        deg,wval,ix,jx,hdeg, totVol);
    % compute again differentily
    [lambda_direct,indvec_direct]=functional_vol_cnstr_ncut_subset_direct(clusters(ind_rest), ... 
        gamma,length(f),kprime,deg(ind_rest),wval_rest,ix_rest,jx_rest,hdeg(ind_rest), totVol_rest,degJ);
    %assert(lambda_direct==lambda);
    indvec=indvec(ind_rest);
    assert(sum(indvec~=indvec_direct)==0);


    if gamma==0
        assert( abs(ncut - lambda) <= 1e-8 );
    else
        assert(  lambda - lambda_thresh >= -1e-8 );
    end

end


% solves the convex inner problem in RatioDCA
function [f_new, lambda_new,indvec_new,converged]= solveInnerProblem(f, ... 
    lambda,subset,k,W_rest,wval_rest,ix_rest,jx_rest,wval,ix,jx,num, ... 
    deg,hdeg,ind_rest,gamma,indvec,degJ,totVol,totVol_rest)

    % set parameters
    MAXITER=5000;
    MAXITER_start=50;

    indmax=zeros(size(f,1),1);
    indmax(f==max(f))=1;

    Deg = spdiags(deg(ind_rest),0,size(f,1),size(f,1));
    Pf = f - (f'*deg(ind_rest)/sum(deg(ind_rest)));
    sg_bal = Deg*sign(Pf) - deg(ind_rest)*sum( Deg*sign( Pf))/sum(deg(ind_rest));

    % sign mistake corrected
    vec=-gamma*indvec + degJ + (lambda/2)*totVol_rest*sg_bal + ... 
        sum(deg(subset))*totVol_rest*lambda*indmax - lambda* sum(deg(subset)) *deg(ind_rest);
    c = sum(deg(subset));   %c is the coefficient of the term: f_n. That is max(f).

    %%% f_new gives the best objective! We are not using f_new1.
    [f_new1,alpha_new1,primalobj1,dualobj1,iter1,f_new,primalobj,dualobj]= ... 
        mex_ip_vol_cnstr_ncut_subset(W_rest,2*vec,zeros(length(wval_rest),1), ... 
        MAXITER,1E-8,4*max(sum(W_rest.^2)),2*c,MAXITER_start);
    f_new = f_new/norm(f_new);

    if primalobj==0
        f_new=f;
        lambda_new=lambda;
        indvec_new=indvec;
    else
        ftemp=zeros(num,1);
        ftemp(ind_rest)=f_new;
        ftemp(subset)=max(f_new);
        [lambda_new, indvec_new]=functional_vol_cnstr_ncut(ftemp,gamma, ... 
            num,k,deg,wval,ix,jx,hdeg, totVol);
        indvec_new=indvec_new(ind_rest);


        % compute again differentily
        kprime=k-sum(hdeg(subset));
        [lambda_direct,indvec_direct]=functional_vol_cnstr_ncut_subset_direct(f_new, ... 
            gamma,length(f_new),kprime,deg(ind_rest),wval_rest,ix_rest, ... 
            jx_rest,hdeg(ind_rest), totVol_rest,degJ);
        %assert(lambda_direct==lambda);
        assert(sum(indvec_direct~= indvec_new)==0);
    end


    % check if converged
    diff=abs(lambda_new-lambda)/lambda;
    if (diff<eps1 || lambda_new == 0 || lambda==0)
        converged=true;
    else
        converged=false;
    end

end

