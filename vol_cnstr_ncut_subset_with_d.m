function [clusters, ncut, objective]= vol_cnstr_ncut_subset_with_d(W,k,gdeg,start,subset,gamma, MAXITER)

    assert(sum(diag(W))==0, 'Diagonal entries of W are non-zero. ');

    if k>sum(gdeg)
        k=sum(gdeg);
        fprintf('Setting k2 to %f\n',k);
    end

    %initialization
    start=abs(start);
    start=start/norm(start,2);
    f=start;

    inner_obj_all=0;

    num=length(f);
    deg=sum(W,2);   % deg is used for finding cuts etc. So this is not general vertex weights/degrees
    [ix,jx,wval]=find(W);
    eps1=1E-12;      

    ind_rest=setdiff((1:num)',subset);  % V\J
    W_rest=W(ind_rest,ind_rest);
    [~,~,wval_rest]=find(W_rest);

    degS=sum(W(ind_rest,subset),2); % in the paper it is  d^(J)

    %assocJ=sum(sum(W(subset,subset)));

    f=f(ind_rest);
    %gdeg_rest=gdeg(ind_rest); % we do not need gdeg_rest because we
    %compute the subgradient using the subgradient of the full problem.

    totVol = sum(sum(W)); % This is used only in the computation of lambda and since lambda is calculated on the whole graph, we use the total volume including the subset.
    totVol_rest = sum(sum(W(ind_rest, :))); % This is used in the formation of the inner problem.
    ftemp=zeros(num,1);
    ftemp(ind_rest)=f;
    ftemp(subset)=max(f);
    [lambda,indvec]=vol_cnstr_ncut_subset_functional(ftemp,gamma,num,k,deg,wval,ix,jx,gdeg, totVol);
    indvec=indvec(ind_rest);

    lambda_all=lambda;


    fprintf('gamma=%.5g \t lambda=%.5g \n',full(gamma), lambda);
    converged=false;
    l=2;

    MAXITER_start = 50;
    while (~converged && max(abs(f)) ~= 0) % avoid the zero starting point (on the subgraph) or solution.

        indmax=zeros(size(f,1),1);
        indmax(f==max(f))=1;

        Deg = spdiags(deg(ind_rest),0,size(f,1),size(f,1));
        Pf = f - (f'*deg(ind_rest)/sum(deg(ind_rest)));
        sg_bal = Deg*sign(Pf) - deg(ind_rest)*sum( Deg*sign( Pf))/sum(deg(ind_rest));

        % sign mistake corrected
        vec=-gamma*indvec + degS + (lambda/2)*totVol_rest*sg_bal + sum(deg(subset))*totVol_rest*lambda*indmax - lambda* sum(deg(subset)) *deg(ind_rest);
        c = sum(deg(subset));   %c is the coefficient of the term: f_n. That is max(f).

        %%% f_new gives the best objective! We are not using f_new1.
        [~,~,~,~,iter1,f_new,primalobj,dualobj]=mex_ip_vol_cnstr_ncut_subset(W_rest,2*vec,zeros(length(wval_rest),1),MAXITER,1E-8,4*max(sum(W_rest.^2)),lambda,1E-8,2*c,MAXITER_start);       
        f_new = f_new/norm(f_new);

        if primalobj==0
            f_new=f;
            lambda_new=lambda;
            indvec_new=indvec;
        else
            ftemp=zeros(num,1);
            ftemp(ind_rest)=f_new;
            ftemp(subset)=max(f_new);
            [lambda_new, indvec_new]=vol_cnstr_ncut_subset_functional(ftemp,gamma,num,k,deg,wval,ix,jx,gdeg, totVol);        
            indvec_new=indvec_new(ind_rest);
        end

        if ~isnan(lambda_new)
            display('checking if there is a descent');
            %&& primalobj < 0 
            %assert( lambda_new <= lambda );                
            %assert(  lambda - lambda_new >= -1e-6 ); % the difference should always be non-negative, but because of numerical accuracies, we allow some eps diff!
            if lambda_new > lambda
                if MAXITER_start >= MAXITER || iter1 >= MAXITER
                %if iter1 >= MAXITER
                    f_new = f;  % We got a positive value and MAXITER exhausted.. so return the last iterate!
		    lambda_new=lambda;
                    indvec_new=indvec;
                    break;
                end
                MAXITER_start = iter1*2;
            end
        end
            
        % check if converged
        diff=abs(lambda_new-lambda)/lambda;
        if (diff<eps1 || lambda_new == 0 || lambda==0)
            converged=true;
        end

        % update variables
        lambda=lambda_new;
        f=f_new;
        indvec=indvec_new;

        inner_obj_all(l)=primalobj;
        lambda_all(l)=lambda_new;

        l=l+1;
    end

    fprintf('  \n');

    if max(abs(f)~=0)
        f_new_full = zeros(num,1);
        f_new_full(subset) = max(f_new)+1e-12;
        f_new_full(ind_rest) = f_new;
        %[clusters, ncut1, ncut2] = opt_thresh_vol_cnstr_ncut_subset(f_new_full, W, gdeg', 1, k, gamma );
        %[clusters, ncut_thres, lambda_thres] = opt_thresh_vol_cnstr_ncut_subset_fused(f_new_full, W, gdeg', k, gamma );
        [clusters, ncut_thres, lambda_thres] = opt_thresh_vol_cnstr_ncut_subset_penalty(f_new_full, W, gdeg', k, [], gamma );
                                                    % we use the opt thres
                                                    % function of penalty
                                                    % version but send
                                                    % empty seed

        ncut = balanced_cut(W, sum(W,2), clusters);
        assert ( abs(ncut - ncut_thres) <= 1e-8 );

        if ~isempty(subset)
        % Check if the paritition (S, S^c) is optimal!
        % This corresponds to allowing empty solution for the new problem on remaining
        % graph.
            clusters_subset = zeros(num,1);
            clusters_subset(subset) = 1;
            ncut_subset = balanced_cut(W, sum(W,2), clusters_subset );

            if ncut > ncut_subset

                clusters = clusters_subset;
                ncut = ncut_subset;

            end
        end
    else
        %[clusters, ncut1, ncut2] = opt_thresh_vol_cnstr_ncut_subset(start, W, gdeg', 1, k, gamma );
        [clusters, ncut_thres, lambda_thres] = opt_thresh_vol_cnstr_ncut_subset_penalty(start, W, gdeg', k, [], gamma );
        ncut = balanced_cut(W, sum(W,2), clusters);
        assert ( abs(ncut - ncut_thres) <= 1e-8 );

    end

    lambda=vol_cnstr_ncut_subset_functional(clusters,gamma,num,k,deg,wval,ix,jx,gdeg, totVol);                  
    objective = lambda*totVol;

    if gamma==0
        assert( abs(ncut - objective) <= 1e-8 );
    else
        assert(  lambda*totVol - lambda_thres >= -1e-8 ); % the difference should always be non-negative, but because of numerical accuracies, we allow some eps diff!
        display('Optimal thresholding did not increase the objective');
    end
    
end
    

