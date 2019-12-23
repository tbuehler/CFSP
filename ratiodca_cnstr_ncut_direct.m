function [clusters, ncut, feasible, lambda] = ratiodca_cnstr_ncut_direct(W, ...
                                            k, h, start, subset, gam, verbosity)
% Solves the normalized cut problem with generalized volume constraint and
% subset constraint using RatioDCA. The subset constraint has been directly 
% incorporated into the objective (leading to a problem of lower dimension) 
% and the volume constraint has been incorporated into the objective as 
% penalty term.
%
% Corresponding paper:
% T. Buehler, S. S. Rangapuram, S. Setzer and M. Hein
% Constrained fractional set programs and their application in local clustering 
% and community detection
% ICML 2013, pages 624-632 (Extended version: http://arxiv.org/abs/1306.3409)
% 
%
% Usage: [clusters, ncut, feasible, lambda] = ratiodca_cnstr_ncut_direct(W, ...
%                                           k, h, start, subset, gam, verbosity)
%
% Input:
% W             Weight matrix (full graph).
% k             Upper bound.
% h             Generalized degree vector used in constraint.
% start         Start vector.
% subset        Indices of seed subset.
% gam           Penalty parameter for volume constraint.
% verbosity     Controls how much information is displayed [0-3]. Default is 1.
%
% Output:
% clusters      Thresholded vector f yielding the best objective.
% ncut          Ncut value of resulting clustering.
% feasible      True if all constraints are fulfilled.
% lambda        Corresponding objective.
%
%
% (C)2012-19 Thomas Buehler, Syama Rangapuram, Simon Setzer and Matthias Hein

    %%  check inputs
    if k>sum(h)
        k=sum(h);
        if (verbosity>1); fprintf('... Setting k to %f\n',k); end
    end
    assert(k>=sum(h(subset)),'Error! Problem is unfeasible.');
    assert(gam>=0,'Error! gam cannot be negative.');

    if (verbosity>1) 
        fprintf("... Solving constrained ncut problem for gam=%.5g:\n", gam);
    end

    %% initialization
    start = abs(start);
    start = start/norm(start,2);
    f = start;
    
    %% some helpers 
    num = length(f);
    deg = full(sum(W,2));
    
    %% restrict graph
    ind_rest = setdiff((1:num)',subset);
    W_rest = W(ind_rest,ind_rest);
    [ix_rest, jx_rest, wval_rest] = find(W_rest);
    f = f(ind_rest);
    kprime = k-sum(h(subset));
    degJ = full(sum(W(ind_rest,subset),2));
    totVol_rest = full(sum(sum(W(ind_rest, ind_rest))));
    Deg = spdiags(deg(ind_rest),0,size(f,1),size(f,1));
    W_triu = triu(W_rest,1); % diagonal entries play no role in inner problem
    L = 4*max(sum(W_rest.^2)); % Lipschitz constant used in inner problem    

    %% evaluate objective
    [lambda, indvec] = lambda_cnstr_ncut_direct(f, gam, num-length(subset), k, ...
       deg(ind_rest), wval_rest, ix_rest, jx_rest, h(ind_rest), totVol_rest, degJ);

    it = 0;
    converged = false;
    eps1 = 1E-4;

    if (verbosity>1) 
        fprintf('... it=%d\tlambda=%.5g\n',it, lambda);
    end

    %% main loop
    while (~converged)
        it = it+1;
 
        % solve inner problem
        [f_new, lambda_new, indvec_new, obj] = solveInnerProblem(f, lambda, ...
           subset, k, W_triu, wval_rest, ix_rest, jx_rest, num, deg, Deg, ...
           h, ind_rest, gam, indvec, degJ, totVol_rest, L, verbosity>2);
        
        % check if converged
        reldiff = abs(lambda_new-lambda)/lambda;
        converged = (reldiff < eps1);       

        if (verbosity>1) 
            fprintf('... it=%d\tlambda=%10.4g\tdiff=%6.4g\tinnerobj=%10.4g\n', ...
                it, lambda_new, reldiff, obj);
        end

        % update variables
        lambda = lambda_new;
        f = f_new;
        indvec = indvec_new;
    end

    %% Perform optimal thresholding
    [clusters_temp, ncut_prime, lambda, feasible] = opt_thresh_cnstr_ncut_direct(...
        f_new, W_rest, deg(ind_rest), h(ind_rest), kprime, gam, sum(deg(subset)), degJ);
    ncut = ncut_prime * full(sum(sum(W)));    

    if (verbosity>1) 
        fprintf('... Result after optimal thresholding: lambda=%.5g ncut''=%.5g ncut=%.5g feasible=%d\n', ... 
                lambda, ncut_prime, ncut, feasible);
    end

    %% construct vector with respect to original graph
    clusters = zeros(num,1);
    clusters(subset) = 1;
    clusters(ind_rest) = clusters_temp;

    % sanity checks
    assert ( abs(balanced_cut(W, sum(W,2), clusters) - ncut) <= 1e-8 );
    assert ( feasible == (sum(h(clusters==1))<=k && sum(clusters(subset))==length(subset)) );

    %% check if partition (S=seed set, S^c) is optimal
    if ~isempty(subset)
        clusters_subset = zeros(num,1);
        clusters_subset(subset) = 1;
        ncut_subset = balanced_cut(W, sum(W,2), clusters_subset);
        feasible_subset = k>=sum(h(subset)); 
        if (verbosity>1) 
            fprintf('... Edge case of empty set: ncut=%.5g feasible=%d\n', ...
                    ncut_subset, feasible_subset);
        end
        if feasible_subset && ncut > ncut_subset
            clusters = clusters_subset;
            ncut = ncut_subset;
            feasible = feasible_subset;
            lambda = ncut_subset;
        end
    end
    if (verbosity>1) 
        fprintf('... Final result: lambda=%.5g ncut=%.5g feasible=%d\n', ...
                lambda, ncut, feasible);
    end
end



%% solves the inner problem in RatioDCA
function [f_new, lambda_new, indvec_new, obj] = solveInnerProblem(f, lambda, ...
         subset, k, W_triu, wval_rest, ix_rest, jx_rest, num, deg, Deg, h, ...
         ind_rest, gam, indvec, degJ, totVol_rest, L, debug)
        
    % set parameters
    MAXITER = 5120;
    MAXITER_start = 40;

    % compute subgradient of subset penalty
    indmax = zeros(size(f,1),1);
    indmax(f==max(f)) = 1 / sum(f==max(f));
    
    % compute subgradient of ncut balancing term
    Pf = f - (f'*deg(ind_rest)/sum(deg(ind_rest)));
    sg_bal = Deg*sign(Pf) - deg(ind_rest)*sum( Deg*sign( Pf))/sum(deg(ind_rest));
    
    % compute constants
    c1 = sum(degJ);
    c2 = gam*indvec - degJ - (lambda/2)*totVol_rest*sg_bal ... 
        - sum(deg(subset))*totVol_rest*lambda*indmax ... 
        + lambda* sum(deg(subset)) *deg(ind_rest);
    
    % sanity check: compute primal obj using old f. should be close to 0
    obj_old = 0.5 * sum(wval_rest.*abs(f(ix_rest) - f(jx_rest))) + c1 * max(f) + sum(c2.*f);
    assert(abs(obj_old) < 1E-8);

    % solve inner problem with FISTA
    if (c1~=0)
        [f_new, obj] = mex_ip_cnstr_ncut_subset(W_triu, c2, zeros(nnz(W_triu),1), ...
                       MAXITER, 1E-8, L, c1, MAXITER_start, debug);
    else
        [f_new, obj] = mex_ip_cnstr_ncut(W_triu, c2, zeros(nnz(W_triu),1), ...
                       MAXITER, 1E-8, L, MAXITER_start, debug);
    end
    assert(obj<=0);    

    % in this case, take the last result
    if abs(obj)<=1E-15
        f_new = f;
        lambda_new = lambda;
        indvec_new = indvec;
    else
        % sanity check: recompute primal obj
        f_new = f_new/norm(f_new);
        obj2 = 0.5 * sum(wval_rest.*abs(f_new(ix_rest) - f_new(jx_rest))) + c1 * max(f_new) + sum(c2.*f_new);
        assert(abs(obj-obj2) < 1E-10 * max(abs(obj),1));
    
        [lambda_new, indvec_new] = lambda_cnstr_ncut_direct(f_new, gam, num-length(subset), ...
            k, deg(ind_rest), wval_rest, ix_rest, jx_rest, h(ind_rest), totVol_rest, degJ);
    end
    assert(lambda_new <= lambda + 1E-15);
end

