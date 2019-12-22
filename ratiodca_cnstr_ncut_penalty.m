function [clusters, ncut, feasible, lambda] = ratiodca_cnstr_ncut_penalty(W, ...
                                    k, h, start, subset, gam1, gam2, verbosity)
% Solves the normalized cut problem with generalized volume constraint and
% subset constraint using RatioDCA. Both subset constraint and volume constraint
% have been incorporated into the objective as penalty term.
%
% Corresponding paper:
% T. Buehler, S. S. Rangapuram, S. Setzer and M. Hein
% Constrained fractional set programs and their application in local clustering 
% and community detection
% ICML 2013, pages 624-632 (Extended version: http://arxiv.org/abs/1306.3409)
%
%
% Usage: [clusters, ncut, feasible, lambda] = ratiodca_cnstr_ncut_penalty(W, ...
%                                    k, h, start, subset, gam1, gam2, verbosity)
%
% Input:
% W             Weight matrix (full graph).
% k             Upper bound.
% h             Generalized degree vector used in constraint.
% start         Start vector.
% subset        Indices of seed subset.
% gam1          Penalty parameter for seed constraint.
% gam2          Penalty parameter for volume constraint.
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

    %% check inputs
    assert(k <= sum(h), 'Error! Upper bound exceeds total volume.');
    assert(k >= sum(h(subset)),  ... 
        'Error! Upper bound is smaller than volume of the seed subset.');
    assert(gam1>=0,'Error! gam1 cannot be negative.');
    assert(gam2>=0,'Error! gam2 cannot be negative.');
    assert(~issparse(h),'Error! h should not be sparse.');
    
    %% initialization
    start = abs(start);
    start = start/norm(start,2);
    f = start;
    
    %% some helpers 
    num = length(f);
    deg = full(sum(W,2));
    totVol = full(sum(sum(W)));
    [ix,jx,wval] = find(W);
    Deg = spdiags(deg, 0, size(f,1), size(f,1));
    W_triu = triu(W,1);
    L = 4*max(sum(W.^2));

    %% evaluate objective
    f(subset) = max(f);
    [lambda, indvec] = lambda_cnstr_ncut_penalty(f, gam1, gam2, num, k, ...
                       subset, deg, wval, ix, jx, h, totVol);

    it = 0;
    converged = false;
    eps1 = 1E-4;
    if (verbosity>1)     
        fprintf('... it=%d\tlambda=%.5g\n',it, lambda);
    end
    
    %% main loop
    while (~converged && max(abs(f)) ~= 0) % avoid zero starting point (on the subgraph).
        it = it+1;
        
        % solve inner problem
        [f_new, lambda_new, indvec_new, obj] = solveInnerProblem(f, lambda, ...
           subset, k, W_triu, wval, ix, jx, num, deg, Deg, h, gam1, gam2, indvec, totVol, L, verbosity>2);
                
        % check if converged
        reldiff = abs(lambda_new-lambda)/lambda;
        converged = (reldiff < eps1 || lambda_new==0 || lambda==0);       

        if (verbosity>1) 
            fprintf('... it=%d\tlambda=%10.4g\tdiff=%6.4g\tinnerobj=%10.4g\n', ...
                it, lambda_new, reldiff, obj);
        end

        % update variables
        lambda = lambda_new;
        f = f_new;
        indvec = indvec_new;
    end
    
    %% perform optimal thresholding
    if (max(abs(f))==0) % in this case we take the start vector
        f = start;
    end
    [clusters, ncut_prime, lambda_thresh, feasible] = opt_thresh_cnstr_ncut_penalty(...
                                    f, W, deg, deg, h, k, subset, gam1, gam2);
    ncut = ncut_prime * full(sum(sum(W)));

    if (verbosity>1) 
        fprintf('... Result after optimal thresholding: lambda=%.5g ncut''=%.5g ncut=%.5g feasible=%d\n', ... 
                lambda_thresh, ncut_prime, ncut, feasible);
    end

    % sanity checks
    assert ( abs(balanced_cut(W, sum(W,2), clusters) - ncut) <= 1e-8 );
    assert ( feasible == (sum(h(clusters==1))<=k && sum(clusters(subset))==length(subset)) );
    
    %% compute objective
    lambda = lambda_cnstr_ncut_penalty(clusters, gam1, gam2, num, k, ...
                                       subset, deg, wval, ix, jx, h, totVol);
    assert (lambda_thresh - lambda <= 1E-8);                  
  
    if gam1==0 && gam2==0
        assert( abs(ncut/sum(sum(W)) - lambda) <= 1e-8 );
    end
    if (verbosity>1) 
        fprintf('... Final result:\tlambda=%.5g ncut=%.5g feasible=%d\n', ...
                lambda, ncut, feasible);
    end
end



%% solves the inner problem in RatioDCA
function [f_new, lambda_new, indvec_new, obj] = solveInnerProblem(f, lambda, ...
         subset, k, W_triu, wval, ix, jx, num, deg, Deg, h, gam1, gam2, indvec, totVol, L, debug)
       
    % set parameters
    MAXITER = 5120;
    MAXITER_start = 40;

    % compute subgradient of ncut balancing term
    Pf = f - (f'*deg/sum(deg));
    sg_bal = Deg*sign(Pf) - deg*sum(Deg*sign(Pf))/sum(deg);

    % subgradient for subset constraint
    gam_subset = zeros(num,1);
    gam_subset(subset) = gam1;
    
    % compute constants
    c1 = gam1*length(subset);
    c2 = gam2*indvec - gam_subset - (lambda/2)*totVol*sg_bal;
    
    % sanity check: compute primal obj using old f. should be close to 0
    obj_old = 0.5 * sum(wval.*abs(f(ix) - f(jx))) + c1 * max(f) + sum(c2.*f);
    assert(abs(obj_old) < 1E-10);

    % solve inner problem with FISTA
    if gam1 ~= 0
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
        obj2 = 0.5 * sum(wval.*abs(f_new(ix) - f_new(jx))) + c1 * max(f_new) + sum(c2.*f_new);
        assert(abs(obj-obj2) < 1E-10 * max(abs(obj),1));

        [lambda_new, indvec_new] = lambda_cnstr_ncut_penalty(f_new, gam1, gam2, ...
                                   num, k, subset, deg, wval, ix, jx, h, totVol);
    end
    assert(lambda_new<=lambda +1E-15);
end
