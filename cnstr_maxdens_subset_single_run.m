function [dc, maxdens, lambda] = ... 
    cnstr_maxdens_subset_single_run(W, k1, k2, gdeg, start, subset, gam, verbosity)
% Performs one run with initialization given by start of the RatioDCA for
% the constrained maximum density problem.
%
% Usage: [dc, maxdens, lambda_new] ... 
%    = cnstr_maxdens_subset_single_run(W, k1, k2, gdeg, start, subset, gam, verbosity)
%
% Input: 
% W         The weight matrix.
% k1        Lower bound
% k2        Upper bound
% gdeg      Generalized degrees
% subset    Index of seed set
% gam       Penalty parameter gamma.
% verbosity Controls how much information is displayed [0-3]. Default is 1.
%
% Output:
% dc        Indicator vector of the found community
% maxdens   The corresponding density
% lambda    The value of the objective

    % do some checks
    if k2>size(W,1)
        k2 = size(W,1);
        if (verbosity>0) fprintf('Setting k2 to %d.\n', k2); end
    end
    assert(k2>length(subset),'Error: k2 has to be larger than the size of the subset');
    assert(k1>=1,'Error: k1 has to be at least 1');
    assert(k1<=k2,'Error:k1 has to be less or equal to k2');

    if (verbosity>1)
        fprintf("... Solving maximum density problem for gamma=%.5g:\n", gam);
    end

    % initialization
    start = abs(start);
    start = start/norm(start,2);
    f = start;

    inner_obj_all = 0;

    num = length(f);
    deg = sum(W,2);
   
    % restrict graph
    ind_rest = setdiff((1:num)', subset);
    W_rest = W(ind_rest,ind_rest);
    [ix_rest, jx_rest, wval_rest] = find(W_rest);

    degJ = sum(W(ind_rest,subset), 2);
    deg_tilde = deg(ind_rest)+degJ;

    assocJ = sum(sum(W(subset, subset)));
    gvolJ = sum(gdeg(subset));

    f = f(ind_rest);
    gdeg_rest = gdeg(ind_rest);
    
    assert(nnz(W_rest)==length(wval_rest));
    
    k1prime = max(k1-length(subset), 1);
    k2prime = k2-length(subset);

    assert(k1prime>=1);
    assert(k2prime>=1);

    % if no gamma is given, we compute it using a subset of V with
    % cardinality k2
    if nargin<7
        ixx = ix_rest(1);
        if ixx-1+k2prime<=length(ind_rest)
            indgam = ixx:1:ixx-1 +k2prime;
        else
            indgam = num-k2prime+1:num;
        end
        indtemp = [subset; ind_rest(indgam)];
        gam = sum(deg)*sum(gdeg(indtemp)) / sum(sum(W(indtemp, indtemp)));
    end

    % compute the value of the objective
    [lambda, indvec, indvec1] = compute_lambda(f, gam, k1, k2, deg(ind_rest), ...
                   gvolJ, degJ, assocJ, wval_rest, ix_rest, jx_rest, gdeg_rest);
    
    lambda_all = lambda;
    
    % make some output
    if (verbosity>1)
        fprintf('... gamma=%.5g \t lambda=%.5g \n', full(gam), lambda);
    end
    converged = false;
    it = 0;
    eps1 = 1E-4;
    while (~converged)
        it = it+1;

        % solve inner problem
        [f_new, lambda_new, indvec_new, indvec_new1, obj] = ... 
            solveInnerProblem(f, indvec1, indvec, gam, deg, W_rest, ...
                wval_rest, ix_rest, jx_rest, gdeg_rest, ind_rest, lambda, ...
                deg_tilde, assocJ, k1prime, k2prime, gvolJ, degJ, verbosity>2);
               
        % check if converged
        reldiff = abs(lambda_new-lambda)/lambda;
        converged = (reldiff < eps1);

        % make some output
        if (verbosity>1)
            fprintf('... it=%d\tlambda=%.5g\tdiff=%6.4g\tinnerobj=%.5g\n', ...
                    it, lambda_new, reldiff, obj);
        end     
      
        % update variables
        lambda = lambda_new;
        f = f_new;
        indvec = indvec_new;
        indvec1 = indvec_new1;

        inner_obj_all(it+1) = obj;
        lambda_all(it+1) = lambda_new;
    end

    % make some output
    if (verbosity>1)
        fprintf('  \n');
    end

    % perform optimal thresholding
    dc_temp = opt_thresh_cnstr_maxdens_subset(W_rest, f_new, ...
                  gdeg(ind_rest), max(1, k1prime), k2prime, gam, ...
                  assocJ, sum(gdeg(subset)), W(ind_rest, subset));
    dc = ones(size(W,1),1);
    dc(ind_rest) = dc_temp>0;
    maxdens = full(sum(sum(W(dc==1, dc==1))) / sum(gdeg(dc>0)));
end



% solves the convex inner problem in RatioDCA
function [f_new, lambda_new, indvec_new, indvec_new1, obj] = ... 
    solveInnerProblem(f, indvec1, indvec, gam, deg, W_rest, wval_rest, ...
         ix_rest, jx_rest, gdeg_rest, ind_rest, lambda, deg_tilde, ...
         assocJ, k1prime, k2prime, gvolJ, degJ, debug)

    % some parameters
    MAXITER = 1600;%800;

    % compute constants
    indmax = zeros(size(f,1),1);
    indmax(f==max(f)) = 1 / sum(f==max(f));
    c2 = gdeg_rest - gam*indvec1 + gam*indvec - lambda*deg_tilde;
    c1 = k1prime*gam + gvolJ - lambda*assocJ;
   
    % Lipschitz constant
    L = 2*lambda^2*max(sum(W_rest.^2));
    num = nnz(W_rest);
     
    % sanity check: compute primal obj using old f. should be close to 0
    obj_old = lambda * 0.5 * sum(wval_rest.*abs(f(ix_rest) - f(jx_rest))) + c1 * max(f) + sum(c2.*f);
    assert(abs(obj_old) < 1E-8);

    % solve inner problem
    [f_new, obj] = mex_IP_maxdens_FISTA(W_rest, c2, zeros(num,1), ...
                                        MAXITER, L, lambda, c1, debug);
    assert(obj<=0);
   
    % recompute dual objective (for testing)    
    if (debug)
        dualval2 = -norm(f_new)^2;
    else 
        dualval2 = inf;
    end
    
    % renormalize
    f_new = f_new/norm(f_new);

    % compute objective
    if (abs(obj)<=1E-15)
        f_new = f;
        lambda_new = lambda;
        indvec_new = indvec;
        indvec_new1 = indvec1;
    else
        % sanity check: compute primal obj
        obj2 = lambda*0.5*sum(wval_rest.*abs(f_new(ix_rest)-f_new(jx_rest))) + ...
               c1*max(f_new) + f_new'*c2;
        assert(abs(obj-obj2) < 1E-10 * max(abs(obj),1));

        [lambda_new, indvec_new, indvec_new1] = ... 
            compute_lambda(f_new, gam, k1prime, k2prime, deg(ind_rest), ...
                gvolJ, degJ, assocJ, wval_rest, ix_rest, jx_rest, gdeg_rest);
    end
end
