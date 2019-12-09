function [clusters, ncut, lambda, feasible] = opt_thresh_cnstr_ncut_penalty(...
                                         f, W, deg, g, h, k, subset, gam1, gam2)
% Performs optimal thresholding of the unconstrained set ratio for ncut where
% the volume constraint and subset constraint have been incorporated into the 
% objective as penalty terms.
% 
% Usage: [clusters, ncut, lambda, feasible] = opt_thresh_cnstr_ncut_penalty(...
%                                        f, W, deg, g, h, k, subset, gam1, gam2)
%
% Input:
% f                 The vector (dimension of full graph).
% W                 Weight matrix (full graph including seed set).
% deg               The degree vector (of the full graph).
% g                 Generalized degree function for the denominator.
% h                 Generalized degree function for the second constraint.
% k                 Upper bound on generalized volume (with respect to h).
% subset            The indices of the seed subset.
% gam1              Penalty parameter for subset constraint.
% gam2              Penalty parameter for volume constraint.
%
% Output:
% clusters          Thresholded vector f yielding the best objective.
% ncut              (Generalized) Ncut value of resulting clustering.
% lambda            Corresponding objective.
% feasible          True if all constraints are fulfilled.
%
% If the result is feasible (volume and seed constraints are fulfilled), 
% the objective value agrees with the ncut objective.
%
% (C)2012-19 Thomas Buehler, Syama Rangapuram, Simon Setzer and Matthias Hein
    
    %% make sure this is used correctly
    assert(size(f,1)==size(W,1),'Input vector has wrong dimension.');
    assert(size(deg,1)==length(f),'Degree vector has wrong dimension.');
    assert(size(g,1)==length(f),'Degree vector has wrong dimension.');
    assert(size(h,1)==length(f),'Degree vector has wrong dimension.');   
    assert(~issparse(g),'Vector should not be sparse.');
    assert(~issparse(h),'Vector should not be sparse.');
    assert(~issparse(deg),'Vector should not be sparse.');
    
    % sort vector f
    [f_sort, ind_sort] = sort(f);
    [f_unique, ind_unique] = unique(f_sort, 'first');

    % sort everything else accordingly
    W_sort = W(ind_sort, ind_sort);
    g_sort = g(ind_sort);
    h_sort = h(ind_sort);
    deg_sort = deg(ind_sort);

    % calculate cuts
    vol_thresh = cumsum(deg_sort);
    triup = triu(W_sort,1);
    temp = 2*cumsum(full(sum(triup,1)))' + cumsum(full(diag(W_sort)));
    cuts_thresh = [0;vol_thresh(1:end-1)] - [0;temp(1:end-1)];

    % those are the terms appearing in the denominator
    gvol_compl_thresh = [0; cumsum(g_sort(1:end-1))];     % gvol(V-C_i)
    gvol_thresh = sum(g) - gvol_compl_thresh;             % gvol(C_i)
   
    % now extract all the unique indices
    cuts_thresh = cuts_thresh(ind_unique);
    gvol_compl_thresh = gvol_compl_thresh(ind_unique);
    gvol_thresh = gvol_thresh(ind_unique);

    % The first possible cluster now consists of n entries, the size of 
    % the rest is n. The last cluster consists of i entries, where i is the 
    % number of occurences of the given value. The size of the rest is i.
    % Thus we have to exclude the first entry. 
    cuts_thresh = cuts_thresh(2:end);
    f_unique = f_unique(2:end);

    % compute balancing term (S)
    balance_thresh = gvol_compl_thresh(2:end).*gvol_thresh(2:end);

    % now compute the value of the penalty term for the seed constraint
    subset_vec = zeros(length(f),1); 
    subset_vec(subset) = 1;
    subset_vec = subset_vec(ind_sort);
    pen_subset_thresh = [0;cumsum(subset_vec(1:end-1))]; 
    pen_subset_thresh = pen_subset_thresh(ind_unique);
    pen_subset_thresh = gam1*pen_subset_thresh(2:end);
    
    % now compute the value of the penalty term for the volume constraint
    cum_hvols = sum(h_sort) - [0; cumsum(h_sort(1:end-1))];
    pen_vol_thresh = gam2 * max(0, cum_hvols-k);
    pen_vol_thresh = pen_vol_thresh(ind_unique);
    pen_vol_thresh = pen_vol_thresh(2:end);
        
    % now compute total objective
    lambda_thresh = (cuts_thresh + pen_vol_thresh ...
                    + pen_subset_thresh)./balance_thresh;
   
    % find best one and compute indicator vector
    [lambda, ind_thresh] = min(lambda_thresh);

    % compute output
    clusters = double(f>=f_unique(ind_thresh));
    feasible = sum(h(clusters==1))<=k && sum(clusters(subset))==length(subset);
    ncut = cuts_thresh(ind_thresh) / balance_thresh(ind_thresh);
    
    % just to be sure
    assert(abs(ncut-balanced_cut(W,g,clusters)/sum(sum(W)))< 1E-12);
    if (feasible); assert(lambda==ncut); end
end
