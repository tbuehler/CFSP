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
    [~, ind_sort] = sort(f);

    % sort everything else accordingly
    W_sort = W(ind_sort, ind_sort);
    g_sort = g(ind_sort);
    h_sort = h(ind_sort);
    deg_sort = deg(ind_sort);

    % calculate cuts
    triup = triu(W_sort,1);
    temp = 2*cumsum(full(sum(triup,1)))' + cumsum(full(diag(W_sort)));
    cuts_thresh = [0; cumsum(deg_sort)] - [0; temp];

    % compute balancing term (S)
    gvol_compl_thresh = [0; cumsum(g_sort)];     % gvol(V-C_i)
    gvol_thresh = sum(g) - gvol_compl_thresh;    % gvol(C_i)
    balance_thresh = gvol_compl_thresh.*gvol_thresh;

    % now compute the value of the penalty term for the seed constraint
    subset_vec = zeros(length(f),1); 
    subset_vec(subset) = 1;
    subset_vec = subset_vec(ind_sort);
    pen_subset_thresh = gam1 * [0; cumsum(subset_vec)]; 
    
    % now compute the value of the penalty term for the volume constraint
    cum_hvols = sum(h_sort) - [0; cumsum(h_sort)];
    pen_vol_thresh = gam2 * max(0, cum_hvols-k);
        
    % now compute total objective
    lambda_thresh = (cuts_thresh + pen_vol_thresh ...
                    + pen_subset_thresh)./balance_thresh;
   
    % find best one (explicitly excluding case of empty / full cluster)
    [lambda, ind_thresh] = min(lambda_thresh(2:end-1));
    ind_thresh = ind_thresh + 1; 

    % compute output
    clusters = zeros(size(f,1), 1);
    clusters(ind_sort(ind_thresh:end)) = 1;
    feasible = sum(h(clusters==1))<=k && sum(clusters(subset))==length(subset);
    ncut = cuts_thresh(ind_thresh) / balance_thresh(ind_thresh);
    
    % just to be sure
    assert(abs(ncut-balanced_cut(W,g,clusters)/sum(sum(W)))< 1E-12);
    if (feasible); assert(lambda==ncut); end
end
