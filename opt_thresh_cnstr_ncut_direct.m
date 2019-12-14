function [clusters, ncut, lambda, feasible] = opt_thresh_cnstr_ncut_direct(...
                       f_rest, W_rest, g_rest, h_rest, kprime, gam, gvolJ, degJ)
% Performs optimal thresholding of the unconstrained set ratio for ncut 
% where the subset constraint has been directly incorporated into the 
% objective (leading to a problem of lower dimension) and the volume 
% constraint has been incorporated into the objective as penalty term.
% 
% Usage: [clusters, ncut, lambda, feasible] = opt_thresh_cnstr_ncut_direct(...
%                      f_rest, W_rest, g_rest, h_rest, kprime, gam, gvolJ, degJ)
%
% Input:
% f_rest            The vector (dimension of restricted graph).
% W_rest            Weight matrix (graph without seed set).
% g_rest            Generalized degree function for the denominator 
%                   (dimension of restricted graph).
% h_rest            Generalized degree function for the second constraint 
%                   (dimension of restricted graph).
% kprime            Upper bound on generalized volume (with respect to h, 
%                   adjusted for seed set).
% gam               Penalty parameter for volume constraint.
% gvolJ             Total generalized volume (gdeg) of subset.
% degJ              For each vertex, the sum of edge weights to the subset.
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
    assert(size(f_rest,1)==size(W_rest,1),'Input vector has wrong dimension.');
    assert(size(g_rest,1)==length(f_rest),'Degree vector has wrong dimension.');
    assert(size(h_rest,1)==length(f_rest),'Degree vector has wrong dimension.');   
    assert(~issparse(g_rest),'Vector should not be sparse.');
    assert(~issparse(h_rest),'Vector should not be sparse.');
    assert(~issparse(degJ),'Vector should not be sparse.');

    %% sort vector f, and everything else accordingly
    [f_sort, ind_sort] = sort(f_rest);
    W_sort = W_rest(ind_sort, ind_sort);
    g_sort = g_rest(ind_sort);
    h_sort = h_rest(ind_sort);
    degJ_sort = degJ(ind_sort);

    %% compute values of cut, gvol and hvol on the restricted graph
    % calculate first the cuts to the rest of the restricted graph
    % note that here we also include the cases where C'_i or V'-C'_i
    % are empty, thus we have one additional threshold
    triup = triu(W_sort,1);
    temp = 2*cumsum(full(sum(triup,1)))' + cumsum(full(diag(W_sort)));
    cuts_thresh = [0; cumsum(sum(W_sort,2))] - [0; temp];

    % compute the gvols for each threshold
    gvol_compl_thresh = [0; cumsum(g_sort)];        % gvol(V'-C'_i)
    gvol_thresh = sum(g_sort) - gvol_compl_thresh;  % gvol(C'_i)
    
    % compute the hvols for each threshold (hvol(C'_i))
    hvol_thresh = sum(h_sort) - [0; cumsum(h_sort)];
        
    %% compute values of cut, gvol and hvol incorporating terms for subset
    % we have to consider two cases: The subset can either be added to 
    % the first cluster (C'_i) or to the second cluster (V'-C'_i). 
    
    % Case 1: C_i = C'_i + J, V-C_i = V'-C'_i
    % In this case the total cut is cut(C'_i, V'-C'_i) + cut(J,V'-C'_i).
    % One has gvol(V-C_i) = gvol(V'-C_i), and gvol(C_i) = gvol(C'_i)+gvol(J).
    % The hvol constraint is given by hvol(C'_i) <=kprime (note that kprime 
    % is already adjusted for the restricted graph).
    cuts_thresh1 = cuts_thresh + [0;cumsum(degJ_sort)];
    gvol_compl_thresh1 = gvol_compl_thresh;
    gvol_thresh1 = gvol_thresh + gvolJ;
    hvol_thresh1 = hvol_thresh;
        
    % Case 2: C_i = C'_i, V-C_i = V'-C'_i + J
    % In this case the total cut is cut(C'_i, V'-C'_i) + cut(C'_i,J).
    % Now gvol(V-C_i) = gvol(V'-C_i)+gvol(J), and gvol(C_i) = gvol(C'_i).
    % The hvol constraint is given by hvol(V'-C'_i) <=kprime (note that 
    % kprime is already adjusted for the restricted graph).
    cuts_thresh2 = cuts_thresh + sum(degJ_sort) - [0;cumsum(degJ_sort)];
    gvol_compl_thresh2 = gvol_compl_thresh + gvolJ;
    gvol_thresh2 = gvol_thresh;
    hvol_thresh2 = sum(h_sort) - hvol_thresh;

    %% compute best threshold
    % Extract all possible thresholds. In order to also allow 
    % the empty set, we need to add an additional n+1th threshold.
    [f_unique, ind_unique] = unique([f_sort; max(f_sort)+1], 'first');
    
    % Compute best threshold for case 1
    [lambda1, ind_thresh1, feasible1, ncut1] = find_best_thresh(cuts_thresh1,...
       gvol_compl_thresh1, gvol_thresh1, hvol_thresh1, ind_unique, kprime, gam);
    
    % Compute best threshold for case 2
    [lambda2, ind_thresh2, feasible2, ncut2] = find_best_thresh(cuts_thresh2,...
       gvol_compl_thresh2, gvol_thresh2, hvol_thresh2, ind_unique, kprime, gam);
        
    % check which is better
    if (lambda1<lambda2)
        lambda = lambda1;
        clusters = double(f_rest>=f_unique(ind_thresh1));
        feasible = feasible1;
        ncut = ncut1;
    else
        lambda = lambda2;
        clusters = double(f_rest<f_unique(ind_thresh2));
        feasible = feasible2;
        ncut = ncut2;
    end    
end


% compute the best threshold according to the ncut objective with subset 
% incorporated and hvol constraint as penalty
function [obj, ind_thresh, feasible, ncut] = find_best_thresh(cuts_thresh, ...
         gvol_compl_thresh, gvol_thresh, hvol_thresh, ind_unique, kprime, gam)

    % extract the values for thresholding at all the unique indices
    cuts_thresh = cuts_thresh(ind_unique);
    gvol_compl_thresh = gvol_compl_thresh(ind_unique);
    gvol_thresh = gvol_thresh(ind_unique);
    hvol_thresh = hvol_thresh(ind_unique);

    % compute objective for different thresholds
    balance_thresh = gvol_compl_thresh.*gvol_thresh;
    penalty_thresh = gam * max(0, hvol_thresh-kprime);
    lambda_thresh = (cuts_thresh + penalty_thresh)./balance_thresh;
   
    % find best one (this will automatically exclude the case where the 
    % set is empty since it cannot the minimum)
    [obj, ind_thresh] = min(lambda_thresh);
    
    % compute output
    feasible = hvol_thresh(ind_thresh)<=kprime;
    ncut = cuts_thresh(ind_thresh) / balance_thresh(ind_thresh);    
end
