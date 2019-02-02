function [clusters, ncut,best_objective, feasible] = ...
    opt_thresh_vol_cnstr_ncut(f,W, deg, gdeg,hdeg, k, gamma)
% Performs optimal thresholding of the unconstrained set ratio for ncut 
% where the volume constraint has been incorporated into the objective as 
% penalty term.
%
% Input:
% f                 The vector (dimension of full graph).
% W                 Weight matrix (full graph including seed set).
% deg               The degree vector (of the full graph).
% gdeg              Generalized degree function for the denominator.
% hdeg              Generalized degree function for the constraint.
% k                 Upper bound on generalized volume (with respect to hdeg).
% gamma             Penalty parameter for volume constraint.
%
% Output:
% clusters          Thresholded vector f yielding the best objective.
% ncut              (Generalized) Ncut value of resulting clustering.
% best_objective    Corresponding objective.
% feasible          True if constraint is fulfilled.
%
% If the result is feasible (volume constraint is fulfilled), 
% the objective value agrees with the ncut objective.

    % make sure this is used correctly
    assert(size(f,1)==size(W,1),'Input vector has wrong dimension.');
    assert(size(deg,1)==length(f),'Degree vector has wrong dimension.');
    assert(size(gdeg,1)==length(f),'Degree vector has wrong dimension.');
    assert(size(hdeg,1)==length(f),'Degree vector has wrong dimension.');   
    assert(~issparse(gdeg),'Vector should not be sparse.');
    assert(~issparse(hdeg),'Vector should not be sparse.');
    assert(~issparse(deg),'Vector should not be sparse.');
    
    % sort vector f
    [f_sorted, index_sorted]=sort(f);
    [f_unique, index_unique]=unique(f_sorted,'first');
    
    % sort everything else accordingly
    W_sorted=W(index_sorted,index_sorted);
    gdeg_sorted = gdeg(index_sorted);
    hdeg_sorted = hdeg(index_sorted);
    deg_sorted = deg(index_sorted);

    % calculate cuts
    % we have to use here the original degrees not generalized ones
    vol_thresh=cumsum(deg_sorted);
    triup=triu(W_sorted,1);
    temp=2*cumsum(full(sum(triup,1)))' + cumsum(full(diag(W_sorted)));
    cuts_thresh=[0;vol_thresh(1:end-1)] -  [0;temp(1:end-1)];

    % those are the terms appearing in the denominator
    gvol_compl_thresh = [0; cumsum(gdeg_sorted(1:end-1))];  % gvol(V-C_i)
    gvol_thresh = sum(gdeg)- gvol_compl_thresh;             % gvol(C_i)
   
    % now extract all the unique indices
    cuts_thresh=cuts_thresh(index_unique);
    gvol_compl_thresh = gvol_compl_thresh(index_unique);
    gvol_thresh = gvol_thresh(index_unique);

    % The first possible cluster now consists of n entries, the size of 
    % the rest is n. The last cluster consists of i entries, where i is the 
    % number of occurences of the given value. The size of the rest is i.
    % Thus we have to exclude the first entry. 
    cuts_thresh = cuts_thresh(2:end);
    f_unique=f_unique(2:end);

    % compute balancing term (S)
    balance_thresh= gvol_compl_thresh(2:end).*gvol_thresh(2:end);

    % now compute the value of all the penalty terms
    cum_hvols = sum(hdeg_sorted) - [0; cumsum(hdeg_sorted(1:end-1))];
    penalty_thresh = gamma * max(0, cum_hvols-k);
    penalty_thresh = penalty_thresh(index_unique);
    penalty_thresh = penalty_thresh(2:end);
        
    % now compute total objective
    objective_thresh = (cuts_thresh + penalty_thresh)./balance_thresh;
   
    % find best one and compute indicator vector
    [best_objective,thresh_index]=min(objective_thresh);

    % compute output
    clusters=double(f>=f_unique(thresh_index));
    feasible=penalty_thresh(thresh_index)==0;
    ncut=cuts_thresh(thresh_index)/balance_thresh(thresh_index);
    
    % just to be sure
    assert(abs(ncut-balanced_cut(W,gdeg,clusters)/sum(sum(W)))< 1E-12);
    if (feasible) assert(best_objective==ncut); end
    
end

