function [clusters, ncut,best_objective, feasible] = ...
    opt_thresh_vol_cnstr_ncut_subset_direct(f_rest,W_rest, gdeg_rest,hdeg_rest,kprime, gamma,gvolJ,degJ)
% Performs optimal thresholding of the unconstrained set ratio for ncut 
% where the subset constraint has been directly incorporated into the 
% objective (leading to a problem of lower dimension) and the volume 
% constraint has been incorporated into the objective as penalty term.
% 
% Input:
% f_rest            The vector (dimension of restricted graph).
% W_rest            Weight matrix (graph without seed set).
% gdeg_rest         Generalized degree function for the denominator 
%                   (dimension of restricted graph).
% hdeg_rest         Generalized degree function for the second constraint 
%                   (dimension of restricted graph).
% kprime            Upper bound on generalized volume (with respect to 
%                   hdeg, adjusted for seed set).
% gamma             Penalty parameter for volume constraint.
% gvolJ             Total generalized volume (gdeg) of subset.
% degJ              For each vertex, the sum of edge weights to the subset.
%
% Output:
% clusters          Thresholded vector f yielding the best objective.
% ncut              (Generalized) Ncut value of resulting clustering.
% best_objective    Corresponding objective.
% feasible          True if all constraints are fulfilled.
%
% If the result is feasible (volume and seed constraints are fulfilled), 
% the objective value agrees with the ncut objective.


    %% make sure this is used correctly
    assert(size(f_rest,1)==size(W_rest,1),'Input vector has wrong dimension.');
    assert(size(gdeg_rest,1)==length(f_rest),'Degree vector has wrong dimension.');
    assert(size(hdeg_rest,1)==length(f_rest),'Degree vector has wrong dimension.');   
    assert(~issparse(gdeg_rest),'Vector should not be sparse.');
    assert(~issparse(hdeg_rest),'Vector should not be sparse.');
    assert(~issparse(degJ),'Vector should not be sparse.');
        
    
    %% sort vector f, and everything else accordingly
    [f_sorted, index_sorted]=sort(f_rest);
    W_sorted=W_rest(index_sorted,index_sorted);
    gdeg_sorted = gdeg_rest(index_sorted);
    hdeg_sorted = hdeg_rest(index_sorted);
    degJ_sorted=degJ(index_sorted);
    
      
    %% compute values of cut, gvol and hvol on the restricted graph
    % calculate first the cuts to the rest of the restricted graph
    % we have to use here the original degrees not generalized ones
    % note that here we also include the cases where C'_i or V'-C'_i
    % are empty, thus we have one additional threshold
    deg_inside_rest=sum(W_sorted,2);
    vol_thresh=cumsum(deg_inside_rest);
    triup=triu(W_sorted,1);
    temp=2*cumsum(full(sum(triup,1)))' + cumsum(full(diag(W_sorted)));
    cuts_thresh=[0;vol_thresh] -  [0;temp];

    % compute the gvols for each threshold
    gvol_compl_thresh = [0; cumsum(gdeg_sorted)];       % gvol(V'-C'_i)
    gvol_thresh = sum(gdeg_sorted)- gvol_compl_thresh;  % gvol(C'_i)
    
    % compute the hvols for each threshold (hvol(C'_i))
    hvol_thresh = sum(hdeg_sorted) - [0; cumsum(hdeg_sorted)];
    
    
    %% compute values of cut, gvol and hvol incorporating terms for subset
    % we have to consider two cases now: The subset can either be added to 
    % the first cluster (C'_i) or to the second cluster (V'-C'_i). 
    % The second case corresponds to thresholding with < instead of >=, or 
    % equivalently, optimal thresholding of -f. 
    
    % Case 1: C_i = C'_i + J, V-C_i = V'-C'_i
    % In this case the total cut is cut(C'_i, V'-C'_i) + cut(J,V'-C'_i).
    % One has gvol(V-C_i)= gvol(V'-C_i), and gvol(C_i)= gvol(C'_i)+gvol(J).
    % The hvol constraint is given by hvol(C'_i) <=kprime (note that kprime 
    % is already adjusted for the restricted graph).
    cuts_thresh_case1=cuts_thresh + [0;cumsum(degJ_sorted)];
    gvol_compl_thresh_case1= gvol_compl_thresh;
    gvol_thresh_case1=gvol_thresh + gvolJ;
    penalty_thresh_case1 = gamma * max(0, hvol_thresh-kprime);
    
    % Case 2: C_i = C'_i, V-C_i = V'-C'_i + J
    % In this case the total cut is cut(C'_i, V'-C'_i) + cut(C'_i,J).
    % Now gvol(V-C_i) = gvol(V'-C_i)+gvol(J), and gvol(C_i) = gvol(C'_i).
    % The hvol constraint is given by hvol(V'-C'_i) <=kprime (note that 
    % kprime is already adjusted for the restricted graph).
    cuts_thresh_case2=cuts_thresh +  sum(degJ_sorted) - [0;cumsum(degJ_sorted)];
    gvol_compl_thresh_case2= gvol_compl_thresh + gvolJ;
    gvol_thresh_case2=gvol_thresh;
    penalty_thresh_case2 = gamma * max(0, sum(hdeg_sorted) - hvol_thresh-kprime);
        
   
    %% compute best threshold
    % Extract all possible thresholds. In order to also allow 
    % the empty set, we need to add an additional n+1th threshold.
    [f_unique, index_unique]=unique([f_sorted; max(f_sorted)+1],'first');
    
    % Compute best threshold for case 1
    [best_objective_case1, thresh_index_case1,feasible_case1, ncut_case1] = ... 
        find_best_thresh(cuts_thresh_case1, gvol_compl_thresh_case1, ... 
        gvol_thresh_case1,penalty_thresh_case1, index_unique);
    
    % Compute best threshold for case 2
    [best_objective_case2, thresh_index_case2,feasible_case2, ncut_case2] = ...
        find_best_thresh(cuts_thresh_case2, gvol_compl_thresh_case2, ... 
        gvol_thresh_case2, penalty_thresh_case2, index_unique);
        
    % check which is better
    if (best_objective_case1<best_objective_case2)
        best_objective=best_objective_case1;
        clusters=double(f_rest>=f_unique(thresh_index_case1));
        feasible=feasible_case1;
        ncut=ncut_case1;
    else
        best_objective=best_objective_case2;
        clusters=double(f_rest<f_unique(thresh_index_case2));
        feasible=feasible_case2;
        ncut=ncut_case2;
    end
       
end

% compute the best threshold according to the ncut objective with subset 
% incorporated and hvol constraint as penalty
function [best_objective,thresh_index,feasible, ncut] = ... 
    find_best_thresh(cuts_thresh, gvol_compl_thresh, ... 
        gvol_thresh, penalty_thresh, index_unique)

    % extract the values for thresholding at all the unique indices
    cuts_thresh=cuts_thresh(index_unique);
    gvol_compl_thresh = gvol_compl_thresh(index_unique);
    gvol_thresh = gvol_thresh(index_unique);
    penalty_thresh = penalty_thresh(index_unique);
   
    % compute balancing term (S)
    % this will be 0 if one of the sets is empty
    balance_thresh= gvol_compl_thresh.*gvol_thresh;
        
    % now compute total objective
    objective_thresh = (cuts_thresh + penalty_thresh)./balance_thresh;
   
    % find best one (this will automatically exclude the case where the 
    % set is empty since it cannot the minimum)
    [best_objective,thresh_index]=min(objective_thresh);
    
    % compute output
    feasible=penalty_thresh(thresh_index)==0;
    ncut=cuts_thresh(thresh_index)/balance_thresh(thresh_index);
    
end

