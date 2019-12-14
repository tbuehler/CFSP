function [dc, min_cnstr_dens] = opt_thresh_cnstr_maxdens(W, f, g, ...
                                k1, k2, gam, assocJ, gvolJ, Wcut)
% Performs optimal thresholding of the constrained objective (with subset) for
% the (generalized) maximum density subgraph problem.
%
% Usage: 
% [dc, min_cnstr_dens]  = opt_thresh_cnstr_maxdens(W, f, g, k1, k2, gam, ...
%                                                  assocJ, gvolJ, Wcut)
%
% Input: 
% W         The weight matrix.
% f         The vector to be thresholded.
% g         Generalized degrees
% k1        Lower bound.
% k2        Upper bound.
% gam       Penalty parameter gamma.
% assocJ    Association in the seed set
% gvolJ     Sum of generalized degrees in seed set
% Wcut      Part of weight matrix representing edges between seed and rest
%
% Output:
% dc        Indicator vector obtained via optimal thresholding
%
% (C)2012-19 Thomas Buehler, Syama Rangapuram, Simon Setzer and Matthias Hein
   
    assert(size(f,1)==size(W,1), 'Input vector has wrong dimension.');
    assert(size(g,1)==length(f), 'Degree vector has wrong dimension.');
    assert(~issparse(g), 'Vector should not be sparse.');

    sizeJ = size(Wcut, 2);
    num = length(f);

    [~,ind_sort] = sort(f);
    g_sort = g(ind_sort);

    cum_cards = num:-1:1;
    cum_cards = cum_cards';

    Wt = W(ind_sort(1:num-1), ind_sort(1:num-1));
    temp = 2*sum(triu(Wt,1),2) + diag(Wt) ...
           + 2*sum(W(ind_sort(1:num-1), ind_sort(num:num)),2);
    all_assoc = sum(sum(W)) - [0;cumsum(temp)];
    all_assoc(all_assoc<0) = 0;

    if sizeJ>0
        Wcut = Wcut(ind_sort, :);
        cutJ = sum(Wcut, 2);
        all_cut = sum(cutJ) - [0;cumsum(cutJ)];
        all_cut = all_cut(1:end-1);

        all_assoc = all_assoc + assocJ + 2*all_cut;
    end

    penalty_thresholds = max(0, k1 - cum_cards) + max(0, cum_cards - k2);

    cumvols = sum(g_sort) - [0; cumsum(g_sort(1:num-1))] + gvolJ;

    cnstr_dens = cumvols + gam * penalty_thresholds;
    cnstr_dens = cnstr_dens./all_assoc;

    [min_cnstr_dens, ind] = min(cnstr_dens);
    
    dc = zeros(num,1);
    dc(ind_sort(ind:end)) = true;
end
