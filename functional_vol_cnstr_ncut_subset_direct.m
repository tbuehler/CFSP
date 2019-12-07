function [lambda, sg] = functional_vol_cnstr_ncut_subset_direct(f, gam, ...
        num, kprime, deg, wval_rest, ix_rest, jx_rest, hdeg, totVol_rest, deg_S)
% Computes the functional associated to the constrained normalized
% cut problem with subset constraint and upper bound on the generalized
% volume. The volume constraint is treated as a penalty term with 
% parameter gam, while the subset constraint is incorporated directly
% into the objective.
%
% Usage: [lambda, sg] = functional_vol_cnstr_ncut_subset_direct(f, gam, ...
%       num, kprime, deg, wval_rest, ix_rest, jx_rest, hdeg, totVol_rest, deg_S)
%
% Input:
% f             Input vector.
% gam           Parameter corresponding to volume constraint.
% num           Length of f (dimension of graph without seed set).
% kprime        Upper bound on generalized volume (adjusted for seed set).
% deg           Degree vector.
% wval_rest     The values of the (restricted) weight matrix.
% ix_rest       The corresponding row indices.
% jx_rest       The corresponding column indices.
% hdeg          Generalized degrees.
% totVol_rest   Total volume of graph without subset.
% totVol        Total volume of graph
% deg_S         For each vertex: sum of weights of edges to seed set.
%
% Output:
% lambda        Value of the functional.
% sg            Subgradient of penalty function for the volume constraint.

    % make sure this is used correctly
    assert(size(f,1)==num,'Input vector has wrong dimension.');
    assert(size(deg,1)==length(f),'Degree vector has wrong dimension.');
    assert(size(hdeg,1)==length(f),'Degree vector has wrong dimension.');   

    % Compute subgradient
    [fsort, sortind] = sort(f);
    shdeg = hdeg(sortind);
    sg = zeros(num,1);
    cumvols = sum(shdeg) - [0; cumsum(shdeg(1:num-1))];
    ii = find(cumvols > kprime);
    if ~isempty(ii)
        sg(sortind(ii)) = shdeg(ii);
        sg(sortind(ii(end)+1:num)) = 0;
        sg(sortind(ii(end))) = cumvols(ii(end))-kprime;
    end
     
    % Compute functional
    Pf = f - (f'*deg/sum(deg));
    R1 = 0.5*sum(wval_rest.*abs(f(ix_rest)-f(jx_rest))) + sum(deg_S)*max(f);
    R2 = deg_S'*f - gam*f'*sg ; % sg is g - t
    S1 = 0.5*totVol_rest* deg'*abs(Pf) + sum(deg_S)*totVol_rest*max(f);
    S2 = sum(deg_S) * deg'*f;  
    lambda = (R1-R2) / (S1-S2); 
end
