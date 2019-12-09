function [lambda, sg] = lambda_cnstr_ncut_penalty(f, gamma1, gamma2, num, k, ...
                        subset, deg, wval, ix, jx, hdeg, totVol)
% Computes the functional associated to the constrained normalized
% cut problem with subset constraint and upper bound on the generalized
% volume. Both constraints are incorportated via separate penalties
% with parameters gamma1 and gamma2.
%
% Usage: [lambda,sg] = lambda_cnstr_ncut_penalty(f, gamma1, gamma2, num, k, ...
%                      subset, deg, wval, ix, jx, hdeg, totVol)
%
% Input:
% f        Input vector.
% gamma1   Parameter corresponding to subset constraint.
% gamma2   Parameter corresponding to volume constraint.
% num      Length of f.
% k        Upper bound on generalized volume.
% subset   Indices of seed subset.
% deg      Degree vector.
% wval     The values of the weight matrix.
% ix       The corresponding row indices.
% jx       The corresponding column indices.
% hdeg     Generalized degrees.
% totVol   Total volume.
%
% Output:
% lambda   Value of the functional.
% sg       Subgradient of the penalty function for the volume constraint.
%
% (C)2012-19 Thomas Buehler, Syama Rangapuram, Simon Setzer and Matthias Hein
 
    % make sure this is used correctly
    assert(size(f,1)==num,'Input vector has wrong dimension.');
    assert(size(deg,1)==length(f),'Degree vector has wrong dimension.');
    assert(size(hdeg,1)==length(f),'Degree vector has wrong dimension.');   

    % Compute subgradient
    [~, sortind] = sort(f);
    shdeg = hdeg(sortind);
    sg = zeros(num,1);
    cumvols = sum(shdeg) - [0; cumsum(shdeg(1:num-1))];
    ii = find(cumvols > k);
    if ~isempty(ii) % otherwise subgradient is zero
        sg(sortind(ii)) = shdeg(ii); % the last one will be overwritten
        sg(sortind(ii(end)+1:num)) = 0;
        sg(sortind(ii(end))) = cumvols(ii(end))-k;
    end
       
    % Evaluate functional
    Pf = f - (f'*deg/sum(deg));
    R1 = 0.5*sum(wval.*abs(f(ix)-f(jx))) +  gamma1*max(f)*length(subset);
    R2 = gamma1* sum(f(subset)) - gamma2*f'*sg ; % sg is hdeg - t(fk)
    S = 0.5*totVol* deg'*abs(Pf);
    lambda = ( R1 - R2  ) / S; 
end
