function [lambda,sg]=functional_vol_cnstr_ncut(f,gamma,num,k,deg,wval,ix,jx,hdeg, totVol)
% Computes the functional associated to the constrained normalized
% cut problem with an upper bound on the generalized volume, treated as a
% penalty term with parameter gamma.
%
% Usage: [lambda,sg]=functional_vol_cnstr_ncut_subset(f, ... 
%       gamma, num,k,deg,wval,ix,jx,hdeg, totVol)
%
% Input:
% f             Input vector.
% gamma         Parameter corresponding to volume constraint.
% num           Length of f.
% k             Upper bound on generalized volume.
% deg           Degree vector.
% wval          The values of the weight matrix.
% ix            The corresponding row indices.
% jx            The corresponding column indices.
% hdeg          Generalized degrees.
% totVol        Total volume of graph.
%
% Output:
% lambda        Value of the functional.
% sg            Subgradient of penalty function for the volume constraint.

    % make sure this is used correctly
    assert(size(f,1)==num,'Input vector has wrong dimension.');
    assert(size(deg,1)==length(f),'Degree vector has wrong dimension.');
    assert(size(hdeg,1)==length(f),'Degree vector has wrong dimension.');   

    % Compute subgradient
    [fsort,sortind]=sort(f);
    shdeg = hdeg(sortind);
    sg = zeros(num,1);
    cumvols = sum(shdeg) - [0; cumsum(shdeg(1:num-1))];
    ii = find(cumvols > k);
    if ~isempty(ii) % otherwise subgradient is zero
        sg(sortind(ii)) = shdeg(ii);
        sg(sortind(ii(end)+1:num)) = 0;
        sg(sortind(ii(end))) = cumvols(ii(end))-k;
    end
       
    % Evaluate objective
    Pf = f - (f'*deg/sum(deg));
    R = 0.5*sum(wval.*abs(f(ix)-f(jx))) + gamma*f'*sg; % sg is g-t
    S =  0.5*totVol* deg'*abs(Pf);
    lambda = R / S;

end
