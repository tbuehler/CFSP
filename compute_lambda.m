function [lambda, indvec2, indvec1] = compute_lambda(f, gam, k1prime, ...
 k2prime, deg_rest, gvolJ, degJ, assocJ, wval_rest, ix_rest, jx_rest, gdeg_rest)
% Computes the objective of the constrained maximum density problem with
% subset and cardinality constraints
%
% Usage: [lambda, indvec, indvec1] = compute_lambda(f, gam, k1prime, ...
% k2prime, deg_rest, gvolJ, degJ, assocJ, wval_rest, ix_rest, jx_rest, gdeg_rest)
%
% Input: 
% f         The current vector
% gam       The penalty parameter gamma
% k1prime   Lower bound (adjusted for seed set)
% k2prime   Upper bound (adjusted for seed set)
% deg_rest  Degrees (dimension of graph without seeds) 
% gvolJ     Generalized volume of seed subset
% degJ      For each vertex: sum of edge weights to seed subset.
% assocJ    Association inside seed.
% wval_rest The values of the reduced weight matrix
% ix_rest   row indices of reduced weight matrix 
% jx_rest   column indices
% gdeg_rest Generalized degrees
%
% Output:
% Lambda    The objective value
% indvec2   1 - indicator vector of k2 largest components
% indvec1   indicator vector of k1 largest components
    
    [~, sortind] = sort(f);
    num = length(f);
    indvec1 = zeros(num,1);
    indvec1(sortind(num-k1prime+1:num)) = 1;
    indvec2 = zeros(num,1);
    indvec2(sortind(num-k2prime+1:num)) = 1;

    R1 = f'*gdeg_rest + gam*f'* ones(length(f),1) + ...
         (gvolJ + gam * k1prime)*max(f);
    R2 = gam*f'* (indvec1 + indvec2);
    S1 = f'*(deg_rest+degJ) + assocJ*max(f);
    S2 = 0.5*sum(wval_rest.*abs(f(ix_rest)-f(jx_rest)));
    
    lambda = (R1-R2) / (S1-S2);
    
    indvec2 = ones(length(indvec2),1)-indvec2;
end
