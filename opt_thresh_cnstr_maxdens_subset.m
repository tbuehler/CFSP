function dc = opt_thresh_cnstr_maxdens_subset(W, f, gdeg, k1, k2, gamma,assocJ,gvolJ,Wcut)
% Performs optimal threwsholding of the constraint objective (with subset)
%
% Usage: dc = opt_thresh_cnstr_maxdens_subset(W, f, gdeg, k1, k2, gamma,assocJ,gvolJ,Wcut)
%
% Input: 
% W         The weight matrix.
% f        	The vector to be thresholded.
% gdeg      Generalized degrees
% k1        Lower bound.
% k2        Upper bound.
% gamma     Penalty parameter gamma.
% assocJ    Association in the seed set
% gvolJ     Sum of generalized degrees in seed set
% Wcut      Part of weight matrix representing edges between seed and rest
%
% Output:
% dc        Indicator vector obtained via optimal thresholding
% 

    sizeJ=size(Wcut,2);

    num = length(f);
    [~,indsort]=sort(f);
    sgdeg = gdeg(indsort);

    cum_cards = num:-1:1;
    cum_cards = cum_cards';

    Wt=W(indsort(1:num-1),indsort(1:num-1));
    Wtu=triu(Wt,1);
    temp=2*sum(Wtu,2)+diag(Wt) +2 * sum(W(indsort(1:num-1),indsort(num:num)),2);
    all_assoc=sum(sum(W(indsort(1:end),indsort(1:end))))-[0;cumsum(temp)];
    all_assoc(all_assoc<0) = 0;


    if sizeJ>0
        Wcut= Wcut(indsort(1:num),:);
        cutJ= sum(Wcut,2);
        all_cut= sum(cutJ)-[0;cumsum(cutJ)];
        all_cut=all_cut(1:end-1);

        all_assoc=all_assoc+assocJ + 2*all_cut;
    end

    penalty_thresholds = max(0, k1 - cum_cards) + max(0, cum_cards-k2);

    cumvols = sum(sgdeg) - [0; cumsum(sgdeg(1:num-1))] + gvolJ;

    cnstr_dens=cumvols+ gamma* penalty_thresholds;
    cnstr_dens = cnstr_dens./all_assoc;

    ind3=find(cnstr_dens==min(cnstr_dens),1);

    dc=zeros(num,1);
    dc(indsort(ind3:end))=true;

end
