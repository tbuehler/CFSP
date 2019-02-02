function [allClustersInClusterM, cutPart1,cutPart2] = ...
    vol_seed_cnstr_opt_thresh_ncut(seed, vmin_comp,W, vertex_weights_comp, ...
    criterion_threshold,vol2, index_comp,index_m,cut_rest,size_rest,size_m)
      
        if nargin < 7
            
            index_comp = 1:size(W,1);
            index_m = 1:size(W,1);
            cut_rest = 0;
            size_rest = 0;
            size_m = size(W,1);
            
        end
        
        n = size(W,1);
        deg=full(sum(W));
        W_comp = W(index_comp, index_comp);
        
        [vminM_sorted, index]=sort(vmin_comp);
        seed_index = find(index == seed);
%         [vminU,indexU]=unique(vminM_sorted); % if you change to unique,
%         then consider checking if the seed does not cause problems.
        vminU = vminM_sorted; indexU = 1:size(W,1);
        
        
        W_sorted=W_comp(index,index);
        vertex_weights_comp = vertex_weights_comp(index);
        balance_thresholds = cumsum(vertex_weights_comp);

        % calculate cuts
        deg_comp=deg(index_m(index_comp));
        volumes_threshold=cumsum(deg_comp(index));
        triup=triu(W_sorted);
        tempcuts_threshold=volumes_threshold - 2*cumsum(full(sum(triup)));
        tempcuts_threshold2=(volumes_threshold(end)-volumes_threshold) - (sum(sum(W_sorted))-2*cumsum(full(sum(triup,2)))');            

        % enforce that seed comes in the subset
        tempcuts_threshold(seed_index: end) = inf;
        tempcuts_threshold2(seed_index: end) = inf;
        
        tempcuts_threshold=tempcuts_threshold(indexU);
        tempcuts_threshold2=tempcuts_threshold2(indexU);

        %volumes_threshold=volumes_threshold(indexU); %balance_thresholds will be used instead of volumes_thresholds
        balance_thresholds = balance_thresholds(indexU);
        
        
        % divide by size/volume
        %if(normalized)
        % Option 1: merge the remaining components and C_t.
        % In this case, we need to count the constraints voilated
        % between the components C_t^prime and the merged part.
        cutparts1_threshold=(tempcuts_threshold(1:end-1)+cut_rest)./(balance_thresholds(1:end-1)+size_rest);
        cutparts2_threshold=tempcuts_threshold2(1:end-1)./(balance_thresholds(end)-balance_thresholds(1:end-1));

        % Option 2: merge the remaining components adn C_t^prime
        % In this case, we need to count the constraints voilated
        % between the components C_t and the merged part.
        cutparts1b_threshold=tempcuts_threshold(1:end-1)./balance_thresholds(1:end-1);
        cutparts2b_threshold=(tempcuts_threshold2(1:end-1)+cut_rest)./((balance_thresholds(end)-balance_thresholds(1:end-1))+size_rest);
            
%         sgvol = vertex_weights_comp(index);
%         sgvol = sgvol(indexU);
%         sgvol = sgvol';

        sgvol = vertex_weights_comp';
        cumvols = sum(sgvol) - [0; cumsum(sgvol(1:n-1))];
        penalty_thresholds = max(0, cumvols-vol2);
        penalty_thresholds(penalty_thresholds > 0) = inf;
        penalty_thresholds = penalty_thresholds(2:end);
        penalty_thresholds = penalty_thresholds';
        
        %calculate total cuts
        if(criterion_threshold==1)

            cuts_threshold=cutparts1_threshold+cutparts2_threshold + penalty_thresholds;
            [cut1,threshold_index]=min(cuts_threshold);
            
            
            cuts_threshold_b=(cutparts1b_threshold+cutparts2b_threshold) + penalty_thresholds;
            [cut1b,threshold_index_b]=min(cuts_threshold_b);
                        
            comp_case=1;
            if (cut1b<cut1) 
                comp_case=2;
            end
        else
            cheegers_threshold=max(cutparts1_threshold,cutparts2_threshold) + penalty_thresholds;
            [cheeger1,threshold_index]=min(cheegers_threshold);
            
            cheegers_threshold_b=max(cutparts1_threshold,cutparts2_threshold) + penalty_thresholds;
            [cheeger1b,threshold_index_b]=min(cheegers_threshold_b);
            
            comp_case=1;
            if (cheeger1b<cheeger1) 
                comp_case=2;
            end
        end

        if(comp_case==1)
            cutPart1=cutparts1_threshold(threshold_index);
            cutPart2=cutparts2_threshold(threshold_index);
            
            allClustersInClusterM_comp= (vmin_comp>vminU(threshold_index));
        
            allClustersInClusterM= zeros(size_m,1);
            allClustersInClusterM(index_comp)=allClustersInClusterM_comp;
        else
            cutPart1=cutparts1b_threshold(threshold_index_b);
            cutPart2=cutparts2b_threshold(threshold_index_b);
            allClustersInClusterM_comp= (vmin_comp>vminU(threshold_index_b));
        
            allClustersInClusterM= ones(size_m,1);
            allClustersInClusterM(index_comp)=allClustersInClusterM_comp;
        end
        
%    cuts_threshold

    if sum(sgvol(allClustersInClusterM(index)>0)) > vol2 % volume constraint is violated while enforcing seed constraint!
        allClustersInClusterM = zeros(n, 1);
        cutPart1 = NaN;
        cutPart2 = NaN;
    end

end

