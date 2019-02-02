function [lambda,sg]=vol_cnstr_ncut_subset_penalty_functional(f,gamma,num,k,subset,deg,wval,ix,jx,gvol, totVol)

    [~,sortind]=sort(f);
    sgvol = gvol(sortind);
    sg = zeros(num,1);
    cumvols = sum(sgvol) - [0; cumsum(sgvol(1:num-1))];
    ii = find(cumvols > k);
    if ~isempty(ii) % otherwise subgradient is zero
        sg(sortind(ii)) = sgvol(ii);
        sg(sortind(ii(end)+1:num)) = 0;
        sg(sortind(ii(end))) = cumvols(ii(end))-k;
    end
   

    Pf = f - (f'*deg/sum(deg));
        
    lambda = ( 0.5*sum(wval.*abs(f(ix)-f(jx))) + gamma*f'*sg + gamma*max(f)*length(subset) - gamma* sum(f(subset))  ) / ...
                    ( 0.5*totVol* deg'*abs(Pf) );

end
