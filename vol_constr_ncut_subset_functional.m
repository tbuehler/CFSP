function [lambda,sg1,sg2]=vol_constr_ncut_subset_functional(f,gamma,num,k1,k2,deg,wval,ix,jx,gvol)


    [~,sortind]=sort(f);
    sgvol = gvol(sortind);
    sg2 = zeros(num,1);
    cumvols = sum(sgvol) - [0; cumsum(sgvol(1:num-1))];
    ii = find(cumvols > k2);
    sg2(sortind(ii)) = sgvol(ii);
    sg2(sortind(ii(end)+1:num)) = 0;
    sg2(sortind(ii(end))) = cumvols(ii(end))-k2;
   
    
    %indvec=zeros(num,1);indvec(sortind(1:(num-k2)))=1;
    sg1=zeros(num,1);sg1(sortind(num-k1+1:num))=1;

    
    lambda= (f'*gvol + gamma*f'*sg2 -gamma*f'*sg1 + k1*gamma*max(f))/(deg'*f - 0.5*sum(wval.*abs(f(ix)-f(jx))));

end
