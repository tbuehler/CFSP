function [v2, l2, ones_C] = eig_std_Laplacian(W,normalized,deg)
    W = triu(W);
    W = W+W';
    
    D=spdiags(sum(W,2),0,size(W,1),size(W,1));
    opts.disp=1;
    opts.tol = 1E-6;    
    if size(W,1) > 50000
        opts.maxit=50;
    end
     opts.issym = 1;
    if (normalized)
        [eigvec,eigval]= eigs(D-W, spdiags(deg,0,size(W,1),size(W,1)),2,'SA',opts);
    else
        [eigvec,eigval]= eigs(D-W, 2,'SA',opts);
    end
    v2 = eigvec(:,2);
    l2 = eigval(2,2);

    ones_C = v2;
    
end
