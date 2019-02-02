function ccut = compute_cheeger_cut(W, gdeg, Y)
    labels = unique(Y);
    
    ix = find(Y == labels(1) );
    ixc = find(Y ~= labels(1) );
    
    cut = sum(sum(W(ix, ixc)));
    vol = sum(gdeg(ix));
    volc = sum(gdeg(ixc));
    
    ccut = cut/ min(vol, volc);

end
