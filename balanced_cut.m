function bcut = balanced_cut( W, gdeg, Y)
% Computes balanced cut where the balance is specified by gdeg
% and the partition by Y. Also works for multiple partitions.
%
% Usage: bcut = balanced_cut( W, gdeg, Y)

    classes = unique(Y);
    k = length(classes);

    bcut = 0;
    n = size(W,1);

    for i=1:k

        jx = find(Y==classes(i));
        jxc = setdiff(1:n, jx);

        cut = sum(sum( W(jx, jxc) ));
        vol = sum(gdeg(jx));

        bcut = bcut + cut/vol;

    end

end
