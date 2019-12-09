function bcut = balanced_cut(W, gdeg, Y)
% Computes the balanced cut where the balancing term uses generalized degrees
% gdeg. The clustering is specified by Y. Also works for multiple partitions.
%
% Usage: bcut = balanced_cut(W, gdeg, Y)
%
% (C)2012-19 Thomas Buehler, Syama Rangapuram, Simon Setzer and Matthias Hein

    classes = unique(Y);
    k = length(classes);
    bcut = 0;
    n = size(W,1);

    for i=1:k
        jx = find(Y==classes(i));
        jxc = setdiff(1:n, jx);

        cut = sum(sum(W(jx, jxc)));
        vol = sum(gdeg(jx));

        bcut = bcut + cut/vol;
    end
end
