function [dc_best, maxdens_best, gam_best, dc_all, maxdens_all, gam_all] ...
    = community_detection(W, k1, k2, seed, gdeg, numRuns, dist_max, verbosity)
% Performs community detection on a graph/network by formulating the task
% as a constrained local (generalized) maximum density subgraph problem.
% The optimization problem is then solved using the method RatioDCA.
%
% Corresponding paper:
% T. Buehler, S. S. Rangapuram, S. Setzer and M. Hein
% Constrained fractional set programs and their application in local clustering 
% and community detection
% ICML 2013, pages 624-632 (Extended version: http://arxiv.org/abs/1306.3409)
%
%
% Usage: [dc_best, maxdens_best, gam_best, dc_all, maxdens_all, gam_all] ... 
%    = community_detection(W, k1, k2, seed, gdeg, numRuns, dist_max, verbosity)
%
% Input: 
% W             The weight matrix.
% k1            Lower bound
% k2            Upper bound
% seed          Index of seed set
%
% Optional Input:
% gdeg          Generalized degrees. Default is all ones vector.
% numRuns       Number of runs. Default is 10.
% dist_max      Maximum distance from the seed vertices. Default is 2.
% verbosity     Controls how much information is displayed [0-3]. Default is 1.
%
% Output:
% dc_best       Indicator vector of the community with largest density.
% maxdens_best  The corresponding density.
% gam_best      The gamma value used in the optimization.
% dc_all        Cell array of indicator vector for each gamma value.
% maxdens_all   The densities of the corresponding communities.
% gam_all       The corresponding gamma values.
%
%
% (C)2012-19 Thomas Buehler, Syama Rangapuram, Simon Setzer and Matthias Hein

    % default parameters
    if (nargin<8); verbosity = 1; end
    if (nargin<7); dist_max = 2; end
    if (nargin<6); numRuns = 10; end
    if (nargin<5); gdeg = ones(size(W,1),1); end

    % some checks 
    assert(issparse(W),'Error! W has to be sparse.');
    assert(k1<=k2,'Error! k1 cannot be larger than k2.');
    assert(k1<=size(W,1),'Error! k1 cannot be larger than the size of the graph.');
    assert(k2>=length(seed),'Error! k2 has to be larger than the seed set.');
    assert(dist_max>=1,'Error! dist has to be at least 1.');
    assert(numRuns>=1,'Error! numRuns has to be at least 1.');
    assert(size(gdeg,1)==size(W,1),'Error! gdeg has wrong format.');
    assert(~isempty(seed),'Error! Seed set has to be non-empty.');
    
    % extract neighbours of distance at most dist_max 
    if (verbosity>0)
        disp(['Restricting to neighborhood around seed set with distance ', ...
              'at most ', num2str(dist_max), '.']);
    end
    
    gdeg = full(gdeg);
    num = size(W,1);
    seed_vec = zeros(num,1);
    seed_vec(seed) = 1; 
    sub_vec = seed_vec';
           
    A = double(W~=0);
    sub_vec_new = sub_vec;
    for l=1:dist_max
        sub_vec_new = sub_vec_new * A;
        sub_vec = sub_vec + sub_vec_new;
    end
    ind_sub = find(sub_vec>0)';

    % extract subgraph
    W_sub = W(ind_sub, ind_sub);
    gdeg_sub = gdeg(ind_sub);
    seed_sub = find(seed_vec(ind_sub));
 
    % throw away big graph to save memory
    clear W;
    clear gdeg;
    
    % compute the optimal maximum density subgraph on the reduced graph
    if (verbosity>0)
        fprintf('Starting optimization.\n');
    end
    [dc_best, maxdens_best, gam_best, dc_all, maxdens_all, gam_all] ...
        = cnstr_maxdens_subset(W_sub, k1, k2, gdeg_sub, seed_sub, numRuns, verbosity);
   
    % check if feasible
    feasible = (sum(dc_best)<=k2 && sum(dc_best)>=k1 && ...
                sum(dc_best(seed_sub))==length(seed_sub));
    assert(feasible);
   
    % convert index vectors of reduced graph (within
    % neighbourhood of seed) to index vectors of original graph
    for l=1:length(dc_all)
        dc = zeros(num,1);
        dc(ind_sub) = dc_all{l};
        dc_all{l} = dc;
    end
    dc_best_temp = dc_best;
    dc_best = zeros(num,1);
    dc_best(ind_sub) = dc_best_temp;
end



function [dc_best, maxdens_best, gam_best, dc_all, maxdens_all, gam_all] ...
    = cnstr_maxdens_subset(W, k1, k2, gdeg, subset, numRuns, verbosity)
% Performs multiple runs with random initializations of the RatioDCA for
% the constrained maximum density problem (subset constraint and 
% cardinality constraints).
%
% Usage: [dc_best, maxdens_best, gam_best, dc_all, maxdens_all, gam_all] ... 
%    = cnstr_maxdens_subset(W, gdeg, k1, k2, subset, numRuns)
%
% Input: 
% W         The weight matrix.
% k1        Lower bound
% k2        Upper bound
% gdeg      Generalized degrees
% subset    Index of seed set
% numRuns   Number of runs.
% verbosity Controls how much information is displayed [0-3]. Default is 1.
%
% Output:
% dc_best       Indicator vector of the community with largest density
% maxdens_best  The corresponding density
% gam_best      The gamma value used in the optimization
% dc_all        Cell array of indicator vector for each initialization
% maxdens_all   The densities of the corresponding communities
% gam_all       The corresponding gamma values

    maxdens_best = 0;
    gamma_loop = true;
    gam = 1;
    l = 1;
    tic1 = tic;
    
    % gamma loop
    while (gamma_loop)
        lambda_all(l) = inf;
        
        for i=1:numRuns
            % run algorithm for current gamma
            start = randn(size(W,1),1);
            [dc, maxdens, lambda] = ratiodca_cnstr_maxdens(W, k1, k2, gdeg, ...
                                    start, subset, gam, verbosity);
            feasible = sum(dc>0)<=k2 && sum(dc>0)>=k1;
            
            % store result if feasible and better
            if (feasible)
                if maxdens>maxdens_best
                    gam_best = gam;
                    dc_best = dc;
                    maxdens_best = maxdens;
                    lambda_best = lambda;
                end
                if (verbosity>0)
                    fprintf('Found feasible solution with density %.4f for gamma=%.3f and starting point %d/%d.\n', ...
                        maxdens, gam, i, numRuns);
                end
                tic1 = tic;
            else
                toc1=toc(tic1);
                if toc1>30
                    if (verbosity>0)
                        fprintf('Current gamma is %.3f. Finished initialization %d/%d.\n', ...
                            gam, i, numRuns);
                    end
                    tic1 = tic;
                end
            end
            
            % store the one with best lambda for current gamma
            if (lambda < lambda_all(l))
                lambda_all(l) = lambda;
                dc_all{l} = dc;
                maxdens_all(l) = maxdens;
                gam_all(l) = gam;
            end
        end
        dc = dc_all{l};
        
        % the decision whether we increase gamma is based on the
        % feasibility of the one with best lambda
        feasible = sum(dc>0)<=k2 && sum(dc>0)>=k1;
        if (feasible) % it might not be the best one
            % now we are done for this run
            gamma_loop = false;
        % otherwise, try to find a better gamma
        else
            if gam<100
                gam = gam*2;
            else
                gamma_loop = false; % we skip this one
            end
        end
        l=l+1;
    end
    toc1 = toc(tic1);
    
    % in case we aborted all runs
    if sum(maxdens_all==0)
        error('Could not find feasible solution. Increase the number of runs.');
    end
     
    % check that solution is feasible
    assert(sum(dc_best)>=k1);
    assert(sum(dc_best)<=k2);
    assert(sum(sum(W(dc_best==1,dc_best==1)))/sum(gdeg(dc_best==1))==maxdens_best);
    
    % report result
    if (verbosity>0)
        fprintf('Best feasible result found has density %.4f.\n', maxdens_best);
    end
end
