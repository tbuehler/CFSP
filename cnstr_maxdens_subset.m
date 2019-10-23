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

    if (nargin<7) verbosity = 1; end

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
            [dc, maxdens, lambda] = cnstr_maxdens_subset_single_run(W, k1, ...
                                    k2, gdeg, start, subset, gam, verbosity);
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
