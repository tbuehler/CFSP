function [clusters, ncut, feasible, lambda, solutions] = constrained_ncut(...
         W, k, h, subset, nRuns, verbosity)
% Solves the constrained normalized cut problem with an upper bound on the
% (generalized) volume as well as a subset constraint.
%
% The subset constraint is directly incorporated into the objective (leading to 
% a problem of lower dimension) and the volume constraint is incorporated into 
% the objective as penalty term. The optimization problem is then solved using 
% the method RatioDCA.
%
% Corresponding paper:
% T. Buehler, S. S. Rangapuram, S. Setzer and M. Hein
% Constrained fractional set programs and their application in local clustering 
% and community detection
% ICML 2013, pages 624-632 (Extended version: http://arxiv.org/abs/1306.3409)
%
%
% Usage: [clusters, ncut, feasible, lambda, solutions] = ...
%        constrained_ncut(W, k, h, subset, nRuns, verbosity)
%
% Input:
% W             Weight matrix (full graph).
% k             Upper bound on the (generalized) volume.
% h             Generalized degree vector used in constraint.
% subset        Indices of seed subset.
% nRuns         Number of runs with random initialization. Default is 10.
% verbosity     Controls how much information is displayed [0-3]. Default is 1.
%
% Output:
% clusters      Thresholded vector f yielding the best objective.
% ncut          Ncut value of resulting clustering.
% feasible      True if all constraints are fulfilled.
% lambda        Corresponding objective value.
% solutions     Obtained clusters for all starting points.
%
%
% (C)2012-19 Thomas Buehler, Syama Rangapuram, Simon Setzer and Matthias Hein

    if (nargin<6); verbosity = 1; end
    if (nargin<5); nRuns = 10; end

    cur_gamma = 0;
    gamma_diff = 0.05;
    num_gammas = 1;
    satisfied = false;
    num = size(W,1);
    perturbation = false;    
    starts = randn(num, nRuns);
        
    lambdas = zeros(nRuns, num_gammas);
    ncuts_temp = zeros(nRuns, num_gammas);
    feasibles = zeros(nRuns, num_gammas);
    solutions = zeros(num, nRuns, num_gammas);
               
    stop_runs = zeros(nRuns,1);
 
    while ~satisfied % while not all constraints are satisfied... gamma loop
        % we have the gamma iteration as outer loop since it allows us to stop 
        % early if we find very good solution for smaller gamma without 
        % exploring the gamma path for other starts
        for it=1:nRuns % we try from 10 different starting points 
            if stop_runs(it)
                clusters = solutions(:, it, num_gammas-1);
                assert(norm(clusters - starts(:, it))==0);
                lambda = lambdas(it, num_gammas-1);
                ncut = ncuts_temp(it, num_gammas-1);
            else
                [clusters, ncut, feasible, lambda] = solve_one_run_with_pert(...
                        W, k, h, starts(:, it), subset, cur_gamma, perturbation, verbosity);
                assert (isempty(setdiff(subset, find(clusters))));
                assert ((sum(h(clusters>0)) <= k) == feasible); 
            end

            % Now save the results for the current gamma... we have nRuns solutions
            lambdas(it, num_gammas) = lambda;
            ncuts_temp(it, num_gammas) = ncut;
            feasibles(it, num_gammas) = feasible;
            solutions(:, it, num_gammas) = clusters;
            assert(isempty(setdiff(subset, find(clusters))));
            
            % sanity check    
            if k==sum(h)
                assert(feasibles(it, num_gammas) == 1);
            end
            
            if (verbosity >0)  
                fprintf("Finished run %d for gamma=%.2f: ncut=%.5g lambda=%.5g, feasible=%d\n", ...
                        it, cur_gamma, ncut, lambda, feasible);
            end
        end

        % Let us discard some of the runs if their objectives are larger than 
        % the current best feasible objective.
        [min_feas_lambda, ~] = min(lambdas(feasibles(:, num_gammas)==1, num_gammas));
        if cur_gamma > 0 && ~isempty(min_feas_lambda)
            stop_runs = lambdas(:, num_gammas) >= min_feas_lambda;
        end
            
        % Now choose the minimum objective among the nRuns
        [lambda, min_ix] = min(lambdas(:, num_gammas));
        feasible = feasibles(min_ix, num_gammas);
        clusters = solutions(:, min_ix, num_gammas);

        if (sum(h(clusters>0)) > k) % size constraint is violated.
            satisfied = false;
            if (verbosity>0)
                fprintf("Size constraint is violated. Increasing gamma. Best lambda=%.5g (run %d).  \n", ...
                        lambda, min_ix);
            end
            starts = solutions(:, :, num_gammas);
            gamma_diff = gamma_diff*2;                
            cur_gamma = cur_gamma + gamma_diff;
            num_gammas = num_gammas + 1;        
        else 
            assert(isempty(setdiff(subset, find(clusters)))) % seed constraint needs to be fulfilled.
            satisfied= true;
            % save the current best result!
            ncut = ncuts_temp(min_ix, num_gammas);
            assert ( abs(ncut - balanced_cut(W, sum(W,2), clusters)) <= 1e-8 ); 
            if (verbosity>0)
                fprintf("Found feasible solution for bound=%d using gamma=%.2f (run %d): ncut=%.5g lambda=%.5g\n",  ... 
                        k, cur_gamma, min_ix, ncut, lambda);
            end
        end
    end 
end



function [clusters, ncut, feasible, lambda] = solve_one_run_with_pert(W, ...
         k, h, start, subset, gam, perturbation, verbosity)
    [clusters, ncut, feasible, lambda] = ratiodca_cnstr_ncut_direct(W, k, h, ...
                                         start, subset, gam, verbosity);                     
    if perturbation && gam > 0
        pert_probs = [0.01:0.01:0.1 0.2:0.1:1];
        for p=1:length(pert_probs) % try different perturbation probabilities 
            pert_counter = 0;
            while pert_counter < 1 % for each prob, try until a perturbation does not improve the result
                pert_start = clusters + pert_probs(p)* range(clusters) * randn(size(clusters,1),1);
                [cluster_pert, ncut_pert, feasible_pert, lambda_pert] = ...
                  ratiodca_cnstr_ncut_direct(W, k, h, pert_start, subset, gam, verbosity);
                if (verbosity>1)
                    fprintf("... Original lambda: %.5g After perturbation with p=%.2f: %.5g\n", ...
                            lambda, pert_probs(p), lambda_pert);                 
                end                
                if lambda_pert < lambda
                    clusters = cluster_pert;
                    ncut = ncut_pert;
                    feasible = feasible_pert;
                    lambda = lambda_pert;
                else
                    pert_counter = pert_counter + 1;
                end
            end
        end
    end
end
