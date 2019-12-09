function [cluster_grid, ncut_grid, vol_grid, subs_grid] = soft_constrained_ncut(...
         W, k, h, subset, gammas_subs, gammas_vol, numRuns, verbosity)
% Solves the constrained normalized cut problem with an upper bound on the
% (generalized) volume as well as a subset constraint.
%
% The subset constraint as well as the volume constraint are treated as soft 
% constraints which are controlled by corresponding parameters gamma_sub and 
% gamma_vol. The optimization problem is solved for all provided values of 
% gamma_sub and gamma_vol using the method RatioDCA.
%
% Corresponding paper:
% T. Buehler, S. S. Rangapuram, S. Setzer and M. Hein
% Constrained fractional set programs and their application in local clustering 
% and community detection
% ICML 2013, pages 624-632 (Extended version: http://arxiv.org/abs/1306.3409)
%
%
% Usage:  
% [cluster_grid, ncut_grid, vol_grid, subs_grid] = soft_constrained_ncut(W, ...
%                     k, h, subset, gammas_subs, gammas_vol, starts, verbosity);
%
% Input:
% W             Weight matrix (full graph).
% k             Upper bound on the (generalized) volume.
% h             Generalized degree vector used in constraint.
% subset        Indices of seed subset.
% gammas_subs   Penalty parameter for subset constraint.
% gammas_vol    Penalty parameter for volume constraint.
% numRuns       Number of runs. Default is 10.
% verbosity     Controls how much information is displayed [0-3]. Default is 1.
%
% Output:
% cluster_grid  Grid containing clusterings for all values of gamma1 and gamma2.
% ncut_grid     Corresponding ncut values.
% vol_grid      Corresponding volumes.
% subs_grid     Corresponding sizes of intersections with subsets.
%
%
% (C)2012-19 Thomas Buehler, Syama Rangapuram, Simon Setzer and Matthias Hein

    if (nargin<8); verbosity = 1; end
    if (nargin<7); numRuns = 10; end

    num = size(W,1);
    starts = randn(num, numRuns);
    perturbation = false;

    lambdas = zeros(numRuns, length(gammas_subs), length(gammas_vol));
    ncuts_temp = zeros(numRuns, length(gammas_subs), length(gammas_vol));
    feasibles = zeros(numRuns, length(gammas_subs), length(gammas_vol));
    solutions = zeros(num, numRuns, length(gammas_subs), length(gammas_vol));
 
    cluster_grid = cell(length(gammas_subs), length(gammas_vol));
    ncut_grid = zeros(length(gammas_subs), length(gammas_vol));
    vol_grid = zeros(length(gammas_subs), length(gammas_vol));
    subs_grid = zeros(length(gammas_subs), length(gammas_vol));
                
    for i=1:length(gammas_subs)
        gamma1 = gammas_subs(i);
        for j=1:length(gammas_vol)
            gamma2 = gammas_vol(j);
            if (verbosity>0)
                fprintf("Starting to solve NCut problem for gamma1=%.2f gamma2=%.2f:\n", ...
                        gamma1, gamma2); 
            end            
            for it=1:numRuns % we try from numRuns different starting points 
                [clusters, ncut, feasible, lambda] = ... 
                  solve_one_run_with_pert(W, k, h, starts(:, it), subset, ...
                                          gamma1, gamma2, perturbation, verbosity); 
                assert(abs(ncut-balanced_cut(W, sum(W,2), clusters))<=1E-12);              
    
                % Save the results for the current gamma1 and gamma2.
                lambdas(it, i, j) = lambda;
                ncuts_temp(it, i,j) = ncut;
                feasibles(it, i,j) =  feasible; 
                solutions(:, it, i,j) = clusters;
                if (verbosity >0)  
                    fprintf("Finished run %d for gamma1=%.2f, gamma2=%.2f: ncut=%.5g lambda=%.5g, feasible=%d\n", ...
                            it, gamma1, gamma2, ncut, lambda, feasible);
                end
            end
                             
            %% Now choose the minimum objective among the numRuns
            [min_obj, min_ix] = min(lambdas(:, i,j));
            cur_feasibles = feasibles(min_ix, i,j);
            clusters = solutions(:, min_ix, i,j);
            ncut = ncuts_temp(min_ix, i,j);
                
            cluster_grid{i,j} = clusters;
            ncut_grid(i,j) = ncut;
            vol_grid(i,j)=  sum(sum(W(clusters==1,:)));
            subs_grid(i,j) = sum(clusters(subset)); %size of intersection with subset 
            
            if (verbosity >0)
                fprintf("Best run: %d lambda=%.5g, feasible=%d\n", ...
                        min_ix, min_obj, cur_feasibles);
            end  
        end
    end            
end

function [clusters, ncut, feasible, lambda] = ...
         solve_one_run_with_pert(W, k, h, start, subset, gamma1, gamma2, perturbation, verbosity)
    % scaling factor for penalty terms
    pen_vol_max=0.5*sum(sum(W))-k;
    pen_subs_max=length(subset);
    scaling=pen_vol_max/pen_subs_max;
    [clusters, ncut, feasible, lambda] = ratiodca_cnstr_ncut_penalty(W, k, ...
                           h, start, subset, gamma1*scaling, gamma2, verbosity);
    if perturbation && gamma1 > 0 && gamma2>0
        pert_probs = [0.01:0.01:0.1 0.2:0.1:1];
        for p=1:length(pert_probs) % try different perturbation probabilities 
            pert_counter = 0;
            while pert_counter < 1 % try until a perturbation does not improve the result
                pert_start = clusters + pert_probs(p)* range(clusters) * randn(n,1);
                [cluster_pert, ncut_pert, feasible_pert, lambda_pert] = ...
                  ratiodca_cnstr_ncut_penalty(W, k, h, pert_start, subset, gamma1, gamma2, verbosity); 
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
