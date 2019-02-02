function [dc_best, maxdens_best, gam_best, dc_all,maxdens_all,gam_all] ...
    = community_detection(W, k1,k2,seed,gdeg,numRuns,dist_max)
% Performs community detection on a graph/network by formulating the task
% as a constrained local (generalized) maximum density subgraph problem.
% The optimization problem is then solved using the method RatioDCA.
%
% Corresponding paper:
% T. Buehler, S. S. Rangapuram, S. Setzer and M. Hein
% Constrained fractional set programs and their application in 
% local clustering and community detection
% ICML 2013, pages 624-632
% (Extended version available at http://arxiv.org/abs/1306.3409)
%
% Usage: [dc_best, maxdens_best, gam_best, dc_all,maxdens_all,gam_all] ... 
%    = community_detection(W, k1,k2,seed,gdeg,numRuns,dist_max)
%
% Input: 
% W         The weight matrix.
% k1        Lower bound
% k2        Upper bound
% seed      Index of seed set
%
% Optional Input:
% gdeg      Generalized degrees. Default is all ones vector.
% numRuns   Number of runs. Default is 10.
% dist_max  Maximum distance from the seed vertices. Default is 2.
%
% Output:
% dc_best       Indicator vector of the community with largest density
% maxdens_best  The corresponding density
% gam_best      The gamma value used in the optimization
% dc_all        Cell array of indicator vector for each gamma value
% maxdens_all   The densities of the corresponding communities
% gam_all       The corresponding gamma values



    %% default parameters
    if nargin< 7
        dist_max = 2;
    end
    if nargin< 6
        numRuns=10;
    end
    if nargin< 5
        gdeg=ones(size(W,1),1);
    end
    

    %% some checks 
    assert(issparse(W),'Error! W has to be sparse.');
    assert(k1<=k2,'Error! k1 cannot be larger than k2.');
    assert(k1<=size(W,1),'Error! k1 cannot be larger than the size of the graph.');
    assert(k2>=length(seed),'Error! k2 has to be larger than the seed set.');
    assert(dist_max>=1,'Error! dist has to be at least 1.');
    assert(numRuns>=1,'Error! numRuns has to be at least 1.');
    assert(size(gdeg,1)==size(W,1),'Error! gdeg has wrong format.');
    assert(~isempty(seed),'Error! Seed set has to be non-empty.');
    
    
    %% extract neighbours of distance at most dist_max 
    disp(['Restricting to neighborhood around seed set with distance ', ...
          'at most ', num2str(dist_max), '.']);
    
    num=size(W,1);
    seed_vec=zeros(num,1);
    seed_vec(seed)=1; 
    sub_vec=seed_vec';
           
    A=double(W~=0);
    for l=1:dist_max
        sub_vec=sub_vec*A;
    end
    ind_sub=find(sub_vec>0)';
  

    % extract subgraph
    W_sub=W(ind_sub,ind_sub);
    gdeg_sub=gdeg(ind_sub);
    seed_sub=find(seed_vec(ind_sub));
   
 
    % throw away big graph to save memory
    clear W;
    clear gdeg;
 
    
    %% compute the optimal maximum density subgraph on the reduced graph
    fprintf('Starting optimization.\n');
    [dc_best, maxdens_best, gam_best, dc_all,maxdens_all,gam_all] ...
        = cnstr_maxdens_subset(W_sub,k1,k2,gdeg_sub, seed_sub,numRuns);
   
    % check if feasible
    feasible=(sum(dc_best)<=k2 && sum(dc_best)>=k1 && sum(dc_best(seed_sub))==length(seed_sub));
    assert(feasible);
   
    % convert index vectors of reduced graph (within
    % neighbourhood of seed) to index vectors of original graph
    for l=1:length(dc_all)
        dc=zeros(num,1);
        dc(ind_sub)=dc_all{l};
        dc_all{l}=dc;
    end
    dc_best_temp=dc_best;
    dc_best= zeros(num,1);
    dc_best(ind_sub)=dc_best_temp;


end
