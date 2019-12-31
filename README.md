# Constrained fractional set programming



This archive contains a Matlab implementation of the CFSP algorithm for the
constrained normalized cut and densest subgraph problems as described in the paper
 
    T. Buehler, S. S. Rangapuram, S. Setzer and M. Hein
    Constrained fractional set programs and their application in local clustering and community detection
    ICML 2013, pages 624-632
    (Extended version available at http://arxiv.org/abs/1306.3409)

The implementation uses mex-files to solve the inner problem of the RatioDCA. 
Compile them by typing 'make' in the matlab command line.



## Documentation


### Constrained densest subgraph

Performs community detection on a graph/network by treating the task
as a constrained local (generalized) maximum density subgraph problem with
a subset constraint as well as upper and lower bounds on the cardinality.

#### Usage: 

    [dc_best, maxdens_best, gam_best, dc_all, maxdens_all, gam_all] = community_detection(... 
                                         W, k1, k2, seed, gdeg, numRuns, dist_max, verbosity)

#### Input: 

    W               The weight matrix.
    k1              Lower bound on size of solution.
    k2              Upper bound on size of solution.
    seed            Index of seed set.
    gdeg            Generalized degrees (vertex weights) used in denominator of 
                    the density. Default is all ones vector.
    numRuns         Number of runs. Default is 10.
    dist_max        Maximum distance from the seed vertices. Default is 2.
    verbosity       How much information is displayed [0-3]. Default is 1.

#### Output:

    dc_best         Indicator vector of the community with largest density.
    maxdens_best    The corresponding density.
    gam_best        The gamma value used in the optimization.
    dc_all          Cell array of indicator vector for each value of gamma.
    maxdens_all     The densities of the corresponding communities.
    gam_all         The corresponding gamma values.



### Constrained normalized cut

Solves the constrained normalized cut problem with a subset constraint and an 
upper bound on the (generalized) volume. Internally, the subset constraint is 
directly incorporated into the objective (leading to a problem of lower 
dimension) and the volume constraint is incorporated as a penalty term. 

#### Usage:

    [clusters, ncut, feasible, lambda, solutions] = constrained_ncut(W, k, h, subset, nRuns, verbosity)

#### Input:
    W               Weight matrix (full graph).
    k               Upper bound on the (generalized) volume.
    h               Generalized degrees (vertex weights) used in constraint.
    subset          Indices of seed subset.
    nRuns           Number of runs with random initializations. Default is 10.
    verbosity       How much information is displayed [0-3]. Default is 1.      


#### Output:

    clusters        Thresholded vector f yielding the best objective.
    ncut            Ncut value of resulting clustering.
    feasible        True if all constraints are fulfilled.
    lambda          Corresponding objective.
    solutions       Obtained clusters for all starting points.


### Soft-constrained normalized cut

Solves the constrained normalized cut problem with a subset constraint and an 
upper bound on the (generalized) volume. Both constraints are treated as soft 
constraints controlled by corresponding parameters gamma_sub and gamma_vol. 
The optimization problem is solved for all given values of gamma_sub and gamma_vol.

#### Usage: 
    
    [cluster_grid, ncut_grid, vol_grid, subs_grid] = soft_constrained_ncut(...
                  W, k, h, subset, gammas_subs, gammas_vol, numRuns, verbosity)

#### Input:
    
    W               Weight matrix (full graph).
    k               Upper bound on the (generalized) volume.
    h               Generalized degres (vertex weights) used in constraint.
    subset          Indices of seed subset.
    gammas_subs     Array of penalty parameters for subset constraint.
    gammas_vol      Array of penalty parameters for volume constraint.
    numRuns         Number of runs. Default is 10.
    verbosity       How much information is displayed [0-3]. Default is 1.

#### Output:

    cluster_grid    Grid containing clusterings for all values of gammas_subs 
                    and gammas_vol.
    ncut_grid       Corresponding ncut values.
    vol_grid        Corresponding volumes.
    subs_grid       Corresponding sizes of intersections with subsets.
 


## License

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


If you use this code for your publication, please include a reference 
to the paper "Constrained fractional set programs and their application in 
local clustering and community detection".
 


