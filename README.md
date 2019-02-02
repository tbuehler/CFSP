# Constrained fractional set programming



This archive contains a Matlab implementation of the CFSP algorithm for 
community detection as described in the paper
 
T. Buehler, S. S. Rangapuram, S. Setzer and M. Hein
Constrained fractional set programs and their application in 
local clustering and community detection
ICML 2013, pages 624-632
(Extended version available at http://arxiv.org/abs/1306.3409)

Current version: V1.0



## Installation

The implementation uses a mex-file to solve the inner problem of the
RatioDCA. Compile by typing 'make' in the matlab command line.



## Documentation

For more information, type 'help functionname' in the Matlab command line.

### Community detection via constrained maximum density:

#### Usage: 

    [dc_best, maxdens_best, gam_best, dc_all,maxdens_all,gam_all] ... 
        = community_detection(W, k1,k2,gdeg,seed,numRuns,dist_max)

#### Input: 

    W         The weight matrix.
    k1        Lower bound
    k2        Upper bound
    seed      Index of seed set

#### Optional Input:

    gdeg      Generalized degrees. Default is all ones vector.
    numRuns   Number of runs. Default is 10.
    dist_max  Maximum distance from the seed vertices. Default is 2.

#### Output:

    dc_best       Indicator vector of the community with largest density
    maxdens_best  The corresponding density
    gam_best      The gamma value used in the optimization
    dc_all        Cell array of indicator vector for each initialization
    maxdens_all   The densities of the corresponding communities
    gam_all       The corresponding gamma values



### Constrained Ncut, direct integration of seed subset, volume constraint via penalty term:

#### Usage:

    [clusters, ncut,feasible,lambda]= vol_cnstr_ncut_subset_direct(W,k,hdeg,start,subset,gamma)

#### Input:
    W                 Weight matrix (full graph).
    k                 Upper bound
    hdeg              Generalized degrees used in constraint.
    start             Start vector
    subset            Indices of seed subset
    gamma             Penalty parameter for volume constraint.

#### Output:

    clusters          Thresholded vector f yielding the best objective.
    ncut              Ncut value of resulting clustering.
    feasible          True if all constraints are fulfilled.
    lambda            Corresponding objective.


### Constrained Ncut, seed constraint and volume constraint via penalty term

#### Usage: 
    
    [clusters, ncut,feasible,lambda]= vol_cnstr_ncut_subset_penalty(W,k,hdeg,start,subset,gamma1,gamma2)

#### Input:
    
    W                 Weight matrix (full graph).
    k                 Upper bound
    hdeg              Generalized degrees used in constraint.
    start             Start vector
    subset            Indices of seed subset
    gamma1             Penalty parameter for seed constraint.
    gamma2             Penalty parameter for volume constraint.

#### Output:

    clusters          Thresholded vector f yielding the best objective.
    ncut              Ncut value of resulting clustering.
    feasible          True if all constraints are fulfilled.
    lambda            Corresponding objective.
 



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
 


