# GiVSA
A novel stochastic search variable selection algorithm in normal linear regression problems, termed as group informed variable selection algorithm (GiVSA), which uses the known group structure efficiently to explore the model space without discarding any covariate based on an initial screening.

# Main Algorithm:
GiVSA uses a MCMC algorithm to perform variable selection. The code for this algorithm is written in mcmc.r and the function GiVSA can be called once this file is loaded. mcmc.r has a dependecy in the form of Hash.cpp hich contains optimized code written in C++ to reduce runtime. These to files are sufficient to run the algorithm.

## Inputs:  
The function GiVSA takes the following inputs. Most of the inputs are optionals and it is the practitioners discretion to use them or not. 

    1. y, X : The response and design matrix.  
    2. M : The paramter to control neighborhood size, default value is infinity which implies that the neighborhhod size only dpend on the size of the group (M plays no role on when assigned this default value).  
    3. N : The number of iterations of MCMC.  
    4. bknot : The burning period.  
    5. start : The starting model for the markov chain, default is the null model.  
    6. groups : The group structure in X.  
    7. min_group_size : Merges groups having size less than this parameter into a single group, default is to only merge singletons in a single group.  
    8. target_size : The target model size (not in use, as per the discussion and our theory left at 0)  
    9. delta : increase delta to decrease K_n.  
    10. c0 : increase c0 to increase g_n.  
    11. c : increase c to increase K_n.  
    l2. al : We consider a convex sum of GSIS and uniform probability, al is the weight given to GSIS (default 2/3, corresponding to 2:1)  

## Outputs: 
Returns a list of all the models visited by the markov chain after burning period, along with the model with highest posterior probability. 

## Analysis 
Output(groups).r : A script which takes the output of the mcmc algorithm and returns the mode and median of the chain, as well as the model with highest posterior probability.
