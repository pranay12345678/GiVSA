# GiVSA
A novel stochastic search variable selection algorithm in normal linear regression problems, termed as group informed variable selection algorithm (GiVSA), which uses the known group structure efficiently to explore the model space without discarding any covariate based on an initial screening.

# Markov Chain:
mcmc.r is our gGiVSA algorithm.  
Hash.cpp is just to reduce runtime for the mcmc algorithm.   
Inputs:  
    a. y, X : The response and design matrix.  
    b. M : The paramter to control neighborhood size, default value is infinity which implies that the neighborhhod size only dpend on the size of the group (M plays no role on when assigned this default value).  
    c. N : The number of iterations of MCMC.  
    d. bknot : The burning period.  
    e. start : The starting model for the markov chain, default is the null model.  
    f. groups : The group structure in X.  
    g. min_group_size : Merges groups having size less than this parameter into a single group, default is to only merge singletons in a single group.  
    h. target_size : The target model size (not in use, as per the discussion and our theory left at 0)  
    i. delta : increase delta to decrease K_n.  
    j. c0 : increase c0 to increase g_n.  
    k. c : increase c to increase K_n.  
    l. al : We consider a convex sum of GSIS and uniform probability, al is the weight given to GSIS (default 2/3, corresponding to 2:1)  

Outputs: Returns a list of all the models visited by the markov chain after burning, along with the model with highest posterior probability.  

# Analysis 
Output(groups).r : A script which takes the output of the mcmc algorithm and returns the mode and median of the chain, as well as the model with highest posterior probability.
