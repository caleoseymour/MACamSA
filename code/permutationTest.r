## Permutation test!
## Test data!
#crickets = read.csv('crickets.csv', header=TRUE)

permute_differences = function(seed, x, groups, seeded=FALSE)
{
    if (seeded)
    {
        saved.seed = .Random.seed
        set.seed(seed)
    }
    sample1 = sample(x, sum(groups == unique(groups)[1]),replace=F)
    notsample1 = x[!x %in% sample1]
    
    if(seeded)
    {
        set.seed(seed)
    }
    sample2 = sample(x, sum(groups == unique(groups)[2]),replace=F)
    notsample2 = x[!x %in% sample2]
    out = c(mean(sample1) - mean(notsample1), mean(sample2) - mean(notsample2))
    
    names(out) = unique(groups)
    
    if(seeded)
    {
        set.seed(saved.seed)
    }
    return(out)    
}

permutation.test = function(x, groups, tails=2, nperm=10000)
{
# Function permutation.test: Run permutation test on data in groups.
# arguments: 

    seeded = (exists('.Random.seed'))
    means=tapply(x,groups,mean)
    n=tapply(x,groups,length)
    combinations = combn(1:length(means),m=2)
    
    out = matrix(ncol=2,nrow=0)
    
    ## Small for loop; number of comparisons. It shouldn't have a major effect
    ## on runtime unless there are thousands of group variables.
    for (i in 1:ncol(combinations)) {
        current_comparison = unique(groups)[combinations[,i]]
        diffmeans = means[combinations[1,i]] - means[combinations[2,i]]
        
        ## Randomly generate nperm random seeds to plug into runif.
        ## It's reproducible AND pseudorandom, at the very least.
        if(seeded)
        {
            set.seed(.Random.seed)
        }
        seeds = ceiling(runif(nperm, -1, .Machine$integer.max))
        
        differences = vapply(X=seeds, FUN=permute_differences, FUN.VALUE=double(2), x=x[groups %in% current_comparison], groups = groups[groups %in% current_comparison], seeded=seeded)
        
        critical_value = (rowSums(differences>=diffmeans)/nperm)
        
        p.values = (0.5 * (diffmeans == 0)) + (critical_value * (diffmeans>0)) + ((1-critical_value) * (diffmeans<0))
        
        names(p.values) = paste(current_comparison, 'is a subset of', paste(current_comparison, collapse='+'))
        
        out = rbind(out,p.values)
    }
    
    ## Remove the random seed that was applied through the function earlier, if
    ## we didn't start out with one.
    if (!seeded)
    {
        rm(.Random.seed, envir = globalenv())
    }
    return(out)
}


