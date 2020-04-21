simulate_data = function(n_species=1, n_abiotic_vars=4, n_timepoints=5, n_plots=5, n_subplots_per_plot=20){
    
    #  ========================================
    #   setup fake dataframe 
    #  ========================================
    fake_data = data.frame(matrix(nrow=n_species*n_timepoints*n_plots*n_subplots_per_plot,ncol=(6+n_abiotic_vars)))
    ab_vars = rep("", n_abiotic_vars)
    for (i in 1:n_abiotic_vars){
        ab_vars[i] = paste("v", i, sep="")
    }
    colnames(fake_data) = c("time", "species", "cover", "presence", "plt", "subplt", ab_vars)
    
    
    #  ========================================
    #   simulate data
    #  ========================================
    
    
    # what is happening here? 
    #
    #  we will model the mean true pres/abs  
    #  
    #  strength of effect, beta_v for each env variable
    #  beta_v from normal(0, sigma_str) where sigma_str is size of effect variable
    #
    #  no interactions, but covariance among environmental variables 
    #
    #  how much covariance in conditions across plots? none here 
    #
    #  
    #  alpha = \sum_v (beta_v * e_vi) 
    #  logit(p_i) = alpha
    #  beta_v ~ normal(0, sigma_effect_str)
    #  ef_covariance_str ~ normal(0, sigma_covariance_str)
    #  v_i ~ mvnorm(n_vars, cov_matrix=diag(ef_covariance_str))
    #  
    #  Occupancy_i(t) ~ bernoulli(p_i) 
    
    
    environmental_factor_variance = 0.5
    environmental_factor_covariance  = 0.4
    effect_variance = 1.0
    ef_variance_over_time = 0.4
    ef_variance_across_plots = 5.0
    ef_variance_across_subplots = 1.0

    ## Environmental Factors
    covariance_mat = matrix(rnorm(n_abiotic_vars*n_abiotic_vars,sd=environmental_factor_covariance), nrow=n_abiotic_vars, ncol=n_abiotic_vars)
    covariance_mat = covariance_mat %*% t(covariance_mat)
    
    for (i in 1:n_abiotic_vars){
        covariance_mat[i,i] = rnorm(1, sd=environmental_factor_variance)
    }
    
    efs = rmvnorm(n_plots, sigma =covariance_mat)
    
    # Effect Strengths 
    beta = rnorm(n_abiotic_vars, sd=effect_variance)
    alpha = rnorm(1, mean=0.4, sd=0.5)
    
    
    row = 1
    for (species in 1:n_species){            
        for (plt in 1:n_plots){  
            mean_ef_this_plot = efs[plt] + rnorm(n_abiotic_vars, sd=ef_variance_across_plots)
            for (subplt in 1:n_subplots_per_plot){
                for (t in 1:n_timepoints){
                    
                    # additive noise for subplot and time
                    ef_vals_in_the_here_and_now = mean_ef_this_plot + rnorm(n_abiotic_vars, sd=ef_variance_over_time) 
                    
                    v = sum(diag(ef_vals_in_the_here_and_now) %*% beta) + alpha
        
                    # compute inverse logit to get probability of occurence
                    p_occurance = exp(v) / (exp(v) + 1)
                    
                    presence = rbinom(1,n=1, p=p_occurance)
                    
                    # obv this is not good but im not using it right now
                    abund = runif(1)
                    
                    # colnames(fake_data) ("time", "species", "cover", "presence", "plot", "subplot", ab_vars)
                    this_row = c(t, species, abund, presence, plt, subplt, ef_vals_in_the_here_and_now)
                    fake_data[row,] = this_row
                    row = row + 1
                }
            }
        }
    }
        
    fake_data$presence = as.logical(fake_data$presence)
    fake_data$plt = as.integer(fake_data$plt)

    return(list(df=fake_data, true_values=list(alpha=alpha, beta=beta, ef_covariance_mat = covariance_mat)))
}