
This repository establishes code for three main purposes: estimating linear
effects of covariates within the framework of a multivariate state-space model,
estimating cognitive factors scores within the framework of a factor state-space 
model, and simulating a randomized controlled trial application of cognitive 
factor scores to demonstrate the improvements in power in a hypothetical 
Alzheimer's disease drug trial using factor models.

Below is a general outline of how to use the provided code to get results from
the randomized controlled trial simulation.

Step 1: Estimate linear effects of covariates (to be removed in factor 
        analysis).
      
      Top level program: BayesKalmJoint_gh.R
                            Sub-program: BayesKalmUneq_gh.R
                                  Sub-sub-program: ffbsUneq_gh.R
                            Sub-program: ffbsJoint_gh.R

Step 2: Subtract linear effects from outcomes and re-standardize the resulting
        adjusted outcome values.
      
Step 3: Estimate cognitive factor scores from adjusted outcomes in step 2.
    
      Top level program: SSM_fac_gh.R
                            Sub-program: ffbsfac_gh.R
                            
Step 4: Using estimated factor scores in step 3, run RCT_simulation_gh.R to 
        simulate a set of randomized controlled clinical trials with the
        cognitive factor scores used in a predictive model to improve power
        of the hypothetical trials.
        
Step 5: Run tables_figures_gh.R to create tables and figures reflecting results
        from the simulation.
        
Training and testing data used in steps 4 and 5 are provided. Power, sample size, hazard ratio, and needed event data are provided for creating tables and figures in tables_figures_gh.R.

This R project incorporates the 'renv' package to maximize reproducibility. After opening the R project run renv::restore() when prompted to download and install the correct versions of all required packages. For more information on renv please see: https://rstudio.github.io/renv/articles/renv.html.

This code was written on R version 4.2.3 and may not work using other versions of R. For information on installing other versions of R: https://github.com/r-lib/rig.
                            

      
      
      
      
      
      
      
      