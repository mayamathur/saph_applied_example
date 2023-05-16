
# IMPORTANT NOTES -----------------------------

# This script is designed to be run in interactive Sherlock session,
#  NOT via an sbatch file. It analyzes ONE meta-analysis.
#  To analyze multiple metas, use the script 2022-3-11 applied doParallel SAPH.R

# To quickly run this script in high-mem interactive session:
# setwd("/home/groups/manishad/SAPH/applied_examples/code"); source("2_analyze_lodder_sherlock.R")


# This is very similar to doParallel.R.
# The only real additions/changes are:
#  - plot_trunc_densities_RTMA


# "show rep res"
# quickly look at results when running locally
srr = function() {
    cat("\n")
    print( rep.res %>%
             mutate_if(is.numeric, function(x) round(x,2)) )
    cat("\n")

}



# ~~ Load Packages -----------------------------------------------

toLoad = c("crayon",
           "dplyr",
           "foreach",
           "doParallel",
           "boot",
           "metafor",
           "robumeta",
           "data.table",
           "purrr",
           "metRology",
           "fansi",
           "MetaUtility",
           "testthat",
           "ICC",
           "cfdecomp",
           "tidyr",
           "truncdist",
           "tibble",
           "tmvtnorm",
           "testthat",
           "truncreg",
           "truncnorm",
           "rstan",
           "optimx",
           "phacking",
           "weightr",
           "here")



# load within installation if needed
for (pkg in toLoad) {
  
  cat( paste("\nAbout to try loading package", pkg) )
  
  tryCatch({
    # eval below needed because library() will otherwise be confused
    # https://www.mitchelloharawild.com/blog/loading-r-packages-in-a-loop/
    eval( bquote( library( .(pkg) ) ) )
  }, error = function(err) {
    install.packages(pkg)
  })
  
}


# ~~ User-Specified Parameters -----------------------------------------------

# which dataset to run?
# "Fig 1" or "Appendix A"
#dataset.to.run = "Appendix A" (NOT the one analyzed in RSM paper)
dataset.to.run = "Fig 1" 


data.dir = here("Data and materials/Prepped data")
local.results.dir = here("Results from R")

# specify which methods to run, as in doParallel
# but obviously can't run gold-std on a non-simulated meta-analysis
all.methods = "naive ; maon ; pcurve ; 2psm ; rtma ; prereg-naive"
# parse methods string
all.methods = unlist( strsplit( x = all.methods,
                                split = " ; " ) )


# load data
if ( dataset.to.run == "Fig 1" ) {
  setwd(data.dir)
  dp = read.csv("lodder_prepped.csv")
  expect_equal( nrow(dp), 287 )
}

# helper code
setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Linked to OSF (SAPH)/Code (git)/Sherlock code")
source("helper_SAPH.R")
 


# MAKE DATA SUBSETS ------------------------------

# published nonaffirmatives only
dpn = dp[ dp$affirm == FALSE, ]

# special dataset for CSM: 
# throws away affirmatives from hacked studies (i.e., all non-prereg studies)
dp.csm = dp %>% filter( Preregistered == TRUE | affirm == FALSE )

dp.csm %>% group_by(Preregistered, affirm) %>%
  summarise(n())

# unhacked only
dp.prereg = dp %>% filter(Preregistered == TRUE)
dpn.prereg = dpn %>% filter(Preregistered == TRUE)



# RUN ANALYSIS ------------------------------

# initialize rep.res st run_method_safe and other standalone estimation fns
#  will correctly recognize it as having 0 rows
rep.res = data.frame()

# # ~~ Start Values ------------------------------
# #@UNLIKE in doParallel, here we start at (0,1)
# #  by default, but if running method-naive, those will be the start values instead
# Mhat.start = 0
# Shat.start = 1
# 
# # in case we're not doing jeffreys-mcmc or it fails
# Mhat.MaxLP = NA
# Shat.MaxLP = NA
# 
# Mhat.MAP = NA
# Shat.MAP = NA

# ~~ Fit Naive Meta-Analysis (All PUBLISHED Draws) ------------------------------

if ( "naive" %in% all.methods ) {
  rep.res = run_method_safe(method.label = c("naive"),
                            method.fn = function() {
                              mod = rma( yi = dp$yi,
                                         vi = dp$vi,
                                         method = "REML",
                                         knha = TRUE )
                              
                              report_meta(mod, .mod.type = "rma")
                            },
                            .rep.res = rep.res )
  
  Mhat.naive = rep.res$Mhat[ rep.res$method == "naive" ]
  Shat.naive = rep.res$Shat[ rep.res$method == "naive" ]
}

srr()

# # ~~ Change Starting Values -----
# if ( !is.na(Mhat.naive) ) Mhat.start = Mhat.naive 
# if ( !is.na(Shat.naive) ) Shat.start = Shat.naive 



# ~~ Fit MAON (Nonaffirmative Published Draws) ------------------------------

if ( "maon" %in% all.methods ) {
  
  rep.res = run_method_safe(method.label = c("maon"),
                            method.fn = function() {
                              mod = robu( yi ~ 1, 
                                          data = dpn, 
                                          studynum = 1:nrow(dpn),
                                          var.eff.size = vi,
                                          small = TRUE )
                              
                              report_meta(mod, .mod.type = "robu")
                            },
                            .rep.res = rep.res )
  
}


srr()

# ~~ Fit 2PSM (All Published Draws) ------------------------------

if ( "2psm" %in% all.methods ) {
  
  rep.res = run_method_safe(method.label = c("2psm"),
                            method.fn = function() {
                              
                              naive = rma( yi = dp$yi,
                                                 vi = dp$vi,
                                                 method = "REML",
                                                 knha = TRUE )
                              
                              mod = metafor::selmodel(naive, type="stepfun", alternative = "greater", steps=c(0.025, 1))
                              
                              stats = report_meta(mod, .mod.type = "rma")
                              
                              stats$stats$Eta = 1/mod$delta[2]
                              stats$stats$EtaLo = 1/mod$ci.ub.delta[2]
                              stats$stats$EtaHi = 1/mod$ci.lb.delta[2]
                              
                              stats
                            },
                            .rep.res = rep.res )
  
}

srr()


# ~~ Fit P-Curve (Published Affirmatives) ------------------------------

if ( "pcurve" %in% all.methods ) {
  # since pcurve.opt internally retains only 
  rep.res = run_method_safe(method.label = c("pcurve"),
                            method.fn = function() {
                              #@later, revisit the decision to use df_obs = 1000
                              #   to effectively treat yi/sei z-scores
                              
                              # pass all significant studies (that's what pcurve.opt does internally):
                              Mhat = pcurve.opt( t_obs = dp$Zi,
                                                 df_obs = rep(1000, length(dp$Zi)),
                                                 dmin = -5, #@HARD-CODED and arbitrary
                                                 dmax = 5)
                              
                              
                              # # using only published affirmatives:
                              # Mhat = pcurve.opt( t_obs = dpa$yi/dpa$sei,
                              #                    df_obs = rep(1000, length(dpa$yi)),
                              #                    dmin = -5, #@HARD-CODED and arbitrary
                              #                    dmax = 5)
                              
                              
                              return( list( stats = data.frame( Mhat = Mhat) ) )
                            },
                            .rep.res = rep.res )
  
}

# ~~ RTMA ------------------------------

if ( "rtma" %in% all.methods ) {
  # # temp for refreshing code
  # path = "/home/groups/manishad/SAPH"
  # setwd(path)
  # source("helper_SAPH.R")
  
  # this one has two labels in method arg because a single call to estimate_jeffreys_mcmc
  #  returns 2 lines of output, one for posterior mean and one for posterior median
  # order of labels in method arg needs to match return structure of estimate_jeffreys_mcmc
  rep.res = run_method_safe(method.label = "rtma",
                            method.fn = function() {
                              
                              mod = phacking_meta(yi = dp$yi,
                                                  sei = sqrt(dp$vi))
                              
                              # follow the same return structure as report_meta
                              list( stats = data.frame( Mhat = mod$stats$mode[1],
                                                        MLo = mod$stats$ci_lower[1],
                                                        MHi = mod$stats$ci_upper[1],
                                                        
                                                        Shat = mod$stats$mode[2],
                                                        SLo = mod$stats$ci_lower[2],
                                                        SHi = mod$stats$ci_upper[2] ) ) 
                              
                              
                            },
                            .rep.res = rep.res )

  cat("\n doParallel flag: Done rtma if applicable")
}

cat("\n")
rep.res
cat("\n")


# ANALYSES OF PREREG STUDIES (SPECIFIC TO LODDER) ------------------

# ~~ Fit Naive (Prereg Only) ------------------------------

if ( "prereg-naive" %in% all.methods ) {
  rep.res = run_method_safe(method.label = c("prereg-naive"),
                            method.fn = function() {
                              mod = rma( yi = dp.prereg$yi,
                                         vi = dp.prereg$vi,
                                         method = "REML",
                                         knha = TRUE )
                              
                              report_meta(mod, .mod.type = "rma")
                            },
                            .rep.res = rep.res )
  
  Mhat.naive = rep.res$Mhat[ rep.res$method == "prereg-naive" ]
  Shat.naive = rep.res$Shat[ rep.res$method == "prereg-naive" ]
}

srr()



# DATASET INFO AND SANITY CHECKS ------------------------------

# Info About Dataset ------------

rep.res$k.pub = nrow(dp)
rep.res$k.pub.affirm = sum(dp$affirm == TRUE)
rep.res$k.pub.nonaffirm = sum(dp$affirm == FALSE)
rep.res$prob.pub.study.affirm = rep.res$k.pub.affirm / rep.res$k.pub


# ~ WRITE LONG RESULTS ------------------------------

setwd(local.results.dir)
fwrite( rep.res, "results_lodder.csv" )

