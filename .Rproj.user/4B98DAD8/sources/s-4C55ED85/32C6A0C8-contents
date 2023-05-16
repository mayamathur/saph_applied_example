
# PRELIMINARIES ---------------------------------------------


# you just saved results from 2_analyze_lodder_sherlock (but should get integrated into this file)
# next clean up this file and have it write the etahats, taus, etc.
# YOU GOT THIS!


toLoad = c("crayon", "dplyr", "foreach", "doParallel", "boot", "metafor", 
            "robumeta", "data.table", "purrr", "metRology", "fansi", "MetaUtility", 
            "ICC", "cfdecomp", "tidyr", "truncdist", "tibble", "tmvtnorm", 
            "testthat", "truncreg", "truncnorm", "rstan", "optimx", "weightr", 
            "here", "RColorBrewer", "phacking")

lapply( toLoad,
        require,
        character.only = TRUE)

# helper fns
general.code.dir = "~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Linked to OSF (SAPH)/Code (git)/Sherlock code"
setwd(general.code.dir)
source("helper_SAPH.R")
source("analyze_sims_helper_SAPH.R")

# get prepped data
setwd( here("Lodder/Data and materials/Prepped data") )
dp = fread("lodder_prepped.csv")
expect_equal( nrow(dp), 287 )  # this is the "Figure 1" analysis, not the "Appendix A" analysis 


# for this script's own results
results.dir = here("Lodder/Results from R")
overleaf.dir.figs = "/Users/mmathur/Dropbox/Apps/Overleaf/P-hacking (SAPH) Overleaf/figures_SAPH/lodder"
overleaf.dir.nums = "/Users/mmathur/Dropbox/Apps/Overleaf/P-hacking (SAPH) Overleaf/results_from_R_SAPH"


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


# specify which methods to run, as in doParallel
# but obviously can't run gold-std on a non-simulated meta-analysis
all.methods = "naive ; maon ; pcurve ; 2psm ; rtma ; prereg-naive"
# parse methods string
all.methods = unlist( strsplit( x = all.methods,
                                split = " ; " ) )



# FIT MODELS  ------------------------------

# initialize rep.res st run_method_safe and other standalone estimation fns
#  will correctly recognize it as having 0 rows
rep.res = data.frame()

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



# ADD BASIC DATASET INFO TO RESULTS ------------

rep.res$k.pub = nrow(dp)
rep.res$k.pub.affirm = sum(dp$affirm == TRUE)
rep.res$k.pub.nonaffirm = sum(dp$affirm == FALSE)
rep.res$prob.pub.study.affirm = rep.res$k.pub.affirm / rep.res$k.pub


update_result_csv( name = "Lodder k total",
                   value = nrow(dp),
                   .overleaf.dir = overleaf.dir.nums )

update_result_csv( name = "Lodder k prereg",
                   value = sum(dp$Preregistered == 1),
                   .overleaf.dir = overleaf.dir.nums )

update_result_csv( name = "Lodder k non-prereg",
                   value = sum(dp$Preregistered == 0),
                   .overleaf.dir = overleaf.dir.nums )

update_result_csv( name = "Lodder k prereg affirm",
                   value = sum(dp$Preregistered == 1 & dp$affirm == 1),
                   .overleaf.dir = overleaf.dir.nums )

update_result_csv( name = "Lodder k prereg nonaffirm",
                   value = sum(dp$Preregistered == 1 & dp$affirm == 0),
                   .overleaf.dir = overleaf.dir.nums )

update_result_csv( name = "Lodder k non-prereg affirm",
                   value = sum(dp$Preregistered == 0 & dp$affirm == 1),
                   .overleaf.dir = overleaf.dir.nums )

update_result_csv( name = "Lodder k non-prereg nonaffirm",
                   value = sum(dp$Preregistered == 0 & dp$affirm == 0),
                   .overleaf.dir = overleaf.dir.nums )

update_result_csv( name = "Lodder k affirm",
                   value = sum(dp$affirm == 1),
                   .overleaf.dir = overleaf.dir.nums )

update_result_csv( name = "Lodder k nonaffirm",
                   value = sum(dp$affirm == 0),
                   .overleaf.dir = overleaf.dir.nums )


# ~ WRITE LONG RESULTS ------------------------------

setwd(results.dir)
fwrite( rep.res, "results_lodder.csv" )




# Z-SCORE DENSITY  ------------------------------


# these are WITHIN-study Z-scores
xmin = floor(min(dp$Zi))
xmax = ceiling(max(dp$Zi))

p = ggplot(data = data.frame(x = c(xmin, 3)),
            aes(x)) +
  
  geom_vline(xintercept = 0,
             lwd = 1,
             color = "gray") +
  
  geom_vline(xintercept = qnorm(0.975),
             lty = 2,
             color = "red") +
  
  # estimated density of estimates
  geom_density( data = dp,
                aes(x = Zi),
                adjust = .3 ) +
  
  ylab("") +
  scale_x_continuous( breaks = seq(-2, 8, 1),
                      limits = c(-2, 8)) +
  
  # put the Z=1.96 label slightly off from that location so it's visible with the line
  annotate(geom = "text", x = 1.75, y = 0.08, label = "Z = 1.96", color = "red",
           angle = 90) +
  
  xlab("Z-score") +
  theme_minimal() +
  scale_y_continuous(breaks = NULL) +
  theme( text = element_text(size = 16, face = "bold"),
         panel.grid.major.y = element_blank(),
         panel.grid.minor.y = element_blank() )



setwd(results.dir)
my_ggsave( name = paste( "lodder_z_density.pdf", sep="_" ),
           .width = 7, 
           .height = 6,
           .results.dir = results.dir,
           .overleaf.dir = overleaf.dir.figs )


# Z-SCORE DENSITIES STRATIFIED BY PREREG STATUS ------------------------------

dp$Preregistered[ dp$Preregistered == 1 ] = "Preregistered"
dp$Preregistered[ dp$Preregistered == 0 ] = "Not preregistered"


# these are WITHIN-study Z-scores
xmin = floor(min(dp$Zi))
xmax = ceiling(max(dp$Zi))

p = ggplot(data = data.frame(x = c(xmin, 3)),
           aes(x,
              color = Preregistered,
              fill = Preregistered) ) +

  scale_color_manual( values = c("black", "orange"),
                      name = "") +
  scale_fill_manual( values = c("black", "orange"),
                     name = "" ) +
  
  geom_vline(xintercept = 0,
             lwd = 1,
             color = "gray") +
  
  geom_vline(xintercept = qnorm(0.975),
             lty = 2,
             color = "red") +
  
  # estimated density of estimates
  geom_density( data = dp,
                aes(x = Zi),
                adjust = .3,
                alpha = 0.3
                ) +
  
  ylab("") +
  scale_x_continuous( breaks = seq(-2, 8, 1),
                      limits = c(-2, 8)) +
  
  # put the Z=1.96 label slightly off from that location so it's visible with the line
  annotate(geom = "text", x = 1.75, y = 0.3, label = "Z = 1.96", color = "red",
           angle = 90) +
  
  xlab("Z-score") +
  theme_minimal() +
  theme(legend.position="bottom") +
  scale_y_continuous(breaks = NULL) +
  theme( text = element_text(size = 16, face = "bold"),
         panel.grid.major.y = element_blank(),
         panel.grid.minor.y = element_blank() )


setwd(results.dir)
my_ggsave( name = paste( "lodder_z_density_by_prereg.pdf", sep="_" ),
           .width = 7, 
           .height = 6,
           .results.dir = results.dir,
           .overleaf.dir = overleaf.dir.figs )



# QQ PLOT FOR RTMA FIT  ------------------

Mhat = rep.res$Mhat[ rep.res$method == "rtma" ]
Shat = rep.res$Shat[ rep.res$method == "rtma" ]

p = yi_qqplot(yi = dpn$yi,
              sei = dpn$sei,
              Mhat = Mhat,
              Shat = Shat) + 
  theme_bw() +
  theme( text = element_text(size = 16, face = "bold"),
         panel.grid.major.y = element_blank(),
         panel.grid.minor.y = element_blank() )

setwd(results.dir)
my_ggsave( name = paste( "lodder_rtma_qqplot.pdf", sep="_" ),
           .width = 7, 
           .height = 7,
           .results.dir = results.dir,
           .overleaf.dir = overleaf.dir.figs )

# FOREST PLOT OF ESTIMATES BY METHOD ------------------------------

# ~ Recode method ------------------------

# only name the methods we're going to show in paper
rep.res$method.pretty = NA
rep.res$method.pretty[ rep.res$method == c("naive") ] = "Uncorrected"
rep.res$method.pretty[ rep.res$method == c("maon") ] = "MAN"
rep.res$method.pretty[ rep.res$method == c("2psm") ] = "SM"
rep.res$method.pretty[ rep.res$method == c("prereg-naive") ] = "Preregistered only"
rep.res$method.pretty[ rep.res$method %in% c("rtma") ] = "RTMA"
table(rep.res$method, rep.res$method.pretty)



# ~ Plotting dataframe ----------------------
rsp = rep.res %>% filter( !is.na(method.pretty) &
                       # I accidentally included 2 rows for naive
                       !duplicated(method.pretty) ) %>%
  droplevels()
table(rsp$method.pretty)


# force ordering of methods
correct.order = c("Uncorrected",
                  "SM",
                  "Preregistered only",
                  
                  "RTMA",
                  "MAN")

rsp$method.pretty = factor(rsp$method.pretty, levels = rev(correct.order))
levels(rsp$method.pretty)

# ~ Make plot ----------------------

my.shapes = c(16, 2)

# SAVE: this is to set colors automatically
# # have colors and shapes match sim study
# # taken from inside analyze_sims_helper::sim_plot_multiple_outcomes
# n.colors.needed = length(unique(rsp$method.pretty))
# .colors = brewer.pal(n = n.colors.needed, name = "Dark2")
# if( length(.colors) > n.colors.needed ) .colors = .colors[1:n.colors.needed]
# # this maps the colors onto levels of the factor
# names(.colors) = levels( factor(rsp$method.pretty) )
# 
# # highlight certain methods
# .colors[ names(.colors) == "RTMA" ] = "red"

# hard-code colors to match simulations
.colors = c(#SMKH = "#1B9E77",
            MAN = "#ff9900",
            RTMA = "red",
            `Preregistered only` = "#3399ff",
            SM = "#00cc00",
            Uncorrected = "black")


# ~~ Set ggplot linetype scale ----
# by default, dotted lines
# but use solid lines for new proposed methods
.lty = rep("dashed", nuni(rsp$method.pretty))
names(.lty) = names(.colors)

newMethods = c("MAN",
               "RTMA")

.lty[ names(.lty) %in% newMethods ] = "solid"



# ~ Make plot ----------------------

p = ggplot( data = rsp,
            aes( x = Mhat,
                y = method.pretty, 
                color = method.pretty,
                shape = (method.pretty == "Uncorrected"),
                linetype = method.pretty ) ) +
  
  geom_vline(xintercept = 0,
             lwd = .8,
             color = "gray") +
  
  geom_point(size = 3) +
  geom_errorbarh( aes(xmax = MHi,
                      xmin = MLo),
                  height = 0,
                  lwd = 0.8 ) +
  
  # manually provided colors
scale_colour_manual(values = .colors ) +
  
  # manually provided linetypes
  scale_linetype_manual(values = .lty,
                        guide = "none") +

  scale_shape_manual(values = my.shapes,
                     guide = "none") +
  
  #coord_cartesian( xlim = c(xmin, xmax) )+
  scale_x_continuous(breaks= seq(-0.05, 0.35, 0.05 ) ) +
  
  xlab( "Pooled estimate with 95% CI" ) +
  ylab("") +
  labs(color  = "Method", linetype = "Method", shape = "Method") +
  # manually set shapes in color legend to match the second legend
  guides(colour = guide_legend( override.aes = list( shape = c(2,16,16,16,16) ),
                                reverse = TRUE) ) +
  
  theme_bw() +

  theme( text = element_text(face = "bold"),
         panel.grid.major.y = element_blank(),
         panel.grid.minor.y = element_blank(),
         legend.position = "right" )
  

setwd(results.dir)
my_ggsave( name = paste( "lodder_forest.pdf", sep="_" ),
           .width = 7, 
           .height = 3,
           .results.dir = results.dir,
           .overleaf.dir = overleaf.dir.figs )





# WRITE INDIVIDUAL NUMERICAL ESTIMATES ------------------

# Mhat and CI limits
update_result_csv( name = paste( "Lodder", rsp$method.pretty, "Mhat", sep = " "),
                   value = round( rsp$Mhat, 2),
                   .overleaf.dir = overleaf.dir.nums )

update_result_csv( name = paste( "Lodder", rsp$method.pretty, "MLo", sep = " "),
                   value = round( rsp$MLo, 2),
                   .overleaf.dir = overleaf.dir.nums )


update_result_csv( name = paste( "Lodder", rsp$method.pretty, "MHi", sep = " "),
                   value = round( rsp$MHi, 2),
                   .overleaf.dir = overleaf.dir.nums )

# Shat and CI limits
# note: MAN has no CI for Shat because robu doesn't provide one
update_result_csv( name = paste( "Lodder", rsp$method.pretty, "Shat", sep = " "),
                   value = round( rsp$Shat, 2),
                   .overleaf.dir = overleaf.dir.nums )

update_result_csv( name = paste( "Lodder", rsp$method.pretty, "SLo", sep = " "),
                   value = round( rsp$SLo, 2),
                   .overleaf.dir = overleaf.dir.nums )

update_result_csv( name = paste( "Lodder", rsp$method.pretty, "SHi", sep = " "),
                   value = round( rsp$SHi, 2),
                   .overleaf.dir = overleaf.dir.nums )

# EtaHat for SM
temp = rsp %>% filter(method == "2psm")
update_result_csv( name = paste( "Lodder", temp$method.pretty, "EtaHat", sep = " "),
                   value = round( temp$Eta, 2),
                   .overleaf.dir = overleaf.dir.nums )

update_result_csv( name = paste( "Lodder", temp$method.pretty, "EtaLo", sep = " "),
                   value = round( temp$EtaLo, 2),
                   .overleaf.dir = overleaf.dir.nums )

update_result_csv( name = paste( "Lodder", temp$method.pretty, "EtaHi", sep = " "),
                   value = round( temp$EtaHi, 2),
                   .overleaf.dir = overleaf.dir.nums )


