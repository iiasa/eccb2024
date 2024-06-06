
## packages
library(terra)
library(prioritizr)
library(foreign)
library(dplyr)
library(highs)
library(tidyr)
library(tidyverse)
library(rcbc)
library(viridisLite)

source("FUNCTIONS.R") ## other home made functions.

### this script reads input data, creates and solves problem, maps the output.

## ## ## ## ## ## ## ##
# Planning units ######
## ## ## ## ## ## ## ##

## planning units = raster of the study area (reference grid)
# If budget should be expressed as area (and not cost), in the PU raster, each grid cell should have equal value 1 (so that the budget can be expressed in terms of number of area/grid cells)
PU <- rast("data/PlanningUnits.tif")

## ## ## ## ## ## ## ##
# distributions #######
## ## ## ## ## ## ## ##

# Features = 67 tree species
spp.list <- list.files(path = "data/SpeciesDistributions/", full.names = T, recursive = T, pattern = "tif$")
spp <- rast(spp.list[grep("current", spp.list)])

# rename feature layers by species names
# this will enable to link the features rasters to the table of feature characteristics, weights, targets, taxonomy...
names(spp) <- gsub("_ens-sdms_cur2005_prob_pot", "" ,names(spp))

# define area budget (unit: grid cells)
budget.area <- round(0.3 * length(cells(PU))) ##  30 percent


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##            0. create basic problem            ####
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


p <- problem(PU, spp)%>%
  add_min_shortfall_objective(budget = budget.area)%>%
  add_relative_targets(targets = 1) %>%  ## target is 100% for all species distributions.
  add_cbc_solver()%>% ## best performance across open solvers.
  add_proportion_decisions()## typically faster/better than binary decisions.

# solve
s <- solve(p)
# plot map
plot(s)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##     1. modify targets: try log linear targets ####
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

p1 <- problem(PU, spp)%>%
  add_min_shortfall_objective(budget = budget.area)%>%
  add_loglinear_targets(10, 1, 10^4, 0.5) %>% ## features that have a range size of 10 grid cells: protect the whole range. Features that are distributed over 10^4 or more grid cells: protect 50% of their range
  add_cbc_solver()%>% ## best performance across open solvers.
  add_proportion_decisions()

s1 <- solve(p1)

# plot map
plot(s1)

# note: targets can be informed by a combination of range size and red list status (see Jung et al 2021)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## 2. add weights that reflect red list status #######
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


## read red list information
redlist.trees <- read.csv('data/species_red_list.csv')

## assign weight based on red list status
redlist.trees$weight <- ifelse(redlist.trees$Global == "Vulnerable" | redlist.trees$Europe == "Vulnerable", 4,
                                      ifelse(redlist.trees$Global == "Near Threatened" | redlist.trees$Europe == "Near Threatened", 2,
                                             ifelse(redlist.trees$Global == "Data Deficient" | redlist.trees$Europe == "Data Deficient", 2, 1)))

## must be in the same order as the features (spp) rasterstack
rownames(redlist.trees) <- redlist.trees$spp_name
redlist.trees <- redlist.trees[names(spp),]

p2 <- p1 %>%
  add_feature_weights(redlist.trees$weight)

s2 <- solve(p2)
plot(s2) ## in this case the solution does not change drastically, since only a few species have a higher weight.


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##      2. bis.  plan for future distributions ######
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# read sdm under climate scenario rcp 8.5
spp.rcp85 <- rast(spp.list[grep("rcp85", spp.list)])

# rename feature layers by species names
names(spp.rcp85) <- gsub("_ens-sdms_rcp85_fut2065_prob_pot.tif", "" ,names(spp.rcp85))

## create problem with future distributions as features:
p2_bis <- problem(PU, spp.rcp85)%>%
  add_min_shortfall_objective(budget = budget.area)%>%
  add_loglinear_targets(10, 1, 10^4, 0.5) %>%
  add_feature_weights(redlist.trees$weight) %>%
  add_cbc_solver()%>%
  add_proportion_decisions()  ## entire grid cells (planning units) will be selected in the solution rather than a proportion

s2_bis <- solve(p2_bis)


## average across the two above solutions, to see what sites emerge as top priorities for these species, in both current and future climate conditions.
mean_s_climate <- mean(s2, s2_bis)
plot(mean_s_climate, col = viridisLite::mako(n = 10, direction = -1))

## note: this ignores existing protected areas!



## ## ## ## ## ## ## ## ## ##
# 3. Add Protected areas ####
## ## ## ## ## ## ## ## ## ##

# protected areas data
PA <- rast("data/protectedareas.tif") ## to find priorities that complement and expand on all protected areas

### try with locked in constraints functionality
p3 <- p2 %>%
  add_locked_in_constraints(stPA)  ## this locks in cells that have non zero and non NA values = binary. This is not suitable for European PA at 10x10k resolution.

s3 <- solve (p3) ## budget cannot be met, because protected areas are present in more than 30% of all planning unit. It would require changing the PA layer to a binary layer with a threshold .




## try again with the manual bounded constraints functionality to incorporate the amount protected per PU. ###
## create manual bounded constraints dataframe with protected area coverage per planning unit##
## this is to include proportion of the planning unit that is protected ##
pa_constraints <- data.frame(pu = cells(PA), ### cell ID
                             lower = unname(PA[!is.na(PA)]), ### lower bound that needs to be included in the solution = proportion of grid cell already protected
                             upper = 1) ## upper bound = the whole planning unit can be selected

p3 <- p2 %>%
  add_manual_bounded_constraints(pa_constraints) ## to lock in proportional PA coverage per planning unit.


s3 <- solve(p3)



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##         3.bis : plan for future distributions ####
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# read sdm under climate scenario rcp 8.5
spp.rcp85 <- rast(spp.list[grep("rcp85", spp.list)])

# rename feature layers by species names
names(spp.rcp85) <- gsub("_ens-sdms_rcp85_fut2065_prob_pot.tif", "" ,names(spp.rcp85))

## create problem with future distributions as features:
p3_bis <- problem(PU, spp.rcp85)%>%
  add_min_shortfall_objective(budget = budget.area)%>%
  add_loglinear_targets(10, 1, 10^4, 0.5) %>%
  add_feature_weights(redlist.trees$weight) %>%
  add_manual_bounded_constraints(pa_constraints)%>% ## to lock in proportional PA coverage per planning units
  add_cbc_solver()%>%
  add_proportion_decisions()  ## entire grid cells (planning units) will be selected in the solution rather than a proportion

s3_bis <- solve(p3_bis)


## Question: what areas emerge as climatically resilient protected area expansion priorities for these 67 species?
## average across the two solutions that expand on protected areas with current and future distributions:
mean_s_climate <- mean(s3, s3_bis)
expansion_climate <- mean_s_climate - PA

plot(expansion_climate, col = viridisLite::mako(n = 10, direction = -1))



## add complexity to the problem ####


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
###            4. Locked out constraints         ####
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## create locked out constraints : for no-go areas or areas that should be left out of the solution.

## here we use layers of high-intensity forests and urban areas as a proxy for no-go areas
## to lock out the planning units that have over 50% of urban/periurban, or over 50% of high intensity forest (reason: in these high-intensity areas, conservation would likely conflict with economic interests).

## from Dou et al., 2021
## aggregated at 10x10 k and aligned with PU
HI.forest <- rast("data/HI_forest_prct.tif")
urban <- rast("data/urban_prct.tif")

locked.out <- sum(HI.forest, urban)
rclmat <- matrix(ncol = 3, nrow = 2, byrow = T,
                 c(0,50, 0,
                   50, 101, 1))

locked.out.bin <- terra::classify(locked.out, rclmat) ## convert to binary : 1 = pu that have more than 50% coverage of urban and/or HI forest. Else 0

p4 <- p3 %>%
  add_locked_out_constraints(locked.out.bin)
  ## important: this MUST come AFTER the manual bounded constraints (if using), otherwise locked out constraints are ignored.
  ## note that locked out constraints can sometimes also conflict with the manual bounded constraints, in other words locked in PA might become locked out...

s4 <- solve(p4)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#####             5. Linear penalties          ######
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## linear penalties can be used to avoid the selection of sites with a high value, for example, socio-economic costs if available.
## here, we use the human modification index as a penalty.

gHM <- rast("data/gHM.tif")

gHM[gHM<0.3] <- 0 ## set threshold so that sites that have GHM index lower than specified threshold are not penalized

p5 <- p4 %>%
  add_linear_penalties(penalty = 1, data = gHM) ## note that when penalty score is set too high, this sometimes prevents the budget area from being met

s5 <- solve (p5)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
###   6. Linear penalties with negative penalty  ####
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


## Linear penalties can also be used with a negative penalty score, to nudge the selection of sites with a high value.
## for example, one may use linear penalties with a negative penalty score to incorporate pre-defined connectivity layers; or known climate refugia...
## here, we use NDVI as an example. NDVI is often interpreted as dense and healthy vegetation, and one may be interested in selecting sites with a higher NDVI

ndvi <- rast("data/ndvi.tif")

p6 <- p5 %>%
  add_linear_penalties(penalty = -1, data = ndvi) ## note that when penalty score is set too high, this sometimes prevents the budget area from being met

s6 <- solve(p6)


### a word of caution: sometimes adding constraints and penalties will tend to drive the solution much more strongly than the biodiversity features themselves.
### example: try running the same problem with the future SDM instead of the current, and compare with the above solution
#
### create problem with future distributions as features:
# p6_bis <- problem(PU, spp.rcp85)%>%
#   add_min_shortfall_objective(budget = budget.area)%>%
#   add_loglinear_targets(10, 1, 10^4, 0.5) %>%
#   add_feature_weights(redlist.trees$weight) %>%
#   add_manual_bounded_constraints(pa_constraints)%>% ## to lock in proportional PA coverage per planning unit.
#   add_locked_out_constraints(locked.out.bin) %>%
#   add_linear_penalties(1, data = gHM) %>% ## note that when penalty score is set too high, this sometimes prevents the budget area from being met
#   add_linear_penalties(-1, data = ndvi) %>% ## negative penalty score can be used if we want to nudge selection of sites with high value in the spatial data layer
#   add_cbc_solver()%>%
#   add_proportion_decisions()  ## entire grid cells (planning units) will be selected in the solution rather than a proportion
#
# s6_bis <- solve(p6_bis)


### to limit the influence of the penalty data layer, you can consider decreasing the penalty value




## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#### 7. change decision variable to binary      #####
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## rewrite problem since we cannot overwrite the previously defined decision variable.
## need to use different constraints for Protected areas since the 30% budget cannot be met with binary decision + manually bounded constraints
PA_large <- PA
PA_large[PA_large<0.5] <- 0

pa_constraints_bin <- data.frame(pu = cells(PA_large), ### cell ID
                               lower = unname(PA_large[!is.na(PA_large)]), ### lower bound that needs to be included in the solution = proportion of grid cell already protected
                               upper = 1) ## upper bound = the whole planning unit can be selected

## create problem with binary decision:
p7 <- problem(PU, spp)%>%
  add_min_shortfall_objective(budget = budget.area)%>%
  add_loglinear_targets(10, 1, 10^4, 0.5) %>%
  add_feature_weights(redlist.trees$weight) %>%
  add_manual_bounded_constraints(pa_constraints_bin)%>% ## to lock in proportional PA coverage per planning unit.
  add_locked_out_constraints(locked.out.bin) %>%
  add_linear_penalties(1, data = gHM) %>% ## note that when penalty score is set too high, this sometimes prevents the budget area from being met
  add_linear_penalties(-1, data = ndvi) %>% ## negative penalty score can be used if we want to nudge selection of sites with high value in the spatial data layer
  add_cbc_solver()%>%
  add_binary_decisions()  ## entire grid cells (planning units) will be selected in the solution rather than a proportion

s7 <- solve(p7)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
####          8. Modify budget area to 10%      #####
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## modify the budget: e.g. 10% top priorities that expand on strict protected areas
budget.area <- round(0.1 * length(cells(PU))) ## 10 percent

stPA <-  rast("data/protectedareas_I_II.tif") ## to find priorities that complement and expand on strictly protected areas only (IUCN i and ii)

stpa_constraints <- data.frame(pu = cells(stPA), ### cell ID
                             lower = unname(stPA[!is.na(stPA)]), ### lower bound that needs to be included in the solution = proportion of grid cell already protected
                             upper = 1) ## upper bound = the whole planning unit can be selected

## create new problem for expansion of strict protected areas: new budget, new manual bounded constraints.
p8 <- problem(PU, spp)%>%
  add_min_shortfall_objective(budget = budget.area)%>%
  add_loglinear_targets(10, 1, 10^4, 0.5) %>%
  add_feature_weights(redlist.trees$weight) %>%
  add_manual_bounded_constraints(stpa_constraints)%>%
  add_locked_out_constraints(locked.out.bin) %>%
  add_linear_penalties(1, data = gHM) %>%
  add_linear_penalties(-1, data = ndvi) %>%
  add_cbc_solver()%>%
  add_proportion_decisions()

s8 <- solve(p8)



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#               II. Compare outputs           #######
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## 1. compare maps #####

# plot all solutions to compare them
plot(c(s, s1, s3, s3_bis, s4, s5, s6, s7, s8),
     main = c("basic problem", "add log linear targets",
              "add protected areas", "plan for future distributions" , "add locked-out constraints", "add gHM", "add NVDI",
              "binary decision", "change area budget to 10 percent strict PA expansion"),
     axes = FALSE)


##### map safe bets for expansion priorities across all variations of the problems expanding on protected areas ##
mean_s <- mean(s3, s3_bis, s4, s5, s6, s7)

exp <- mean_s - PA
# set negative values to zero (these correspond to planning units that were locked out by urban/forestry layer, but currently protected)
exp[exp<0]<-0
plot(exp, col = viridisLite::mako(n = 10, direction = -1))

plot(c(expansion_climate, exp), col = viridisLite::mako(n = 10, direction = -1))

## compare this map with the one obtained without considering additions in solutions 4-7 : how different are they?


## remember that the solutions are highly dependent on methodological choices, and specifically on the input data (features), constraints, and the objective function used, as well as the tool. Yet, the tool (zonation/prioritizr/marxan) leads to a smaller difference than the input data (features), the constraints, and the chosen planning objective
## also note that, the more features are included, the more the solution converges, hence it is always better to include as many high-quality features data as possible, to have an ecologically robust output.


## 2. compare performance #####

## compare the performance of different solutions in terms of representation:
## create rasterstack of solutions for which you want to compare performance
## here, we compare the scenarios that optimize for current distributions within 30% budget area
solutions <- c(s, s1, s2, s3, s4, s5, s6, s7)
names(solutions) <- c("basic problem", "add log linear targets", "add weights",
                      "add protected areas",  "add locked-out constraints", "add gHM", "add NVDI",
                      "binary decision" )

## analyse representation gains in the different solutions with a given budget
## for individual species
scenarios_performance_species <- data.frame(solution = character(),
                              feature = character(),
                              class = character(),
                              order = character(),
                              relative_held = numeric(),  ## representation: percentage of distribution held in the solution
                              relative_shortfall = numeric()) ## shortfall to target: how far from the area target for each species

## loop across solutions to extract representations for species and target shortfall
for (i in 1:nlyr(solutions)){
  cat(paste0(i, " \n"))  # keep track
  rpz.s_i <- eval_target_coverage_summary(p1, solutions[[i]]) ## for each species. Note that here, we assess the target shortfall based on the targets defined in p1, i.e. log linear targets.
  rpz.s_i$order <- redlist.trees$order[match(rpz.s_i$feature, redlist.trees$spp_name)]
  rpz.s_i$class <- redlist.trees$class[match(rpz.s_i$feature, redlist.trees$spp_name)]
  rpz.s_i$solution <- names(solutions)[i]
  rpz_i <- as.data.frame(rpz.s_i)
  scenarios_performance_species <- rbind(scenarios_performance_species,
                           rpz_i[, c("solution", "feature",  "class","order", "relative_held", "relative_shortfall")]
                           )
}


scenarios_performance_species$solution <- factor(scenarios_performance_species$solution, levels = names(solutions)) ## to plot in the right order.

## compare performance of different solutions in terms of representation
ggplot(scenarios_performance_species, aes(x = solution, y = relative_held)) +
  geom_boxplot()+
  theme_bw()

# subdivide per groups of species to be more ecologically informative
ggplot(scenarios_performance_species, aes(x = solution, y = relative_held)) +
  geom_boxplot(aes(fill = order), alpha = 0.2, outlier.size = 0)+
  theme_bw()

## add jitter points to see individual species representations
ggplot(scenarios_performance_species, aes(x = solution, y = relative_held)) +
  geom_boxplot(aes(fill = class), alpha = 0.2, outlier.size = 0)+
  geom_point(aes(x = solution, y = relative_held, colour = class), position = position_jitterdodge())+
  theme_bw()


## compare performance of different solutions in terms of target shortfall
ggplot(scenarios_performance_species, aes(x = solution, y = relative_shortfall)) +
  geom_boxplot(aes(fill = class), alpha = 0.2, outlier.size = 0)+
  geom_point(aes(x = solution, y = relative_shortfall, colour = class), position = position_jitterdodge())+
  theme_bw()

## what do these two performance metrics tell us?



## sometimes, one might be interested in the average representation as a barplot, aggregated for groups of species.
## to do so, we can use the home made function (see FUNCTIONS.R script) summarize_rpz


# evaluate representation under different scenarios for current distributions, grouped by taxonomic class
mean_rpz_solutions_by_class <- summarize_rpz(problem = p, solutions = solutions, feature_info_table = redlist.trees, group.by = "class")

ggplot(data = mean_rpz_solutions_by_class, aes(x = solutions, y = mean_rpz, fill = group.by)) +
  geom_bar(stat = "identity", position = "dodge")



## 3. Bonus: create a spatial ranking for a given problem   #####
##### produce a ranking map with increasing for loop budget #####
# for this, we will build on solution #3: that expands on protected areas for current distributions, but does not include other constraints
# we will start with the existing protected area and incrementally add budget until the whole study area is reached
# then, we can average across all solutions to obtain the ranking.

## initialise a raster stack with existing PA to store solutions as budget area increases.
incremental.solutions <- PA

protected.land <-  round(sum(PA[PA>0]))
total.land <- sum(PU[PU>0])

steps <- c(seq(from =protected.land, to = total.land, by = 5000 )[-1], total.land-1)

## skip the first as this is the initial PA layer + add the total land amount
## the argument "by" can be decreased for finer ranking.

## Note: this will take a while (1-2 minutes per run)
for (budget.area in steps){
  p_i <- problem(PU, spp)%>%
    add_min_shortfall_objective(budget = budget.area)%>%
    add_loglinear_targets(10, 1, 10^4, 0.5) %>%
    add_feature_weights(redlist.trees$weight) %>%
    add_manual_bounded_constraints(pa_constraints)%>%
    add_cbc_solver()%>%
    add_proportion_decisions()

  s_i <- solve(p_i)
  incremental.solutions <- c(incremental.solutions, s_i)
}


ranking.expansion.priorities <- mean(incremental.solutions) - PA

plot(ranking.expansion.priorities, col = viridisLite::magma(n = 100, direction =-1))









#### Other useful functionalities ####

## To use a different solver :
# add_gurobi_solver() %>%  ## requires installing gurobi first (free with academic email)
# add_highs_solver() %>% ## not the best open solver.
# add_cbc_solver()%>% ## best performance across open solvers.

