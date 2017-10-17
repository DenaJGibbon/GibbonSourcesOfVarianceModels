### R Code for Clink et al. A multivariate analysis of female Bornean gibbon calls reveals substantial inter-individual variation and site-level variation in trill rate 
### MANOVA Model Code 
### Load required libraries 
library(corrplot)
library(bayesplot)
library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores())
library(stringr)

## Load data 
lda.data <- read.csv("/Users/denajaneclink/STAN MANOVA/Clink.et.al.feature.data.csv")
str(lda.data)

## Isolate relevant features from data set
d.manova <- lda.data[,c("Delta.Time..s..1", "Freq.95...Hz..1","Delta.Time..s..2", 
	"Freq.95...Hz..2","Delta.Time..s..3", "Freq.95...Hz..3","Delta.Time..s..4", 
	"Freq.95...Hz..4","Delta.Time..s..5", "Freq.95...Hz..5", 
	"rest.dur.comb", "intro.dur", "trill.dur", "unlist.trill.rate.list.")]

## Assign more informative names to features 
names(d.manova) <- c("note.1.dur", "note.1.f95", "note.2.dur", "note.2.f95", 
	"note.3.dur", "note.3.f95", "note.4.dur", "note.4.f95", "note.5.dur", "note.5.f95",
	"rest.dur", "intro.dur", "trill.dur", "trill.rate") 

## Check the structure of the data
str(d.manova)

## Log transform data
d.manova <- log(d.manova)

## Isolate features to use for modelling
d.manova <- d.manova[, c("note.1.dur", "note.1.f95","note.2.dur", "note.2.f95",
	"rest.dur", "intro.dur", "trill.dur", "trill.rate")] 

## Check the structure of the data
str(d.manova)


### Set-up data to pass to Stan. 
# Integer-coded vector of group IDs
group.int <- as.numeric(lda.data$group)

# Check structure 
table(group.int)

# Integer-coded vector of site IDs
site.int <- as.numeric(as.factor(lda.data$site))

# Check structure 
table(site.int)

# Center data matrix at feature means
col.means <- apply(d.manova, MARGIN=2, FUN="mean")
y.centered <- sweep(d.manova, MARGIN=2, STATS=col.means)

# Create a data list to pass to Stan
data_list <- list(
  K = dim(d.manova)[2],
  J= length(unique(lda.data$group)),
  M= length(unique(lda.data$site)),
  N= dim(d.manova)[1],
  y= as.matrix(y.centered), ## features centered at zero
  group= group.int,
  site= site.int
)


# # Code to run the STAN mdoel
# # NOTE: the .stan file must be linked in the code below
# mfinal.stan = stan(file="/Users/denajaneclink/Downloads/MANOVA.stan", 
# 	#file="M5.stan", 
# 	model_name = "Mfinal.redo", 
# 	data=data_list, iter=3000, warmup=2500, chains=2, 
# 	cores=2, 
# 	control = list(stepsize = 0.5, adapt_delta = 0.99, max_treedepth = 20))
# 
# # Optional code to save the output
# save(mfinal.stan, file = "stan.model.output.rda")

# Load stan output for diagnostic plots
load("/Users/denajaneclink/Desktop/stan.model.output.rda")

## Check model output
# Create traceplots to check for mixing for site level variance
draws <- as.array(mfinal.stan, pars="ICC_site")
mcmc_trace(draws)
stan_dens(mfinal.stan, pars=c("ICC_site"))
round(summary(mfinal.stan, pars=c("ICC_site"))$summary, 3)

# Create traceplots to check for mixing for group level variance
draws <- as.array(mfinal.stan, pars="ICC_group")
mcmc_trace(draws)
stan_dens(mfinal.stan, pars=c("ICC_group"))
round(summary(mfinal.stan, pars=c("ICC_group"))$summary, 3)

# Check degrees of freedom parameter
stan_dens(mfinal.stan, pars=c("DF_obs"))
round(summary(mfinal.stan, pars=c("DF_obs"))$summary, 3)


stan_dens(mfinal.stan, pars=c("DF_group"))
round(summary(mfinal.stan, pars=c("DF_group"))$summary, 3)


stan_dens(mfinal.stan, pars=c("DF_site"))
round(summary(mfinal.stan, pars=c("DF_site"))$summary, 3)


