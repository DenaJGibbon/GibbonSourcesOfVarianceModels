### R Code for Clink et al. A multivariate analysis of female Bornean gibbon calls reveals substantial inter-individual variation and site-level variation in trill rate 
### Code to re-create figures
### Load required libraries 
library(rstan)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(gplots)
library(corrplot)
library(viridis)

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

# Load stan output
# NOTE: Must run MANOVA model code first
load("/Users/denajaneclink/Desktop/stan.model.output.rda")

## Extract the posterior samples
## These are permuted samples; we are talking permuted because output structure is more simple
## Once extracted can't use to check mixing because no longer in order along markov chain
post.samples <- extract(mfinal.stan,
                        pars= c("ICC_site", "ICC_group", "Maha_sqd", "DF_site","DF_group",
                                "DF_obs", "Vcov_site", "Vcov_group", "Vcov_obs"),
                        permuted=TRUE)

# Create vector with feature names
features <- c("note.1.dur", "note.1.f95","note.2.dur", "note.2.f95",
              "rest.dur", "intro.dur", "trill.dur", "trill.rate")

# Add column names to samples
colnames(post.samples$ICC_group) <- features
colnames(post.samples$ICC_site) <- features

# Convert to dataframe
icc.site <- as.data.frame(post.samples$ICC_site)
icc.group <- as.data.frame(post.samples$ICC_group)

# Calculate ICC for observation-level and add column names
icc.obs <- as.data.frame(1 - icc.group - icc.site)
colnames(icc.obs) <- features


### Code to make Figure 3. Posterior densities for the intraclass correlation coefficients for the three levels in our dataset 
### (call, female and site) for each feature of the Bornean gibbon great call. 
## Note 1 duration
note.1.dur.site <- cbind.data.frame(icc.site$note.1.dur, rep("Site",length(icc.site$note.1.dur)))
note.1.dur.group <- cbind.data.frame(icc.group$note.1.dur,rep("Group",length(icc.group$note.1.dur)))
note.1.dur.obs <- cbind.data.frame(icc.obs$note.1.dur,rep("Obs",length(icc.obs$note.1.dur)))

colnames(note.1.dur.site) <- c("note.1.dur.samples","icc")
colnames(note.1.dur.group) <- c("note.1.dur.samples","icc")
colnames(note.1.dur.obs) <- c("note.1.dur.samples","icc")

note.1.dur.dens.df <- rbind.data.frame(note.1.dur.site,note.1.dur.group,note.1.dur.obs)
head(note.1.dur.dens.df)

note.1.dur.plot <- ggplot(note.1.dur.dens.df, aes(x=note.1.dur.samples, fill=icc)) + geom_density()+xlim(0,1)+#+ylim(0,4)+
  scale_fill_manual(name = "icc",
                    values = c(alpha("red", .3),
                               alpha("forestgreen", .3),
                               alpha("blue", .3)))+  theme_classic() +
  guides(fill=F)+xlab("")+ylab("") +theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                                          axis.text.x  = element_text(size=24))
#note.1.dur.plot

## Note 1 max freq

note.1.f95.site <- cbind.data.frame(icc.site$note.1.f95, rep("Site",length(icc.site$note.1.f95)))
note.1.f95.group <- cbind.data.frame(icc.group$note.1.f95,rep("Female",length(icc.group$note.1.f95)))
note.1.f95.obs <- cbind.data.frame(icc.obs$note.1.f95,rep("Call",length(icc.group$note.1.f95)))

colnames(note.1.f95.site) <- c("note.1.f95.samples","icc")
colnames(note.1.f95.group) <- c("note.1.f95.samples","icc")
colnames(note.1.f95.obs) <- c("note.1.f95.samples","icc")

note.1.f95.dens.df <- rbind.data.frame(note.1.f95.site,note.1.f95.group,note.1.f95.obs)
head(note.1.f95.dens.df)

note.1.f95.plot <- ggplot(note.1.f95.dens.df, aes(x=note.1.f95.samples, fill=icc)) + geom_density()+xlim(0,1)+ #ylim(0,90)+
  scale_fill_manual(name = "ICC",
                    values = c(alpha("red", .3),
                               alpha("forestgreen", .3),
                               alpha("blue", .3)))+  theme_classic() +
  guides(fill=F)+xlab("")+ylab("") +theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                                          axis.text.x  = element_text(size=24),legend.text=element_text(size=24))+
  xlab("")+ylab("") + theme(legend.position = c(.85,.7))+
  guides(fill = guide_legend(keywidth = 4, title="",keyheight = 2))
#note.1.f95.plot

## Note 2 duration
note.2.dur.site <- cbind.data.frame(icc.site$note.2.dur, rep("Site",length(icc.site$note.2.dur)))
note.2.dur.group <- cbind.data.frame(icc.group$note.2.dur,rep("Group",length(icc.group$note.2.dur)))
note.2.dur.obs <- cbind.data.frame(icc.obs$note.2.dur,rep("Obs",length(icc.group$note.2.dur)))

colnames(note.2.dur.site) <- c("note.2.dur.samples","icc")
colnames(note.2.dur.group) <- c("note.2.dur.samples","icc")
colnames(note.2.dur.obs) <- c("note.2.dur.samples","icc")

note.2.dur.dens.df <- rbind.data.frame(note.2.dur.site,note.2.dur.group,note.2.dur.obs)
head(note.2.dur.dens.df)

note.2.dur.plot <- ggplot(note.2.dur.dens.df, aes(x=note.2.dur.samples, fill=icc)) + geom_density()+xlim(0,1)+ #ylim(0,90)+
  scale_fill_manual(name = "ICC",
                    values = c(alpha("red", .3),
                               alpha("forestgreen", .3),
                               alpha("blue", .3)))+  theme_classic() +
  guides(fill=F)+xlab("")+ylab("") +theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                                          axis.text.x  = element_text(size=24))
#note.2.dur.plot


## Note 2 max freq
note.2.f95.site <- cbind.data.frame(icc.site$note.2.f95, rep("Site",length(icc.site$note.2.f95)))
note.2.f95.group <- cbind.data.frame(icc.group$note.2.f95,rep("Group",length(icc.group$note.2.f95)))
note.2.f95.obs <- cbind.data.frame(icc.obs$note.2.f95,rep("Obs",length(icc.group$note.2.f95)))

colnames(note.2.f95.site) <- c("note.2.f95.samples","icc")
colnames(note.2.f95.group) <- c("note.2.f95.samples","icc")
colnames(note.2.f95.obs) <- c("note.2.f95.samples","icc")

note.2.f95.dens.df <- rbind.data.frame(note.2.f95.site,note.2.f95.group,note.2.f95.obs)
head(note.2.f95.dens.df)

note.2.f95.plot <- ggplot(note.2.f95.dens.df, aes(x=note.2.f95.samples, fill=icc)) + geom_density()+xlim(0,1)+ #ylim(0,90)+
  scale_fill_manual(name = "ICC",
                    values = c(alpha("red", .3),
                               alpha("forestgreen", .3),
                               alpha("blue", .3)))+  theme_classic() +
  guides(fill=F)+xlab("")+ylab("") +theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                                          axis.text.x  = element_text(size=24))
#note.2.f95.plot

## Rest duration
rest.dur.site <- cbind.data.frame(icc.site$rest.dur, rep("Site",length(icc.site$rest.dur)))
rest.dur.group <- cbind.data.frame(icc.group$rest.dur,rep("Group",length(icc.group$rest.dur)))
rest.dur.obs <- cbind.data.frame(icc.obs$rest.dur,rep("Obs",length(icc.group$rest.dur)))

colnames(rest.dur.site) <- c("trill.samples","icc")
colnames(rest.dur.group) <- c("trill.samples","icc")
colnames(rest.dur.obs) <- c("trill.samples","icc")

rest.dur.df <- rbind.data.frame(rest.dur.site,rest.dur.group,rest.dur.obs)
head(rest.dur.df)

rest.dur.plot <- ggplot(rest.dur.df, aes(x=trill.samples, fill=icc)) +xlim(0,1)+
  geom_density()+
  scale_fill_manual(name = "ICC",
                    values = c(alpha("red", .3),
                               alpha("forestgreen", .3),
                               alpha("blue", .3)))+  theme_classic() +
  guides(fill=F)+xlab("")+ylab("") +theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                                          axis.text.x  = element_text(size=24))
#rest.dur.plot

## Intro duration
intro.dur.site <- cbind.data.frame(icc.site$intro.dur, rep("Site",length(icc.site$intro.dur)))
intro.dur.group <- cbind.data.frame(icc.group$intro.dur,rep("Group",length(icc.group$intro.dur)))
intro.dur.obs <- cbind.data.frame(icc.obs$intro.dur,rep("Obs",length(icc.group$intro.dur)))

colnames(intro.dur.site) <- c("trill.samples","icc")
colnames(intro.dur.group) <- c("trill.samples","icc")
colnames(intro.dur.obs) <- c("trill.samples","icc")

intro.dur.df <- rbind.data.frame(intro.dur.site,intro.dur.group,intro.dur.obs)
head(intro.dur.df)

intro.dur.plot <- ggplot(intro.dur.df, aes(x=trill.samples, fill=icc)) + xlim(0,1)+
  geom_density()+
  scale_fill_manual(name = "ICC",
                    values = c(alpha("red", .3),
                               alpha("forestgreen", .3),
                               alpha("blue", .3)))+  theme_classic() +
  guides(fill=F)+xlab("")+ylab("") +theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                                          axis.text.x  = element_text(size=24))
#intro.dur.plot


## Trill duration
trill.dur.site <- cbind.data.frame(icc.site$trill.dur, rep("Site",length(icc.site$trill.dur)))
trill.dur.group <- cbind.data.frame(icc.group$trill.dur,rep("Group",length(icc.group$trill.dur)))
trill.dur.obs <- cbind.data.frame(icc.obs$trill.dur,rep("Obs",length(icc.group$trill.dur)))

colnames(trill.dur.site) <- c("trill.samples","icc")
colnames(trill.dur.group) <- c("trill.samples","icc")
colnames(trill.dur.obs) <- c("trill.samples","icc")

trill.dens.df <- rbind.data.frame(trill.dur.site,trill.dur.group,trill.dur.obs)
head(trill.dens.df)

trill.dur.plot <- ggplot(trill.dens.df, aes(x=trill.samples, fill=icc)) + geom_density()+
  scale_fill_manual(name = "ICC",
                    values = c(alpha("red", .3),
                               alpha("forestgreen", .3),
                               alpha("blue", .3)))+  theme_classic() +
  guides(fill=F)+xlab("")+ylab("") +theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                                          axis.text.x  = element_text(size=24))

#trill.dur.plot


## Trill rate
trill.rate.site <- cbind.data.frame(icc.site$trill.rate, rep("Site",length(icc.site$trill.rate)))
trill.rate.group <- cbind.data.frame(icc.group$trill.rate,rep("Group",length(icc.group$trill.rate)))
trill.rate.obs <- cbind.data.frame(icc.obs$trill.rate,rep("Obs",length(icc.group$trill.rate)))

colnames(trill.rate.site) <- c("trill.samples","icc")
colnames(trill.rate.group) <- c("trill.samples","icc")
colnames(trill.rate.obs) <- c("trill.samples","icc")

trill.dens.df <- rbind.data.frame(trill.rate.site,trill.rate.group,trill.rate.obs)
head(trill.dens.df)

trill.rate.plot <- ggplot(trill.dens.df, aes(x=trill.samples, fill=icc)) + geom_density()+
  scale_fill_manual(name = "ICC",
                    values = c(alpha("red", .3),
                               alpha("forestgreen", .3),
                               alpha("blue", .3)))+  theme_classic() +
  guides(fill=F)+xlab("")+ylab("") +theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                                          axis.text.x  = element_text(size=24))
#trill.rate.plot


## Combine all plots
plot_grid(hjust=-0.75, label_size = 20, vjust=1,note.1.dur.plot,note.1.f95.plot,note.2.dur.plot,note.2.f95.plot,rest.dur.plot,intro.dur.plot,trill.dur.plot,trill.rate.plot,labels=c("Note 1 Duration","Note 1 Max Freq","Note 2 Duration","Note 2 Max Freq", "Rest Duration","Intro Duration","Trill Duration", "     Trill Rate"),ncol=2)



### Code to make Figure 4. Posterior densities for site-specific intercepts for trill rate in the Bornean gibbon female great call. 

## Extract site-specific random intercepts for trill rate
site.intercepts <- extract(mfinal.stan, pars="site_rand_intercept")$site_rand_intercept
trill.intercepts <- data.frame(site.intercepts[, , 8]) # Intercepts for the 8th feature. 
names(trill.intercepts) <- c("CR", "DK", "DV", "IC", "KB", "MB", "SAF")
str(trill.intercepts)

## Convert data into dataframe to pass to ggplot
CR <- cbind.data.frame(trill.intercepts$CR, rep("Crocker",length(trill.intercepts$CR)))
colnames(CR) <- c("samples","site")
DK <- cbind.data.frame(trill.intercepts$DK, rep("Dermakot",length(trill.intercepts$DK)))
colnames(DK) <- c("samples","site")
DV <- cbind.data.frame(trill.intercepts$DV, rep("Danum",length(trill.intercepts$DV)))
colnames(DV) <- c("samples","site")
IC <- cbind.data.frame(trill.intercepts$IC, rep("Imbak",length(trill.intercepts$IC)))
colnames(IC) <- c("samples","site")
KB <- cbind.data.frame(trill.intercepts$KB, rep("Kinabatangan",length(trill.intercepts$KB)))
colnames(KB) <- c("samples","site")
MB <- cbind.data.frame(trill.intercepts$MB, rep("Maliau",length(trill.intercepts$MB)))
colnames(MB) <- c("samples","site")
SAF <- cbind.data.frame(trill.intercepts$SAF, rep("Kalabakan",length(trill.intercepts$CR)))
colnames(SAF) <- c("samples","site")
random.intercept.df <- rbind.data.frame(CR, DK, DV, IC, KB, MB, SAF)

# Check structure of data frame
head(random.intercept.df)
dim(random.intercept.df)

## Create color palette for figure 
my.cols <- viridis(n=7)


## Create Figure 4 in ggplot
ggplot(random.intercept.df, aes(x=samples, fill=site))+ geom_density(alpha=.45, bw=.015)+
  xlab("")+ylab("")+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                          axis.text.x  = element_text(size=20))+
  scale_fill_manual(values = my.cols)+
  theme_classic()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                        axis.text.x  = element_text(size=20))+
  guides(fill = guide_legend(title="Site"))+
  theme(legend.text=element_text(size=20))+
  theme(legend.title =element_text(size=20))


### Code to create Figure 5. Correlation matrices for the eight features at each level of analysis: call, female and site. 

# Extract the variance/co-variance matrix for each level of analysis
Vcov.site <- extract(mfinal.stan,
                     pars= c("Vcov_site"),
                     permuted=TRUE)$Vcov_site


Vcov.female <- extract(mfinal.stan,
                       pars= c("Vcov_group"),
                       permuted=TRUE)$Vcov_group

Vcov.obs <- extract(mfinal.stan,
                    pars= c("Vcov_obs"),
                    permuted=TRUE)$Vcov_obs


### Turn each variance/co-variance matrix into correlation matrix

Vcov.site <- apply(X=Vcov.site, MAR=1, FUN="cov2cor")
postmean.vcov.cor.site <- apply(X=Vcov.site, MAR=1, FUN ="mean")
postmean.vcov.cor.site <- matrix(postmean.vcov.cor.site, nrow= 8, ncol=8, byrow=T)
row.names(postmean.vcov.cor.site) <- features
colnames(postmean.vcov.cor.site) <- features


Vcov.female <- apply(X=Vcov.female, MAR=1, FUN="cov2cor")
postmean.vcov.cor.female <- apply(X=Vcov.female, MAR=1, FUN ="mean")
postmean.vcov.cor.female <- matrix(postmean.vcov.cor.female, nrow= 8, ncol=8, byrow=T)
row.names(postmean.vcov.cor.female) <- features
colnames(postmean.vcov.cor.female) <- features


Vcov.cor.obs <- apply(X=Vcov.obs, MAR=1, FUN="cov2cor")
postmean.vcov.cor.obs <- apply(X=Vcov.cor.obs, MAR=1, FUN ="mean")
postmean.vcov.cor.obs <- matrix(postmean.vcov.cor.obs, nrow= 8, ncol=8, byrow=T)
row.names(postmean.vcov.cor.obs) <- features
colnames(postmean.vcov.cor.obs) <- features


## Call-level matrix
corrplot(postmean.vcov.cor.obs, method=c("circle"), 
         col=viridis(10), cl.pos = "n",
         type="lower",
         #title= "Call", 
         tl.cex = 1.2,                                        
         tl.srt = 45, diag = F)
title(main="Call", cex.main=2)


## Female-level matrix
corrplot(postmean.vcov.cor.female,   method=c("circle"), 
         col=viridis(10), cl.pos = "n",
         type="lower",
         #title= "Site",
         tl.cex =  1.2, tl.srt = 45, diag = F)
title(main="Female", cex.main=2)


## Site-level matrix
corrplot(postmean.vcov.cor.site,method=c("circle"), 
         col=viridis(10), cl.pos = "n",
         type="lower",
         #title= "Site",
         tl.cex =  1.2, tl.srt = 45, diag = F)
title(main="Site", cex.main=2)


### Code to create Figure 6. Posterior mean Mahalanobis distances, squared and scaled by the number of features, 
### versus F distribution quantiles.

# Extract the average Mahalanobis squared distances over realizations.
# Divide by K = dim(d.manova)[2] to scale for comparison with the F distribution. 
Maha_sqd_scaled <- summary(mfinal.stan, pars=c("Maha_sqd"))$summary[,"mean"]/dim(d.manova)[2]
# Extract the average multivariate-t degrees of freedom over realizations. 
t.df <- summary(mfinal.stan, pars=c("DF_obs"))$summary[,"mean"]
# Plot the distances -vs- F-distribution (df1= K, df2=t.df) quantiles. 
plot(qf(ppoints(Maha_sqd_scaled), df1=dim(d.manova)[2], df2=t.df), sort(Maha_sqd_scaled), 
     xlab="F-distribution quantile", ylab="Scaled Mahalanobis Distance", 
     bty="n")
abline(a=0, b=1)


