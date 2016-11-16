##################################################
# This code takes MLE for SNPs (betas and se.betas)
# and performs a hiearchical model for bias reduction
# the code also calcualtes the proportion of FRR explained by a set of SNPS
# incorporating:
# the uncertainty of the SNP estimates
# bias reduction
# uncertainty in the pre-specified FRR
# David Conti
#################################################

library(R2jags)
set.seed(123)
setwd("/Users/davidconti/Google Drive/Collaborations/CORECT/Manuscripts/CORECT OncoArray Main Effects/BiasCorrection/CORECT.BiasCorrection")

d <- read.table("GRS_dataset.txt", header=T, sep="\t")
Y <- d$status
W <- d[,c("SEX",
          "PC1","PC2","PC3","PC4",
          "EstherVerdi_Fire","PuertoRico","Galeon_Spain_Tribe","SWEDEN_Wolk","MECC_Sephardi","ColoCare","MECC_Jew_unknown",
          "MECC_Ashkenazi_MSKCC","MECC_nonJew_nonArab","SWEDEN_Lindblom","MCCS","ATBC",
          "MEC","USC_HRT_CRC","UK_SEARCH","NHS2")]
X <- 	d[, c("all.GRS.1","all.GRS.2","all.GRS.3","all.GRS.5","all.GRS.6","all.GRS.7")]  # GRS #4 is baseline

reg <- glm(Y~as.matrix(X)+as.matrix(W), family=binomial)
coef <- summary(reg)$coef[2:7,]

##### from data
beta.hat <- coef[,1]
M <- length(beta.hat)
se.beta <- coef[,2]
prec.beta <- 1/(se.beta^2)

M <- length(beta.hat)
zeros <- rep(0, M)

###### priors on effects
sigma <- c(0.35, 0.25, 0.15, 0.15, 0.25, 0.35)
sigma2 <- sigma^2
for(i in 1:M) { print(c(exp(0-1.96*sigma[i]),exp(0+1.96*sigma[i]))) }


###### Run JAGs hierarchical model
# define JAGS model
model.string <-
"model {
  C <- 10000 #this just has to be large enough to ensure all phi[m]'s > 0
  for(m in 1:M) {
    beta[m] ~ dnorm(beta.hat[m], prec.beta[m])  # generate the MLE effect estimates
    OR[m] <- exp(beta[m])

    # normal prior on beta using the zeroes trick
    phi[m] <- C-l[m]
    l[m] <- -0.5*log(2*3.14) - 0.5*log(sigma2[m]) - 0.5 * pow((beta[m]-0),2)/sigma2[m]
    zeros[m] ~ dpois(phi[m])

  }
}"

jdata <- list(beta.hat=beta.hat,  prec.beta=prec.beta,
              M=M, sigma2=sigma2,
              zeros=zeros)

var.s <- c("OR", "beta")
model.fit <- jags.model(file=textConnection(model.string), data=jdata, n.chains=1, n.adapt=5000)
update(model.fit, n.iter=10000, progress.bar="text")
model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=20000, thin=2, quiet=F)

model.coda <- as.data.frame(model.fit[[1]])
est <- apply(model.coda, 2, mean)

beta.est <- est[grep("beta", names(est))]
OR.est <- est[grep("OR", names(est))]


# output results
# plot of SNP bias reduction
pdf("SumStat.BiasReduced.RiskScores.Estimates.pdf")
plot(beta.hat, beta.est, pch=16, xlab="MLE estimate", ylab="Biased Reduced Estimate", ylim=c(-1.5,1.5), xlim=c(-1.5,1.5))
abline(a=0, b=1)
dev.off()

r <- summary(model.fit)
write.table(r$statistics, file="Summary.SumStat.BiasReduced.RiskScores.Estimates.txt", sep="\t")
write.table(r$quantiles, file="Summary.SumStat.BiasReduced.RiskScores.Estimates.txt", sep="\t", append=T)

