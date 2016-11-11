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

# investigating ranges of shrinkage
beta.hat <- log(seq(from=1.02, to=1.5, by=0.01))
M <- length(beta.hat)
se.beta <- rep(0.015, M)
prec.beta <- 1/(se.beta^2)

zeros <- rep(0, M)

###### priors on effects
sigma <- rep(0.05, M)
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
se <- apply(model.coda,2,sd)
lower <- apply(model.coda,2,function(x) {quantile(x,0.025)})
upper <- apply(model.coda,2,function(x) {quantile(x,0.975)})

beta.est <- est[grep("beta", names(est))]
OR.est <- est[grep("OR", names(est))]


# output results
# plot of SNP bias reduction
pdf("SumStat.BiasReduced.RiskScores.Performance.pdf")
plot(exp(beta.hat), exp(beta.est), typ="l", lwd=2, xlab="MLE estimate", ylab="Biased Reduced Estimate", ylim=c(1,max(exp(beta.hat))), xlim=c(1,max(exp(beta.hat))))
abline(a=0, b=1, lty=2)
dev.off()

r <- summary(model.fit)
write.table(round(r$quantiles, 3), file="Summary.SumStat.BiasReduced.RiskScores.Performance.txt", sep="\t")

