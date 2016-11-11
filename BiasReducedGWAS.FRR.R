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
setwd("/Users/davidconti/Documents/Collaborations/CORECT/Manuscripts/CORECT OncoArray Main Effects/BiasCorrection")

d <- read.table("known_novel_hits_risk_allele.txt", header=T, sep="\t")

beta.hat <- d$logOR
prec.beta <- 1/(d$SE^2)
p <- d$FRQ.EUR.p
q <- 1-p
M <- length(beta.hat)
zeros <- rep(0, M)
novel <- d$NOVEL+1  # novel indexed with 2, known with 1
include.FRR <- ifelse(d$PRUNE=="NO", 1, 0)

###### priors on effects
sigma2 <- c(0.25^2, 0.05^2)
# prior on novel SNPs
# 95% prior certainty interval =
c(exp(0-1.96*sqrt(sigma2[2])),exp(0+1.96*sqrt(sigma2[2])))
# prior on known SNPs
# 95% prior certainty interval = 
c(exp(0-1.96*sqrt(sigma2[1])),exp(0+1.96*sqrt(sigma2[1])))

# uncertainty on known FRR
lambda0.m <- 1.9
se0 <- 1/7
c(lambda0.m-(1.96*se0),(lambda0.m+ 1.96*se0))
lambda0.prec <- 1/se0^2


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
    l[m] <- -0.5*log(2*3.14) - 0.5*log(sigma2[novel[m]]) - 0.5 * pow((beta[m]-0),2)/sigma2[novel[m]]
    zeros[m] ~ dpois(phi[m])

    # calculate components for FRR
    lambda[m] <- (p[m]*OR[m]*OR[m] + q[m])/((p[m]*OR[m] + q[m])*(p[m]*OR[m] + q[m]))
    log.lambda[m] <- log(lambda[m])
  }

  # calculate proportion of FRR
  FRR <- inprod(log.lambda[], include.FRR[])/log(lambda0)

  # prior on lambda0 (known familial relative risk
  lambda0 ~ dnorm(lambda0.m, lambda0.prec)
}"

jdata <- list(beta.hat=beta.hat,  prec.beta=prec.beta,
              p=p, q=q, M=M, sigma2=sigma2,
              lambda0.m=lambda0.m, lambda0.prec=lambda0.prec,
              novel=novel, include.FRR=include.FRR, zeros=zeros)

var.s <- c("OR", "FRR", "beta")
model.fit <- jags.model(file=textConnection(model.string), data=jdata, n.chains=1, n.adapt=5000)
update(model.fit, n.iter=10000, progress.bar="text")
model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=20000, thin=2, quiet=F)

model.coda <- as.data.frame(model.fit[[1]])
est <- apply(model.coda, 2, mean)
se <- apply(model.coda,2,sd)
lower <- apply(model.coda,2,function(x) {quantile(x,0.025)})
upper <- apply(model.coda,2,function(x) {quantile(x,0.975)})

lambda.est <- est[grep("lambda", names(est))]
beta.est <- est[grep("beta", names(est))]
OR.est <- est[grep("OR", names(est))]
FRR.est <- est[grep("FRR", names(est))]


# output results
# plot of SNP bias reduction
pdf("BiasReducedGWAS.Estimates.pdf")
plot(exp(beta.hat), exp(beta.est), pch=16, col=novel, xlab="MLE estimate", ylab="Biased Reduced Estimate")
abline(a=0, b=1)
legend(1, 1.15, legend=c("Known", "Novel"), col=c(1,2), pch=16)
dev.off()

r <- summary(model.fit)
write.table(round(r$quantiles, 3), file="SummaryBiasReducedGWAS.FRR.Estimates.txt", sep="\t")

