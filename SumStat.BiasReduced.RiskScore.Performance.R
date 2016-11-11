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
beta.hat <- log(seq(from=1.01, to=1.1, by=0.01))
M <- length(beta.hat)
se.beta <- seq(from=0.01, to=0.02, by=0.001)
N <- length(se.beta)
beta.hat <- rep(beta.hat, N)
se.beta <- rep(se.beta, each=M)
prec.beta <- 1/(se.beta^2)
p.value <- 2*(1-pnorm(beta.hat/se.beta))

zeros <- rep(0, M*N)

###### priors on effects
sigma <- rep(0.05, M*N)
sigma2 <- sigma^2
for(i in 1:M*N) { print(c(exp(0-1.96*sigma[i]),exp(0+1.96*sigma[i]))) }


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
              M=M*N, sigma2=sigma2,
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
pdf("SumStat.BiasReduced.RiskScores.SE.Impact.pdf")
plot(exp(beta.hat[se.beta==se.beta[1]]), exp(beta.est[se.beta==se.beta[1]]), typ="l", lwd=2, xlab="MLE estimate", ylab="Biased Reduced Estimate", ylim=c(1,max(exp(beta.hat))), xlim=c(1,max(exp(beta.hat))))
l.se.betas <- unique(se.beta)
index <- 2
for(s in l.se.betas[l.se.betas!=se.beta[1]]) {
  lines(exp(beta.hat[se.beta==s]), exp(beta.est[se.beta==s]), col=index)
  index <- index + 1
}
legend(1, max(exp(beta.est)), legend=c(paste("SE=", l.se.betas, sep="")), lty=1, col=1:N, cex=.5)
dev.off()
pdf("SumStat.BiasReduced.RiskScores.Beta.Impact.pdf")
plot(se.beta[beta.hat==beta.hat[1]], exp(beta.est[beta.hat==beta.hat[1]]), typ="l", lwd=2, xlab="SE of the MLE estimate", ylab="Biased Reduced Estimate", ylim=c(1,max(exp(beta.hat))))
abline(a=0, b=1, lty=2)
l.betas <- unique(beta.hat)
index <- 2
for(s in l.betas[l.betas!=beta.hat[1]]) {
  lines(se.beta[beta.hat==s], exp(beta.est[beta.hat==s]), col=index)
  index <- index + 1
}
legend(se.beta[1], max(exp(beta.est)), legend=c(paste("OR.hat=", round(exp(l.betas),3), sep="")), lty=1, col=1:M, cex=.5)
dev.off()


pdf("SumStat.BiasReduced.RiskScores.P.Value.Impact.pdf")
l.betas <- unique(beta.hat)
plot(-log(p.value), exp(beta.est), pch=16, col=match(beta.hat, l.betas), xlab="-log(P.Value) of the MLE estimate", ylab="Biased Reduced Estimate", ylim=c(1,max(exp(beta.hat))))
index <- 1
for(s in l.betas) {
  abline(h=exp(s), col=index, lty=2, lwd=2)
  index <- index + 1
}
legend(se.beta[1], max(exp(beta.est)), legend=c(paste("OR.hat=", round(exp(l.betas),3), sep="")), lty=1, col=1:M, cex=.5)
dev.off()

r <- summary(model.fit)
write.table(round(r$quantiles, 3), file="Summary.SumStat.BiasReduced.RiskScores.Performance.txt", sep="\t")

