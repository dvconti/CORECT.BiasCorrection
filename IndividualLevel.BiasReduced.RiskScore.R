##################################################
# This code performs a hiearchical model for bias reduction
# for the calculation of risk score comparisons using individual level data
# David Conti
#################################################

library(R2jags)
set.seed(123)
setwd("/Users/davidconti/Documents/Collaborations/CORECT/Manuscripts/CORECT OncoArray Main Effects/BiasCorrection")

# get and define data variables
d <- read.table("GRS_dataset.txt", header=T, sep="\t")
N <- nrow(d)
Y <- d$status
W <- d[,c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
				"UK_SEARCH","EstherVerdi_Fire","MECC_Ashkenazi_MSKCC","MECC_Sephardi","MECC_nonJew_nonArab","MECC_Jew_unknown","Galeon_Spain_Tribe","PuertoRico","CGN_CFR_Kiel_Mavericc_Moffitt","SWEDEN_Lindblom","ColoCare","NHS2","MEC","ATBC","SWEDEN_Wolk","MCCS","USC_HRT_CRC")]
P <- ncol(W)
X <- 	d[, c("all.GRS.1","all.GRS.2","all.GRS.3","all.GRS.5","all.GRS.6","all.GRS.7")]  # GRS #4 is baseline
Q <- ncol(X)

# Define prior variance for each risk group - as comparred to the IQR group
sigma <- c(0.35, 0.25, 0.15, 0.15, 0.25, 0.35)
tau <- 1/(sigma^2)

for(i in 1:Q) { print(c(exp(0-1.96*sigma[i]),exp(0+1.96*sigma[i]))) }

######## Define the JAGS model
model.string <-
"model {
    for(i in 1:N) {
      Y[i] ~ dbern(pi[i])
      logit(pi[i]) <- alpha + inprod(X[i,], beta[1:Q]) + inprod(W[i,], gamma[1:P])
    }
    for(p in 1:P) { gamma[p] ~ dnorm(0, 0.00001) }
    for(q in 1:Q) {
      beta[q] ~ dnorm(0, tau[q])
      OR[q] <- exp(beta[q])
    }
    alpha ~ dnorm(0, 0.00001)
}"

jdata <- list(Y=Y, X=X, W=W, Q=Q, N=1000, P=P, tau=tau)
var.s <- c("beta", "OR")
model.fit <- jags.model(file=textConnection(model.string), data=jdata, n.chains=1, n.adapt=1000)
update(model.fit, n.iter=1000, progress.bar="text")
model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=4000, thin=2, quiet=F)

model.coda <- as.data.frame(model.fit[[1]])
est <- apply(model.coda, 2, mean)
se <- apply(model.coda,2,sd)
lower <- apply(model.coda,2,function(x) {quantile(x,0.025)})
upper <- apply(model.coda,2,function(x) {quantile(x,0.975)})


beta.est <- est[grep("beta", names(est))]
OR.est <- est[grep("OR", names(est))]

r <- summary(model.fit)
write.table(round(r$quantiles, 3), file="Summary.IndLevel.BiasReduced.RiskScore.txt", sep="\t")

