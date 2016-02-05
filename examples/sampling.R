library(plyr)
library(dplyr)
library(circular)
library(ggplot2)
library(doMC)
library(caTools)
registerDoMC(8)

angles = as.numeric(circular::rvonmises(n=2000, mu=circular::circular(pi), kappa=3) + circular::rvonmises(n=2000, mu=circular::circular(0), kappa=3)) %% (pi * 2)
angles = barycoords$angle %% (pi * 2)
hist(angles)

# why is the z-value distribution different for n?
# because with increasing n, the mean angle of a sample tends to go to the overall mean of all angles
# eg. at n=1 the distributions are the same
# with more n, the distribution goes to the overall mean

samples = ldply(runif(5, 5, n=10000), function(n) {
  n = floor(n)
  anglesample = sample(angles, n)

  data.frame(
    z=sqrt(sum(sin(anglesample))**2 + sum(cos(anglesample) )**2)/length(anglesample),
    angle=as.numeric(atan2(mean(sin(anglesample)), mean(cos(anglesample)))) %% (2*pi),
    prayleigh=as.numeric(testRayleigh(anglesample))
  )
}, .parallel=T)

# relation between sample angle and sample z
plot(samples$angle, samples$z)

samples$anglebin = as.factor(findInterval(samples$angle %% pi*2, seq(0, 2*pi, pi/12)))
ggplot(samples) + geom_violin(aes(x=anglebin, y=z))

plot(samples$anglebin, samples$z)

library(Rcpp)
cppFunction('int add(int x, int y, int z) {
  int sum = x + y + z;
  return sum;
}')

##
# these are p-values for the angle dimensions
# if your distributions of angles is uniform, every angle should have p-value 1, because no angle is more extreme than some other
# if not uniform, you need some way to tell whether a certain angle is extreme
# it is hard to define a true test statistic here, as this depends on your distributions
# i can however directly calculate a p-value by starting from the kernel density distribution of the angles
# and calculate the density under the chosen angle
# it is easy to see that in a perfectly uniform distribution this p-value will be 1
# otherwise, it tells you the chance you expect this or an even more extreme angle
# note however that combining this p-value with the z-score p-value is certainly possible
# but that testing whether this p-value works is not trivial becuase this requires samples
# and you sampling itself is already biased towards certain angles
anglesoi = seq(0, 2*pi, pi/36)
angledens = density.circular(circular(samples$angle), circular(anglesoi), bw=10)$y
pangle = trapz(anglesoi, angledens)
print(pangle)
angle = 2
pangle = trapz(anglesoi, pmin(angledens, angledens[[round(angle / (2*pi) * (length(angledens) - 1)) + 1]]))
print(pangle)
samples$pangle = sapply((samples$angle %% (2*pi)), function(angle) {trapz(anglesoi, pmin(angledens, angledens[[round(angle / (2*pi) * (length(angledens) - 1)) + 1]]))})

plot(anglesoi, angledens, ylim= c(0, max(angledens)))
##

##
# here I calculate for a particular angle its distribution of z-values under null hypothesis
# I can either calculate the exact angle density for every angle
# or, approximate it by first calculating the density at certain fixed angles and then
# getting for every sampled angle the nearest density (estimating density is very expensive)
# these angle densities are necessary to do a correct estimation of z-value density
# as otherwise regions with a lot of angles will have a disproportionate effect on the z-value distribution
# I thus calculate a p-value at different zcutoffs by comparing the number of samples with a z higher than
# the z-cutoff and then sum over their associated (normalized) weights to calculate the probability to
# get an equal or more extreme z-value
angle = circular(0)
angledens = density.circular(circular(samples$angle), circular(anglesoi), bw=3)$y
#samples$angledens = density.circular(circular(samples$angle), circular(samples$angle), bw=3)$y
samples$angledens = sapply(circular(samples$angle, modulo="2pi"), function(angle) {angledens[[round(angle / (2*pi) * (length(angledens) - 1)) + 1]]})
weights = dvonmises(circular(samples$angle, modulo="2pi"), angle, 10) / samples$angledens
cutoffs = seq(0, 1, 0.001)
pvals = sapply(cutoffs, function(zcutoff){
  higher = samples$z > zcutoff
  sum(weights[samples$z > zcutoff]) / sum(weights)
})
plot(cutoffs,pvals, type="l")
plot(as.numeric(circular(samples$angle, modulo="2pi")), weights)

##
# cutoffs for significance at different angles
anglesoi =seq(0, pi*2, pi/24)
zcutoffs = sapply(anglesoi, function(angle) {
  weights = dvonmises(circular(samples$angle, modulo="2pi"), circular(angle), 3) / samples$angledens
  cutoffs = seq(min(samples$z), max(samples$z), 0.001)
  pvals = sapply(cutoffs, function(zcutoff){
    higher = samples$z > zcutoff
    sum(weights[higher]) / sum(weights)
  })

  cutoffs[pvals < 0.1][[1]]
})
plot(anglesoi, zcutoffs)
##

##
# distribution of pvalues
samples$pempirical = mapply(samples$angle[1:2000], samples$z[1:2000], FUN=function(angle, z) {
  weights = dvonmises(circular(samples$angle, modulo="2pi"), circular(angle, modulo="2pi"), 10) / samples$angledens
  higher = samples$z > z
  sum(weights[higher]) / sum(weights)
})
hist(samples$pempirical)
##

##
# comparison between rayleigh and empirical
ggplot(samples[1:2000,]) + geom_point(aes(log10(pempirical * pangle), log10(prayleigh), color=angle)) + coord_equal() + geom_abline() + geom_hline(yintercept=log10(0.1)) + geom_vline(xintercept=log10(0.1))
##

##
# generate a background model given an n and a set of angles
#
library(Rcpp)
backmodel = backgroundModel(angles, 10000000, 10, 10, 5)
anglesoi = runif(10, 1, 2)
#anglesoi = runif(10, 0, 2*pi)
angle = circularMean(anglesoi)
angle2 = round((angle %% (2*pi)) * 2 * pi / 10)+1
z = circularZ(anglesoi)
higher = backmodel$z > z
sum(backmodel$weights[angle2,higher])/sum(backmodel$weights[angle2,])

samples = data.frame(angle=backmodel$angles, z=backmodel$z)



subz = samples$z[between(samples$angle, 1.9, 2.1)]
hist(subz, freq = F)
lambda = vglm(subz ~ 1, rayleigh)@coefficients[[1]]
plot(drayleigh(seq(0, max(samples$z)), scale = lambda))

sampledata = bind_rows(lapply(seq(0, max(samples$z), 1), function(zcutoff) {
  subsample = samples[samples$z>=zcutoff, ]
  subsample$zcutoff = zcutoff
  subsample
}))
sampledata$zcutoff = factor(sampledata$zcutoff)

ggplot() +
  geom_freqpoly(aes(angle, color=zcutoff), data=sampledata, binwidth=0.4, position="identity")

ggplot() +
  geom_histogram(aes(x=angle), data=samples[samples$z>=0,], binwidth=0.2, alpha=0.1, fill="black") +
  geom_histogram(aes(x=angle), data=samples[samples$z>1,], binwidth=0.2, alpha=0.1, fill="black") +
  geom_histogram(aes(x=angle), data=samples[samples$z>2,], binwidth=0.2, alpha=0.1, fill="black") +
  geom_histogram(aes(x=angle), data=samples[samples$z>3,], binwidth=0.2, alpha=0.1, fill="black") +
  geom_histogram(aes(x=angle), data=samples[samples$z>4,], binwidth=0.2, alpha=0.1, fill="black") +
  geom_histogram(aes(x=angle), data=samples[samples$z>5,], binwidth=0.2, alpha=0.1, fill="black") +
  geom_histogram(aes(x=angle), data=samples[samples$z>6,], binwidth=0.2, alpha=0.1, fill="black") +
  geom_histogram(aes(x=angle), data=samples[samples$z>7,], binwidth=0.2, alpha=0.1, fill="black") +
  geom_histogram(aes(x=angle), data=samples[samples$z>8,], binwidth=0.2, alpha=0.1, fill="black")


hist(samples$angle[samples$z > 1])
