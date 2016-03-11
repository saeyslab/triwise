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
angledens = density.circular(circular(samples$angle), circular(anglesoi), bw=20)$y
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
angle = circular(2.5)
angledens = density.circular(circular(samples$angle), circular(anglesoi), bw=20)$y
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
backmodel = triwise::backgroundModel(angles, 10000, 10, 10, 5)
anglesoi = runif(10, 0, 1)
#anglesoi = runif(10, 0, 2*pi)
angle = circularMean(anglesoi)
angle2 = round((angle %% (2*pi)) * 2 * pi / 10)+1
z = circularZ(anglesoi)
higher = backmodel$z > z
sum(backmodel$weights[angle2,higher])/sum(backmodel$weights[angle2,])

samples = data.frame(angle=backmodel$angles, z=backmodel$z)

##
# generate several background models for different n given a set of angles
#
barycoords = transformBarycentric(Eoi)
angles = barycoords[Gdiffexp, "angle"]
noi = seq(5, 50, 3)
nanglesoi = 36
backmodels = lapply(noi, function(n) {
  backmodel = triwise::backgroundModel(angles, 100000, n, nanglesoi, 20)
})

empiricalPvalue = function(angles, backmodels) {
  nid = min(max(round(length(angles)/
                        3), 1), length(noi))

  backmodel = backmodels[[nid]]

  z = circularZ(angles)

  angle = circularMean(angles) %% (2*pi)
  angleid = ((round(angle / (2 * pi) * nanglesoi)) %% (nanglesoi-1)) + 1

  higher = backmodel$z > z
  sum(backmodel$weights[angleid,higher])# * backmodel$anglep[[angleid]] # weights were already normalized
}

sum(sapply(c(1:1000), function(i) {empiricalPvalue(runif(5, -1, 1), backmodels) < 0.01}))
sum(sapply(c(1:1000), function(i) {empiricalPvalue(runif(5, 1, 3), backmodels) < 0.01}))
sum(sapply(c(1:1000), function(i) {empiricalPvalue(runif(5, 2, 4), backmodels) < 0.01}))
sum(sapply(c(1:1000), function(i) {empiricalPvalue(runif(5, 3, 5), backmodels) < 0.01}))
sum(sapply(c(1:1000), function(i) {empiricalPvalue(runif(5, 4, 6), backmodels) < 0.01}))

names(angles) = Gdiffexp
background = names(angles)
testUnidirectionality = function(angles, gsets, Gdiffexp=NULL, backmodels=NULL, minknown=2, minfound=2, maxknown=500, weight=T, angleweights = NULL) {
  background = Gdiffexp
  scores = bind_rows(mclapply(names(gsets), function(gsetid){
    gset = gsets[[gsetid]]
    if (length(gset) < minknown | length(gset) > maxknown) {
      return(NULL)
    }

    gset_filtered = intersect(gset, background)

    if (length(gset_filtered) < minfound){
      return(NULL)
    }

    angles_gset = circular::circular(angles[gset_filtered], type="angles", units="radians")
    angle = as.numeric(triwise::circularMean(angles_gset)) %% (2*pi)

    data.frame(pval2=triwise::testRayleigh(angles_gset), pval=empiricalPvalue(angles_gset, backmodels), angle=angle, n=length(gset_filtered), gsetid=gsetid)
  }, mc.cores = 1))
  if (nrow(scores) > 0) {
    scores$qval = p.adjust(scores$pval, method="fdr")
  }
  scores
}
scores = testUnidirectionality(angles, gsets, Gdiffexp, backmodels, minknown = 4, minfound=4)

ggplot(scores) + geom_point(aes(log10(pval), log10(pval2), color=angle)) + coord_equal() + geom_abline() + geom_hline(yintercept=log10(0.1)) + geom_vline(xintercept=log10(0.1)) + scale_color_gradientn(colours=rainbow(8))

rownames(scores) = scores$gsetid
scores = scores[order(scores$qval),]
scores$name = dplyr::left_join(scores, gsetindex, "gsetid")$name
scores$redundancy = 0.5
interactivePvalplot(scores, setNames(gsetindex$name, gsetindex$gsetid), Coi)

##
# check out an individual gene set
gset = gsets[["GO:0008527"]]
pheatmap(Eoi[gset,], scale="row")
interactiveDotplot(Eoi, Gdiffexp, Glabels, gset)
##

##
# visualize how simple kernel density estimation biases the weights towards high angle density regions because there are more angles there (which then influences the determination of the z cutoff)
# if the bandwidth is small and the number of samples high, this has only a minimal impact
angle = 1
angleid = ((round(angle / (2 * pi) * nanglesoi)) %% (nanglesoi-1)) + 1
plotdata = data.frame(angle=(backmodels[[1]]$angles)%%(2*pi), weight=backmodels[[1]]$weights[angleid,])
plotdata$interval = cut(plotdata$angle, 128, labels=F)
plotdata = ddply(plotdata, ~interval, summarise, sumw=sum(weight), meanw=mean(weight), angle=min(angle) + (max(angle)-min(angle))/2)
ggplot(plotdata) +
  geom_line(aes(angle, sumw)) +
  geom_line(aes(angle, meanw * 200), color="red") +
  geom_vline(xintercept=backmodels[[1]]$anglesoi[[angleid]] %% (2*pi))
#

##
# check whether the distribution of p-values is now uniform
#angles2 = setNames(as.numeric(circular::rvonmises(n=length(angles), mu=circular::circular(pi), kappa=1)), names(angles))
angles2 = angles
gsets2 = lapply(gsets, function(gset) {sample(Gdiffexp, length(intersect(Gdiffexp, gset)))})
#names(angles2) = sample(names(angles2), length(angles2))
backmodels2 = lapply(noi, function(n) {
  backmodel = triwise::backgroundModel(angles, 100000, n, nanglesoi, 20)
})
scores = testUnidirectionality(angles2, gsets2, Gdiffexp,backmodels2, minknown = 5, minfound=5)
hist(scores$pval)
##


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
