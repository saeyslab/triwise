#angles = circular::circular(barycoords$angle, modulo = "2pi")
#dens = circular::density.circular(angles, bw=10000)
#plot(dens)

library(ggplot2)
library(parallel)

weight_bw = 12

scorez = function(angles) {sqrt(sum(sin(angles))**2 + sum(cos(angles))**2)/length(angles) }


gsets = lapply(runif(10, 10, n=2000), function(s) sample(Gdiffexp, floor(s)))

barycoords = transformBarycentric(Eoi)

gsetangles = circular::circular(sapply(gsets, function(gset) {circularMean(barycoords[gset, "angle"])}))
gsetzs = sapply(gsets, function(gset) scorez(barycoords[gset, "angle"]))
meanz = sapply(seq(0, pi*2, pi/12), function(mu) {
  probs = circular::dvonmises(gsetangles, circular(mu), weight_bw)
  sum(gsetzs * probs) / sum(probs)
})
angleweights = 1/(meanz/0.21)
plot(meanz, ylim=c(0,0.5))
plot(gsetzs, gsetangles)

scores = testUnidirectionality(barycoords, gsets, Gdiffexp, weight=T, angleweights = angleweights)
scores2 = testUnidirectionality(barycoords, gsets, Gdiffexp, weight=F, angleweights = angleweights)

scores$weighted = T
scores$biased = 1
scores2$weighted = F
scores2$biased = 1

barycoords2 = barycoords
barycoords2$angle = runif(0, pi*2, n = nrow(barycoords2))
gsetangles = circular::circular(sapply(gsets, function(gset) {circularMean(barycoords2[gset, "angle"])}))
gsetzs = sapply(gsets, function(gset) scorez(barycoords[gset, "angle"]))
meanz = sapply(seq(0, pi*2, pi/12), function(mu) {
  probs = circular::dvonmises(gsetangles, circular(mu), weight_bw)
  sum(gsetzs * probs) / sum(probs)
})
plot(meanz, ylim=c(0,3))
angleweights = 1/(meanz/0.21)

scores3 = testUnidirectionality(barycoords2, gsets, Gdiffexp, weight=T, angleweights = angleweights)
scores4 = testUnidirectionality(barycoords2, gsets, Gdiffexp, weight=F, angleweights = angleweights)

scores3$weighted = T
scores3$biased = 0
scores4$weighted = F
scores4$biased = 0

barycoords3 = barycoords
barycoords3$angle = as.numeric(circular::rvonmises(n=nrow(barycoords3), mu=circular::circular(pi), kappa=1))
gsetangles = circular::circular(sapply(gsets, function(gset) {circularMean(barycoords3[gset, "angle"])}))
gsetzs = sapply(gsets, function(gset) scorez(barycoords3[gset, "angle"]))
meanz = sapply(seq(0, pi*2, pi/12), function(mu) {
  probs = circular::dvonmises(gsetangles, circular(mu), 50)

  print(probs)
  sum(gsetzs * probs) / sum(probs)
})
plot(meanz, ylim=c(0,2))
angleweights = 1/(meanz/0.21)

scores5 = testUnidirectionality(barycoords3, gsets, Gdiffexp, weight=T, angleweights = angleweights)
scores6 = testUnidirectionality(barycoords3, gsets, Gdiffexp, weight=F, angleweights = angleweights)

scores5$weighted = T
scores5$biased = 2
scores6$weighted = F
scores6$biased = 2

histdata = rbind(scores, scores2, scores3, scores4, scores5, scores6)
ggplot(histdata) + geom_histogram(aes(x=pval, fill=weighted+(biased+1)**2), alpha=1, binwidth=0.05, position="identity") + facet_grid(weighted~biased) + theme_minimal()

distanceToUniformity = function(pvals) sqrt(sum((hist(pvals, plot=F, breaks = seq(0, 1, by = 0.1))$counts - length(pvals)/10)**2))

distanceToUniformity(scores$pval)
distanceToUniformity(scores2$pval)
distanceToUniformity(scores3$pval)
distanceToUniformity(scores4$pval)

##

scorebw = function(weight_bw) {
  gsetangles = circular::circular(sapply(gsets, function(gset) {circularMean(barycoords[gset, "angle"])}))
  gsetzs = sapply(gsets, function(gset) scorez(barycoords[gset, "angle"]))
  meanz = sapply(seq(0, pi*2, pi/12), function(mu) {
    probs = circular::dvonmises(gsetangles, circular(mu), weight_bw)
    sum(gsetzs * probs) / sum(probs)
  })
  angleweights = 1/(meanz/3)

  scores = testUnidirectionality(barycoords, gsets, Gdiffexp, weight=T, angleweights=angleweights)
  distanceToUniformity(scores$pval)
}

x = seq(1, 50, 4)
distances = mclapply(x, scorebw, mc.cores=8)
plot(x, distances)

optimization_result = optimize(scorebw, interval = c(1, 50/pi), tol=0.1)
weight_bw = optimization_result$minimum

##
pvals = sapply(c(1:5000), function(i) {
  angles = c(rep(1, 10), runif(500, 0, pi*2))
  weights = c(rep(5, 10), rep(1, 500))
  testRayleigh(angles, weights)
})
hist(pvals, xlim=c(0, 1), breaks=seq(0, 1, 0.05))


# check if the weighting really involves "adding" and "substracting" extra samples
testRayleigh(c(1, 1, 1), c(1,1,1))
testRayleigh(c(1, 1), c(1, 2))
testRayleigh(c(1, 1), c(0.5, 2.5))

##
angles = circular::circular(barycoords$angle, modulo = "2pi")
dens = circular::density.circular(angles, bw=weight_bw*2)
plot(dens$y)
