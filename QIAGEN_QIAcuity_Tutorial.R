## QIAcuity User Manual Extension QIAGEN: Statistics of nanoplate dPCR
## (c) 2023 Winter Maxwell Thayer, Johns Hopkins University,
##        all rights reserved.
## 2022-04-13


##---- QIAcuity User Manual Extension document ----
#  HB-2839-003_UM_QIAcuity_UM_Extension_0621_WW first section Nanoplate Digital
#  PCR contains subsection Statistics of nanoplate dPCR on page 7. This R file
#  describes appropriate procedures for replicating all analyses contained
#  therein with the R statistical environment. Note that you may want to run the
#  function rm(list=ls()) to clear your environment prior to running this
#  script. Be sure no valuable data is stored in your R environment first, maybe
#  by running ls() and looking at the output, just to be safe.
#
#  This script develops several functions in a step-by-step way. The finished
#  versions are included here for ease.
#' \code{computeLambda()} is a function that takes the number of valid
#' partitions and the numbe of positive partitions and returns lambda from the
#' Poisson distribution
computeLambda <- function(valid, positive, total=26000) {
  #! add a warning if the proportion valid is too small
  return(-log((valid - positive) / valid))
}

#' \code{computeLambdaCI()} is a function that takes a molecule count value 
#' value and a desired confidence level and computes the confidence interval for
#' a Poisson-distributed random variable
computeLambdaCI <- function (count, conf.level=0.95) {
  alpha <- 1 - conf.level
  LB <- 0.5 * qchisq(alpha / 2, 2 * count)
  UB <- 0.5 * qchisq(1 - alpha / 2, 2 * count + 2)
  return(c(LB, UB))
}


computeVolume <- function(lambda, volume=.91, inputReactionVolume=F) {
  #  if the user specified an inputReactionVolume, return it with the reaction
  #  volume
  if(inputReactionVolume) {
    return(c((lambda / volume) * 1000,
             (lambda / volume) * 1000) * inputReactionVolume)
  }
  #  if the user didn't specify an input volume, return just the reaction volume
  return((lambda / volume) * 1000)
}

##---- Poisson distribution ----
#  The first portion of the manual describes the natural distribution of the
#  dPCR quantification and suggests the use of the theoretical distribution
#  Poisson to approximate the copies of the target molecule per positive
#  partition. The manual provides an example of a theoretical Poisson
#  distribution with lambda ranging from 0.1 to 5. The formula for Poisson
#  is described and the following parameters are defined:
# 
#  e: Euler's constant, can be obtained with exp(1) in R
#  Lambda: can be thought of as the concentration of the target molecule, or,
#            the average count per microliter
#  k: the copies of the molecule per partition
#  k!: k factorial, obtained with factorial(k) in R
#
#  plots of the theoretical Poisson distribution are then presented, replicate
#  these plots with the following code.
#
#  investigate the help documentation for the poisson distribution family of
#  functions
?dpois() # dpois(x, lambda, log=FALSE)
#  change the plot layout to display five plots side-by-side in one row
layout(matrix(1:5, nrow=1))
#  store the lambda values provided in the example
myLambdas <- c(.1, .5, 1, 2, 5)
#  first try at approximating the color values used in the plots
# myColors <- c("medium blue", "light goldenrod", "tomato3", "grey40",
              # "light green")
#  try for a closer approximate using rgb values obtained with a color picker
myColors <- list(rgb(0.2, 0.349, 0.545),
                 rgb(0.969, 0.69, 0),
                 rgb(0.69, 0.251, 0.106),
                 rgb(0.365, 0.345, 0.345),
                 rgb(0.733, 0.851, 0.631))
#  plot the each lambda using a for-loop, use the same barplot function call
#  five times, each using myLambdas[i] and myColors[[i]]
#  The manual uses commas rather than periods to indicate decimal points on
#  the plot, set that option before and after the plot. Depending on where you
#  are in the world, these commands may not be needed.
options(OutDec=",")
for(i in 1:5) {
  #  create a null plot
  barplot(col = NA, border = NA, axes = FALSE, dpois(x=0:10, myLambdas[i]),
          ylim=c(0, 1), xlim=c(0, 10))
  #  add the horizontal lines first so they are below the bars
  segments(rep(0, 6), seq(0, 1, .2), rep(12, 6), seq(0, 1, .2), col="#cccccc")
  #  plot the Poisson distributions
  barplot(dpois(x=0:10, myLambdas[i]), ylim=c(0, 1), xlim=c(0, 10),
          col=eval(myColors[[i]]), border = 'transparent', add=T, las=1,
          names.arg=0:10, main=bquote(lambda==.(myLambdas[i])))
  #  plot the left hand y-axis label and the bottom x-axis label as the first
  #  and third plots come up. check par(xpd), it will change your device
  #  clipping area, i.e., where R stops drawing on your plot.
  if(i == 1) {
    par(xpd=NA) # clip to the device region
    text(-5, .3, labels="Probability P(X=k)", srt=90, adj=0, xpd=NA)
    par(xpd=F) # reset to clip to plot region
  }
  if(i == 3) {
    par(xpd=NA)
    text(0, -0.15, labels="k = copies per partition", adj=0)
    par(xpd=F)
  }
}
options(OutDec=".") ## set the decimal back to a period
##!! could lighten the blue and the grey
#  This plot shows the theoretical Poisson distribution with five different
#  lambda values. Note that as lambda increases, the the counts become
#  increasingly normally distributed.
# 
#  The last paragraph of this section describes the calculation of copies per
#  partition, which is:
#
#    copies of molecule per partition * number of valid partitions
#
#  The following section describes how to calculate the copies of the molecule
#  per partition.


##----  Absolute quantification - copies per partition ----
#  This section of the manual describes how to calculate lambda, the average
#  count of the target molecule in the well, based on empirical data. The
#  calculation is:
#
#   lambda = -ln( (valid partitions - positive partitions) /  valid partitions )
#
#  The following example is provided in the manual.
validPartitions <- 8000
poisitivePartitions <- 4000
#  note that the first and last parentheses on the next line tells R to print
#  the result in addition to storing it, rather than storing it silently
(thisLambda <- -log((validPartitions - poisitivePartitions) / validPartitions))
#  this output can be rounded to the third decimal, as in the manual
round(thisLambda, 3)

#  we can create a formula based on this so that it is easier to reuse
computeLambda <- function(valid, positive) {
  # the return is what R gives back after executing the function
  return(-log((valid - positive) / valid))
}
#  check that our output is correct
( exampleLambda <- round(computeLambda(8000, 4000), 3) )
#  we can also check that we get what we expect by using the double equals sign
#  to get a boolean value
exampleLambda == 0.693 # TRUE

#  Given this lambda, the total number of copies in all valid partitions is:
validPartitions * exampleLambda # 5544
#  note that if we used the lambda prior to rounding we obtain
validPartitions * thisLambda # 5545.177
#  Wait to round until the final step.

#  Next, the manual describes how to calculate the 95% confidence interval
#  assuming the data are Poisson-distributed. Any calculation of a confidence
#  interval can be broken into the following two parts:
# 
#   someEstimate +or- errorOfSomeEstimate
# 
#  Another way to write this is:
# 
#   lowerBound < estimate < upperBound
# 
#  The writers of this manual have used a different form of the Poisson
#  distribution to express the confidence interval than they used to express
#  lambda in the previous section:
#
#  (1/2)*chiSquare(alpha/2, 2*k) ≤ Mu ≤ (1/2)*chiSquare(1 - (alpha/2), 2*k + 2)
#
#  This expression relies on the relationship between the Poisson,
#  Chi-square, and Gamma distributions. In this application, it breaks down to,
#  we expect the average number of copies per partition to be between our
#  estimate of lambda (the average number of copies per partition) plus or minus
#  our estimate of the amount of variability in lambda based on what we know
#  about the theoretical distribution with this many degrees of freedom, in this
#  case 8000. Do not worry about degrees of freedom right now \U+1F9D8.
#  
#  The manual provides the following calculations for their example
lambdaLow <- 0.693 + ((1.96^2) / (2 * 8000)) -
                sqrt(0.693 * ((1.96^2) / 8000) + ((1.96^4) / (4 * 8000^2)))
lambdaHigh <- 0.693 + ((1.96^2) / (2 * 8000)) +
                sqrt(0.693 * ((1.96^2) / 8000) + ((1.96^4) / (4 * 8000^2)))

#  replacing this with the values we have stored, we get
lambdaLow.2 <- exampleLambda + ((1.96^2) / (2 * validPartitions)) -
  sqrt(exampleLambda * ((1.96^2) / validPartitions) +
         ((1.96^4) / (4 * validPartitions^2)))
lambdaHigh.2 <- exampleLambda + ((1.96^2) / (2 * validPartitions)) +
  sqrt(exampleLambda * ((1.96^2) / validPartitions) +
         ((1.96^4) / (4 * validPartitions^2)))
#  check that they are the same
testLambdas <- c(lambdaLow == lambdaLow.2, lambdaHigh == lambdaHigh.2)
all(testLambdas) # TRUE

#  to make that more readable, break it into smaller parts and give each a name
chiSquareFor95CI <- (1.96^2)
degreesOfFreedom <- (2 * validPartitions)
constantToAddToMean <- (chiSquareFor95CI / degreesOfFreedom)
errorEstimate <- sqrt(exampleLambda * (chiSquareFor95CI / validPartitions) +
                        ((chiSquareFor95CI^2) / (4 * validPartitions^2)))
#  calculate the confidence interval bound
lambdaLow.3 <- exampleLambda + constantToAddToMean - errorEstimate
lambdaHigh.3 <- exampleLambda + constantToAddToMean + errorEstimate
testLambdas <- c(testLambdas, lambdaLow == lambdaLow.3,
                 lambdaHigh == lambdaHigh.3)
all(testLambdas)

#  You may have noticed that these confidence intervals do not match the
#  manual. Can you tell why?
lambdaLow.4 <- thisLambda + constantToAddToMean - errorEstimate
lambdaHigh.4 <- thisLambda + constantToAddToMean + errorEstimate

#  We can write this into a function this way
computeLambdaCI <- function(lambda, valid) {
  lowerBound <- lambda + ((1.96^2) / (2 * valid)) -
    sqrt(lambda * ((1.96^2) / valid) + ((1.96^4) / (4 * valid^2)))
  upperBound <- lambda + ((1.96^2) / (2 * valid)) +
    sqrt(lambda * ((1.96^2) / valid) + ((1.96^4) / (4 * valid^2)))
  return(data.frame(LB=lowerBound, UB=upperBound))
}
exampleCIs <- computeLambdaCI(thisLambda, 8000)
round(exampleCIs, 6) # LB: 0.675142, UB: 0.711633
round(exampleCIs * 8000) # suggests a 95% CI of 5401 to 5693 for the total count
#  We can do better by, for example, replacing 1.96^2 with the appropriate
#  distribution function in R. There are multiple ways to calculate the
#  confidence interval for a Poisson process (i.e., a thing in nature that looks
#  to us like a Poisson-type variable) because Poisson may be applied to
#  different situations. We will use the exact method here, as opposed to, for
#  example, the normal approximation method. This formulation relies on the
#  relationship between the chi-square distribution and the Poisson
#  distribution, as does the manual. There are also functions that have been
#  created by multiple R developers to compute Poisson confidence intervals. You
#  may want to look at the the poisson.test() function in base R and the
#  survival and epitools packages. For pedagogical purposes, we create our own
#  function here, for pedantic purposes there is a tangent in the addendum
#  about Poisson confidence intervals. One important note is that the manual
#  computes the confidence interval on lambda, which you will recall is the
#  negative natural log of the proportion of valid partitions. That is part of
#  the explanation for the square-root symbol and the powers in the calculation.
#  For ease, we calculate the confidence interval on the (target molecule)
#  count, so the result must be back-transformed to get the confidence interval
#  for lambda.
# 
#  computeLambdaCI() is a function that takes a molecule count value and a
#  desired confidence level and computes the confidence interval for a
#  Poisson-distributed random variable
computeLambdaCI <- function (count, conf.level=0.95) {
  alpha <- 1 - conf.level
  LB <- 0.5 * qchisq(alpha / 2, 2 * count)
  UB <- 0.5 * qchisq(1 - alpha / 2, 2 * count + 2)
  return(c(LB, UB))
}
computeLambdaCI(4000)

#  let's try to match this up with the Poisson distribution function, which we
#  read about earlier using ?dpois
moleculeCount <- 0:10
( countProbabilities <- dpois(moleculeCount, thisLambda) )
plot(moleculeCount, countProbabilities, type="h", bty="l",
     main=bquote(paste("Poisson density distribution for ",
                       lambda, " = ", .(thisLambda))),
     ylab="Molecule Count Probabilities", xlab="Count Per Partition")
#  This plot shows the probability (y-axis) of find each count on the x-axis in
#  one of the partitions. We can look at them like this
round(countProbabilities, 3)
#  This suggests that about half the time we will have 0 of the target molecule
#  in a partition, about 35% of the time we will have 1, a little over 10% will
#  have 2, about 2% will have 3, and a small fraction will have 4 or more. If we
#  wanted to be more literal, we could multiply by our number of partitions to
#  get the expected counts, as in
round(8000 * countProbabilities)
#  This is how QIAGEN is telling us to estimate the counts of the target
#  molecule


##----  Absolute quantification - copies per ml ----
#  The manual then provides the calculation for for the number of copies per
#  microliter based on the volume of each partition
# 
#  lambda_volume = lambda / volume in microliters
# 
#  The manual extends the example to the copies per ml calculation with
#  lambda as calculated previously, a volume of .34 nl is used for the 8k plate
V <- 0.34
( exampleCopiesPerML <- thisLambda / V * 1000 ) # 2038.668
#  The manual truncates this value without reporting it. Round will not work
#  here because it will round up to 2039
round(exampleCopiesPerML) # 2039
#  use floor() instead, this always round down to the nearest integer
floor(exampleCopiesPerML) # 2038

#  write this into a function for reusability
computeVolume <- function(lambda, volume) {
  return((lambda / volume) * 1000)
}

#  Next, the manual demonstrates how to calculate the copies of the target
#  molecule in the reaction volume as follows,
inputReactionVolume <- 12
( copiesInReaction <- computeVolume(thisLambda, V) * inputReactionVolume )

#  Let's amend our \code{computeVolume()} function to calculated the input
#  reaction volume if needed. Also, we'll start with a default volume of 
#  .91, which is what the 26000 nanoplate uses
computeVolume <- function(lambda, volume=.91, inputReactionVolume=F) {
  #  if the user specified an inputReactionVolume, return it with the reaction
  #  volume
  if(inputReactionVolume) {
    return(c((lambda / volume) * 1000,
             (lambda / volume) * 1000) * inputReactionVolume)
  }
  #  if the user didn't specify an input volume, return just the reaction volume
  return((lambda / volume) * 1000)
}
computeVolume(thisLambda, V, inputReactionVolume) # 2038 24456

#  The next section prior to concentration range reviews how to convert the
#  copies per ml to copies per ml sample. The manual reviews the example above,
#  explaining that the input reaction volume is calculated as:
#
#    5 microliter DNA sample
#  + 3 microliter 4x QIAcuity Probe Master Mix
#  + 4 microliter dPCR Primer/Probe Assay
#  = 12 microliter dPCR input reaction volume
# 
#  From this, we can calculate the number of copies of the target molecule per
#  microliter as follows
sampleVolume <- 5
( copiesPerML <- inputReactionVolume /
                  sampleVolume * computeVolume(thisLambda, V) ) # 4891.2
copiesPerML * sampleVolume # 24456

#  Next, the manual demonstrates how to calculate a percent confidence value
percentCV <- ((lambdaHigh.4 - lambdaLow.4) / (2 * thisLambda)) * 100 # 2.632024
#  we can format this as a percent and print it to the console as follows
print(paste0("The confidence value is: ", round(percentCV, 2), "%"))


#----  Concentration range ----
## The concentration range section reviews how to estimate the number of copies
#  of the target molecule per partition using the Poisson distribution.

#  Table 2. Percentage of expected target copies per partition for low, medium,
#  and high concentrations
#  This table shows where we expect the count values to be based on the
#  concentration. One way to recreate this table is as follows
table2.lambdas <- c(0.1, 0.5, 1, 1.5, 2, 5)
table2 <- data.frame(t(sapply(table2.lambdas, function(x) { dpois(0:10, x) })))
#  here we are using the sapply function to calculate the Poisson distribution
#  at each concentration in one function call. We then use the t() function to
#  transpose the result to match the manual.
#  Now store the column names, this is called left-handed assignment because the
#  function on the left calls the value of the column name and then stores the
#  values on the right.
colnames(table2) <- paste(0:10)
rownames(table2) <- paste(table2.lambdas)
print(round(table2, 4) * 100)

#  Note that this table does not match exactly. I suspect this is a rounding
#  error in the manual, given that the subsequent tables match.

#  Table 3. Expected number of partitions with different copies per partition
#  count, for 8500 partitions
table3 <- table2 * 8500
print(round(table3))

#  Table 4. Expected number of partitions with different copies per partition
#  count, for 26,000 partitions
table4 <- table2 * 26000
print(round(table4))

#  Write a function to plot the distribution of counts
plot.lambda <- function() {
  
}

##---- Conclusion ----
## This concludes our replication the QIAcuity User Manual Extension document.
#  I hope that this has been informative, and I welcome feedback and questions
#  at wthayer1@jh.edu

##---- References ----
#  QIAGEN (2021) QIAcuity: Nanoplate Digital PCR.
#     HB-2839-003_UM_QIAcuity_UM_Extension_0621_WW.pdf. www.qiagen.com

##---- Addendum ----
#  This is a tangent regarding Poisson confidence intervals, which is an
#  interesting area of research.
#  Here is how the survival package calculates an exact confidence interval for
#  the Poisson distribution, leaving out some of their error prevention and
#  other method goodness, k=number of successes, p=alpha
#  type this into the console to see: survival::cipoisson
#  and type this to read about it: ?survival::cipoisson
survExactPoisCI <- function(k, p=.05) {
  p <- p / 2
  dummy1 <- ifelse(k == 0, 1, k)
  lower <- ifelse(k == 0, 0, qgamma(p, dummy1))
  upper <- qgamma(1 - p, k + 1)
  return(c(lower, upper))
}
survExactPoisCI(computeLambda(8000, 4000))
computeLambdaCI(thisLambda) == survExactPoisCI(thisLambda) # FALSE TRUE
#  let's try the epitools package, ?epitools::pois.exact
epitools::pois.exact(4000)
#  from the exactci package
exactci::poisson.exact(4000)
exactci::exactpoissonPlot(4000)

#  from McMaster University, "Exact" 95% Confidence Intervals
#  https://ms.mcmaster.ca/peter/s743/poissonalpha.html
#  last accessed March 14, 2023
#  This is an interesting visualization of exact 95% confidence intervals for
#  the Poisson distribution, and a comparison to Pearson confidence intervals
#  For each observation x, calculate the 95% CI for mu as follows
x <- seq(0, 20, .1)
mu <- seq(0, 20, .1)
round(cbind(( qchisq(0.025, 2 * x) / 2 ),
            ( qchisq(0.975, 2 * (x + 1)) / 2 ) ), 4)
#  The probability each interval will miss mu, if xPrime satisfies
#    qchisq(0.975, 2*(x'+1))/2 < µ < qchisq(0.975, 2*(x'+2))/2
#  which is
cbind( qchisq(0.975, 2 * (x + 1)) / 2,
      ( qchisq(0.975, 2 * (x + 2)) / 2),
      mu)
#  then the probability of a miss on the right is the Poisson probability that
#  x <= x`
#  for the lower limit, if xPrime satisfies
#    qchisq(0.025, 2*(x'-1))/2 < µ < qchisq(0.025, 2*x')/2
#  which is
round(cbind( qchisq(0.025, 2 * (xPrime - 1)) / 2,
             qchisq(0.025, 2 * xPrime) / 2
             , mu), 4)
#  then the probability the confidence interval will miss on the left is the
#  Poisson probability x >= x`
#  to calculate these, use
xPrime.right <- qpois(0.025, mu, F)
xPrime.left <- qpois(0.975, mu, F)
leftMiss <- ppois(xPrime.right, mu, F)
rightMiss <- 1 - ppois(xPrime.left - 1, mu, F)
totalMiss <- rightMiss + leftMiss
plot(mu, rightMiss, type="l", col="blue", bty="l", ylab=expression(alpha),
     xlab=expression(mu), ylim=c(0, .05))
lines(mu, leftMiss, type="l", col="red")
lines(mu, totalMiss, col="green4")
abline(h=c(.025, .05), col="deepskyblue")
#  reproduce with mu up to 60
mu <- seq(0, 60, .1)
xPrime.right <- qpois(0.025, mu, F)
xPrime.left <- qpois(0.975, mu, F)
rightMiss <- ppois(xPrime.right, mu, F)
leftMiss <- 1 - ppois(xPrime.left - 1, mu, F)
totMiss <- rightMiss + leftMiss
plot(mu, rightMiss, type="l", col="red", bty="l", ylab="alpha",
     ylim=c(0, .05))
lines(mu, leftMiss, type="l", col="blue")
lines(mu, totMiss, col="green4")
abline(h=c(.025, .05), col="lightblue")
#  Now try with Pearson confidence intervals
#  The Pearson confidence limits are the roots of the quadratic equation
#    ( (x - mu)^2) / mu = a
#  which can shown to be
#    (x + a / 2 ) + or - sqrt(a)*sqrt(x + a/4)
#  where a = qchisq(.95, 1)
##! This isn't correct ##!
a <- qchisq(.95, 1)
round(cbind((x + a / 2 ) - sqrt(a) * sqrt(x + a / 4),
            (x + a / 2 ) + sqrt(a) * sqrt(x + a / 4)), 4)
xPrime.left <- qchisq(.975, mu)
xPrime.right <- qchisq(.025, mu)
rightMiss <- ppois(xPrime.left, mu, F)
leftMiss <- 1 - ppois(xPrime.right - 1, mu, F)
totalMiss <- rightMiss + leftMiss
plot(leftMiss, type="l", col="red", bty="l", ylim=c(0, .05))
lines(rightMiss, col="blue")
lines(totalMiss, col="green4")
abline(h=c(.025, .05), col="deepskyblue")

#  Holladay (2014) reviews confidence intervals in a Master's thesis, including
#  Wald, Garwood, and the Scores method, before proposing his own procedure.
#  recreate his example here.
holl.wald <- function(x) {
  lx <- x - qnorm(.975) * sqrt(x)
  ux <- x + qnorm(.975) * sqrt(x)
  return(cbind(lx, ux))
}
round(holl.wald(18), 2) # 9.68 26.32
#  Garwood
holl.garwood <- function(x) {
  lx <- .5 * qchisq(.025, 2*x)
  ux <- .5 * qchisq(.975, 2*(1+x))
  return(cbind(lx, ux))
}
round(holl.garwood(18), 2)
#  Scores
holl.scores <- function(x) {
  lx <- x + .5 * qnorm(.975)^2 - qnorm(.975) * sqrt(x + .25 * qnorm(0.975)^2)
  ux <- x + .5 * qnorm(.975)^2 + qnorm(.975) * sqrt(x + .25 * qnorm(.975)^2)
  return(cbind(lx, ux))
}
round(holl.scores(18), 2)
footTraffic <- data.frame(Day=1:14, Traffic=c(86, 96, 129, 146, 107, 138, 95,
                                              90, 132, 149, 118, 143, 125, 100))
mean(footTraffic$Traffic) * 14 # average foot traffic in 2-week period: 1654
holl.wald(sum(footTraffic$Traffic))
holl.garwood(sum(footTraffic$Traffic))
holl.scores(sum(footTraffic$Traffic))
#  these match, try them with the QIAcuity :registered: example
holl.wald(4000)
holl.garwood(4000)
holl.scores(4000)
#  this is a simplistic approximation of the confidence interval, that might be
#  appropriate with large sample sizes
thisLambda - 1.96 * sqrt(thisLambda / 8000)
thisLambda + 1.96 * sqrt(thisLambda / 8000)


#  Because we know it is related, let's look at the chi-square distribution,
#  just for fun
#  store a large random sample of chi-squre values, then plot it and add a curve
x <- rchisq(10000, df = 8000)
hist(x, breaks = 'Scott', freq = FALSE, xlab = '',
     main = expression(paste(chi^2,
                             "-distribution with 8000 degrees of freedom")))
curve(dchisq(x, df=8000), from=7000, to=9000, n=5000,
      col='yellow3', lwd=2, add=T)

x <- rpois(1000, mu[1])
