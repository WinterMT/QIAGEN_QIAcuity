## QIAcuity User Manual Extension QIAGEN: Statistics of nanoplate dPCR
## (c) 2023 Winter Maxwell Thayer, Johns Hopkins University,
##        all rights reserved.
## 2022-04-13


##---- QIAcuity User Manual Extension document ----
#  HB-2839-003_UM_QIAcuity_UM_Extension_0621_WW first section Nanoplate Digital
#  PCR contains subsection Statistics of nanoplate dPCR on page 7. This R file
#  describes appropriate procedures for replicating all analyses contained
#  therein in the R statistical environment. Note that you may want to run the
#  function rm(list=ls()) to clear your environment prior to running this
#  script. Be sure no valuable data is stored in your R environment first, maybe
#  by running ls() and looking at the output, just to be safe.


##---- Poisson distribution ----
#  The first portion of the manual describes the natural distribution of the
#  dPCR quantification and suggests the use of the theoretical distribution
#  Poisson to approximate the copies of the target molecule per positive
#  partition. The manual provides an example of a theoretical Poisson
#  distribution with lambda ranging from 0.1 to 5. The formula for Poisson
#  is described and the following parameters are defined:
# 
#  e: Euler's constant, can be obtained with \code{exp(1)} in R
#  Lambda: can be thought of as the concentration of the target molecule, or,
#            the average count per microliter
#  k: the copies of the molecule per partition
#  k!: k factorial, obtained with \code{factorial(k)} in R
#
#  plots of the theoretical Poisson distribution are then presented, replicate
#  these plots with the following code
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
  #  clipping area, ie., where R stops drawing on your plot.
  if(i == 1) {
    par(xpd=NA) # clip to the device region
    text(-5, .3, labels="Probability P(X=k)", srt=90, adj=0, xpd=NA)
    par(xpd=F) # reset to clip to plot regoion
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
#  lambda values. Note that as lambda increases, the distribution of the counts
#  becomes increasingly normal.
# 
#  The last paragraph of this section describes the calculation of copies per
#  partition, which is:
#
#    copies of molecule per partition * number of valid partitions
#
#  The following section describes how to calculate the copies of the molecule
#  per partition.


##----  Absolute quantification - copies per partition ----
#  this section of the manual describes how to calculate lambda, the average
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
#  I suggest waiting to round until the final step.

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
#  distribution to express the confidence interval than in the previous section.
#  As of this writing, the Poisson distribution entry for Wikipedia 
#  parameterizes the distribution this way for the confidence interval as well:
# 
#  (1/2)*chiSquare(alpha/2, 2*k) ≤ Mu ≤ (1/2)*chiSquare(1 - (alpha/2), 2*k + 2)
#
#  This expression relies on the relationship between the Poisson,
#  Chi-square, and Gamma distributions. In this application, it breaks down to,
#  we expect the average number of copies per partition to be between our
#  estimate of lambda (the average number of copies per partition) plus or minus
#  our estimate of the amount of variability in lambda based on what we know
#  about the theoretical distribution with this many degrees of freedom, in this
#  case 8000. Do not worry about degrees of freedom right now.
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

#  to make that more readable, break it into smaller parts and give names
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

#  To make this code reusable, we should write a function.
computeLambdaCI <- function(lambda, valid) {
  lowerBound <- lambda + ((1.96^2) / (2 * valid)) -
    sqrt(lambda * ((1.96^2) / valid) + ((1.96^4) / (4 * valid^2)))
  upperBound <- lambda + ((1.96^2) / (2 * valid)) +
    sqrt(lambda * ((1.96^2) / valid) + ((1.96^4) / (4 * valid^2)))
  return(data.frame(LB=lowerBound, UB=upperBound))
}
exampleCIs <- computeLambdaCI(thisLambda, 8000)
round(exampleCIs, 6) # LB: 0.675142, UB: 0.711633


##----  Absolute quantification - copies per ml ----
#  The manual then provides the calculation for for the number of copies per
#  microliter based on the volume of each partition
# 
#  lambda_volume = lambda / volume in microliters
# 
#  The manual extends the example to the copies per ml calculation with
#  lambda as calculated previously, a volume of .34 nl
V <- 0.34
( exampleCopiesPerML <- thisLambda / V * 1000 ) # 2038.668
#  The manual truncates this value without reporting it. Round will not work
#  here because it will round up to 2039
round(exampleCopiesPerML) # 2039
#  use floor() instead, this always round down to the nearest integer
floor(exampleCopiesPerML) # 2038

#  write this into a function for reusability
computeVolume <- function(lambda, volume) {
  return((lambda/volume) * 1000)
}

#  Next, the manual demonstrates how to calculate the copies of the target
#  molecule in the reaction volume as follows,
inputReactionVolume <- 12
( copiesInReaction <- computeVolume(thisLambda, V) * inputReactionVolume )

#  Let's amend our \code{computeVolume()} function to calculated the input
#  reaction volume if needed.
computeVolume <- function(lambda, volume, inputReactionVolume=F) {
  #  if the user specified an inputReactionVolume, return it with the reaction
  #  volume
  if(inputReactionVolume) {
    return(c((lambda/volume) * 1000,
             (lambda/volume) * 1000) * inputReactionVolume)
  }
  #  if the user didn't specify an input volume, return just the reaction volume
  return(floor((lambda/volume) * 1000))
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

#  Table 3. Expected number of partitions with different copies per partition
#  count, for 8500 partitions
table3 <- table2 * 8500
print(round(table3))

#  Table 4. Expected number of partitions with different copies per partition
#  count, for 26,000 partitions
table4 <- table2 * 26000
print(round(table4))


##---- Conclusion ----
## This concludes our replication the QIAcuity User Manual Extension document.
#  I hope that this has been informative, and I welcome feedback and questions
#  at wthayer1@jh.edu