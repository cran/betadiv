########################################################
# R functions to estimate dissimilarity indices.
# Author: Mathieu Fortin, Canadian Forest Service, Canadian Wood Fibre Centre
# Copyright; Her Majesty the Queen in right of Canada
# Date: April 2019
########################################################


#'
#' Data from Urban Environments in Nancy, France
#'
#' A subset of a two-phase inventory that was carried out
#' in Nancy, Eastern France, in 2015
#'
#' @docType data
#'
#' @usage data(subsetUrbanEnvironmentNancy)
#'
#' @keywords datasets
#'
#' @examples
#' data(subsetUrbanEnvironmentNancy)
"subsetUrbanEnvironmentNancy"



.createDissimilarityIndicesEstimator <- function(memSize = NULL) {
  .connectToBetadivLibrary(memSize)
#  print("Instantiating the estimator of dissimilarity indices...")
  msi <- J4R::createJavaObject("biodiversity.indices.MultipleSiteIndex")
#  print("Done.")
  return(msi)
}


.createSample <- function(dataSet, plotIdField, speciesIdField) {
#  .connectToBetadivLibrary()
  sample <- J4R::createJavaObject("java.util.HashMap")

  for (plotId in unique(dataSet[,plotIdField])) {

    releve.i <- dataSet[which(dataSet[,plotIdField] == plotId),]
    java.plot <- J4R::createJavaObject("java.util.ArrayList")
    speciesListBefore <- releve.i[,speciesIdField]
    speciesListAfter <- speciesListBefore[which(!is.na(speciesListBefore))]
    if (length(speciesListAfter) < length(speciesListBefore) && length(speciesListBefore) > 1) {
      stop(paste("There seems to be some NA in the species list of plot ", plotId, sep=""))
    } else {
      if (length(speciesListAfter) > 0) {
        J4R::callJavaMethod(java.plot, "add", as.character(speciesListAfter))
      } else {
        print(paste("This plot is empty ", plotId, sep=""))
      }
      J4R::callJavaMethod(sample, "put", as.character(plotId), java.plot)
    }
  }
  return(sample)
}



#'
#' Estimate Dissimilarity Indices
#'
#' This function computes estimates of the adapted dissimilarity indices of Simpson, Sorensen and
#' nestedness from a sample.
#'
#' The dissimilarity indices were adapted from those of Baselga (2010). These adapted indices are population size independent
#' so that it is possible to compare the dissimilarity of two populations of unequal sizes.
#'
#' This function implements estimators of these adapted indices. A sample of plots with species observations must be passed to the function
#' as well as the population size, that is the number of plots that fit in this population. The variance estimation is based on the Jackknife
#' method. The function returns a data.frame object with the estimates of the multiple-site version of Simpson, Sorensen and nestedness as well
#' as their associated standard errors. In addition, the function also provides an estimate of the alpha and gamma diversity. The
#' gamma diversity estimate is based on the Chao2 estimator (Chao and Lin 2012).
#'
#' @param dataset a data.frame object that contains at least two fields: one for the sample plot ids and the other
#' for the species. Each row is actually an observation of a species in a particular plot.
#' @param plotIdField the name of the field that contains the sample plot id in the dataset.
#' @param speciesIdField the name of the field that contains the species in the dataset.
#' @param populationSize the number of units in the population. That is the total number of sample plots that
#' could fit in the population. Under the assumption that the plot size is constant, the population size is calculated
#' as the area of the population divided by the area of a single sample plot.
#' @param memSize the size of the Java Virtual Machine in Mg (if not specified the JVM is instantiated with the default memory size, which depends on the available RAM)
#'
#' @return a data.frame object with the estimated dissimilarity indices and their standard errors
#'
#' @examples
#'
#' ### An example using the subsetUrbanEnvironmentNancy dataset ###
#'
#' dataReleves <- betadiv::subsetUrbanEnvironmentNancy
#'
#' strataList <- unique(dataReleves$Stratum)
#' output <- NULL
#' baselga <- NULL
#' stratum <- strataList[1]
#'
#' for (stratum in strataList) {
#'   releve.s <- dataReleves[which(dataReleves$Stratum == stratum),]
#'   if (stratum == "forest") {
#'     populationSize <- 3089 * 10000 / (pi * 5^2)
#'   } else if (stratum == "parking") {
#'     populationSize <- 501 * 10000 / (pi * 5^2)
#'   } else {
#'     populationSize <- 100000
#'   }
#'
#'   indices <- getDissimilarityEstimates(releve.s, "CODE_POINT", "Espece",
#'                                          populationSize, memSize = 500)
#'   indices$stratum <- stratum
#'   output <- rbind(output, indices)
#' }
#'
#' @references Fortin, M., A. Kondratyeva, and R. Van Couwenberghe. 2020. Improved Beta-diversity estimators
#' based on multiple-site dissimilarity: Distinguishing the sample from the population. Global Ecology and
#' Biogeography 29: 1073-1084. \url{https://doi.org/10.1111/geb.13080}
#'
#' Baselga, A. 2010. Partitioning the turnover and nestedness components of beta diversity.
#' Global Ecology and Biogeography 19:134-143.
#'
#' Chao A., and C.-W Lin. 2012. Nonparametric lower bounds for species richness and shared species richness
#' under sampling without replacement. Biometrics 68: 912-921.
#'
#'
#' @export
getDissimilarityEstimates <- function(dataset, plotIdField, speciesIdField, populationSize, memSize = NULL) {
  dissimilarityEstimator <- .createDissimilarityIndicesEstimator(memSize)
  sample <- .createSample(dataset, plotIdField, speciesIdField)
  n <- J4R::callJavaMethod(sample, "size")
  messageToBeDisplayed <- paste("Estimating dissimilarity from a sample of", n, "plots...")
  if (n > 500) {
    messageToBeDisplayed <- paste(messageToBeDisplayed, "This may take a while!")
  }
  message(messageToBeDisplayed)
  leaveOneOutMode <- J4R::createJavaObject("biodiversity.indices.MultipleSiteIndex$Mode", "LeaveOneOut")
  dissimilarityEstimates <- J4R::callJavaMethod(dissimilarityEstimator, "getDissimilarityIndicesMultiplesiteEstimator", sample, as.integer(populationSize), TRUE, leaveOneOutMode)

  simpsonEnum <- J4R::createJavaObject("biodiversity.indices.DiversityIndices$BetaIndex", "Simpson")
  sorensenEnum <- J4R::createJavaObject("biodiversity.indices.DiversityIndices$BetaIndex", "Sorensen")
  nestednessEnum <- J4R::createJavaObject("biodiversity.indices.DiversityIndices$BetaIndex", "Nestedness")

  alphaEstimate <- J4R::callJavaMethod(dissimilarityEstimates, "getAlphaDiversity")
  gammaEstimate <- J4R::callJavaMethod(dissimilarityEstimates, "getGammaDiversity")
  simpsonEstimate <- J4R::callJavaMethod(dissimilarityEstimates, "getBetaDiversity", simpsonEnum)
  sorensenEstimate <- J4R::callJavaMethod(dissimilarityEstimates, "getBetaDiversity", sorensenEnum)
  nestednessEstimate <- J4R::callJavaMethod(dissimilarityEstimates, "getBetaDiversity", nestednessEnum)

  Alpha <- J4R::callJavaMethod(J4R::callJavaMethod(alphaEstimate, "getMean"), "getSumOfElements")
  varAlpha <- J4R::callJavaMethod(J4R::callJavaMethod(alphaEstimate, "getVariance"), "getSumOfElements")
  stdErrAlpha <- varAlpha^.5

  Gamma <- J4R::callJavaMethod(J4R::callJavaMethod(gammaEstimate, "getMean"), "getSumOfElements")
  varGamma <- J4R::callJavaMethod(J4R::callJavaMethod(gammaEstimate, "getVariance"), "getSumOfElements")
  stdErrGamma <- varGamma^.5

  Simpson <- J4R::callJavaMethod(J4R::callJavaMethod(simpsonEstimate, "getMean"), "getSumOfElements")
  varSimpson <- J4R::callJavaMethod(J4R::callJavaMethod(simpsonEstimate, "getVariance"), "getSumOfElements")
  stdErrSimpson <- varSimpson^.5

  Sorensen <- J4R::callJavaMethod(J4R::callJavaMethod(sorensenEstimate, "getMean"), "getSumOfElements")
  varSorensen <- J4R::callJavaMethod(J4R::callJavaMethod(sorensenEstimate, "getVariance"), "getSumOfElements")
  stdErrSorensen <- varSorensen^.5

  Nestedness <- J4R::callJavaMethod(J4R::callJavaMethod(nestednessEstimate, "getMean"), "getSumOfElements")
  varNestedness <- J4R::callJavaMethod(J4R::callJavaMethod(nestednessEstimate, "getVariance"), "getSumOfElements")
  stdErrNestedness <- varNestedness^.5

  return(data.frame(n, Simpson, stdErrSimpson,
                    Sorensen, stdErrSorensen,
                    Nestedness, stdErrNestedness,
                    Alpha, stdErrAlpha,
                    Gamma, stdErrGamma))
}



