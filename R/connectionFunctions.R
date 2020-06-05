########################################################
# Basic R function for the package.
# Author: Mathieu Fortin, Canadian Forest Service, Canadian Wood Fibre Centre
# Copyright; Her Majesty the Queen in right of Canada
# Date: April 2019
########################################################

.welcomeMessage <- function() {
  packageStartupMessage("Welcome to betadiv package! This package implements estimators of beta diversity indices.")
  packageStartupMessage("For more information, visit https://sourceforge.net/p/mrnfforesttools/divindices/wiki/Home/ .")
}


.onAttach <- function(libname, pkgname) {
  .welcomeMessage()
}

# .onUnload <- function(libpath) {
#   shutdownJava()
# }
#
# .onDetach <- function(libpath) {
#   shutdownJava()
# }


.addToArray <- function(refArray, array) {
  if (length(refArray) != length(array)) {
    stop("Incompatible array length!")
  } else {
    for (i in 1:length(array)) {
      refArray[[i]] <- c(refArray[[i]], array[[i]])
    }
  }
  return(refArray)
}

.convertJavaDataSetIntoDataFrame <- function(dataSetObject) {
  refArray <- NULL
  observations <- J4R::callJavaMethod(dataSetObject, "getObservations")
  observations <- J4R::getAllValuesFromListObject(observations)
  for (obs in observations) {
    array <- J4R::callJavaMethod(obs, "toArray")
    array <- as.list(J4R::getAllValuesFromArray(array))
    if (is.null(refArray)) {
      refArray <- array
    } else {
      refArray <- .addToArray(refArray, array)
    }
  }
  dataFrame <- NULL
  for (i in 1:length(refArray)) {
    dataFrame <- as.data.frame(cbind(dataFrame, refArray[[i]]))
  }
  colnames(dataFrame) <- J4R::getAllValuesFromListObject(J4R::callJavaMethod(dataSetObject, "getFieldNames"))
  return(dataFrame)
}


.loadREpicea <- function() {
  if (!J4R::checkIfClasspathContains("repicea.jar")) {
    J4R::addUrlToClassPath("repicea.jar", packageName = "betadiv")
  }
}

.connectToBetadivLibrary <- function(memSize = NULL) {
  if (!J4R::isConnectedToJava()) {
    J4R::connectToJava(memorySize = memSize)
  }
  .loadREpicea()
  .loadBetadiv()
}


.loadBetadiv <- function() {
  if (!J4R::checkIfClasspathContains("betadiversityindices.jar")) {
    J4R::addUrlToClassPath("betadiversityindices.jar", packageName = "betadiv")
  }
}

