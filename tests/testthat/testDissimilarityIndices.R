############################################
# Simple tests for dissmilarity indices
# Mathieu Fortin - April 2019
############################################

library(betadiv)

dataReleves <- betadiv::subsetUrbanEnvironmentNancy

strataList <- unique(dataReleves$Stratum)

output <- NULL
baselga <- NULL
stratum <- strataList[1]

for (stratum in strataList) {
  releve.s <- dataReleves[which(dataReleves$Stratum == stratum),]
  if (stratum == "forest") {
    populationSize <- 3089 * 10000 / (pi * 5^2)
  } else if (stratum == "parking") {
    populationSize <- 501 * 10000 / (pi * 5^2)
  } else {
    populationSize <- 100000
  }

  indices <- getDissimilarityEstimates(releve.s, "CODE_POINT", "Espece", populationSize, memSize = 500)
  indices$stratum <- stratum
  output <- rbind(output, indices)
}

test_that("Testing forests", {
  expect_equal(abs(output[which(output$stratum == "forest"), "Simpson"] - 0.2603477) < 1E-6, TRUE)
  expect_equal(abs(output[which(output$stratum == "forest"), "stdErrSimpson"] - 0.02001768) < 1E-8, TRUE)
  expect_equal(abs(output[which(output$stratum == "forest"), "Sorensen"] - 0.4025105) < 1E-6, TRUE)
  expect_equal(abs(output[which(output$stratum == "forest"), "stdErrSorensen"] - 0.01579914) < 1E-8, TRUE)
  expect_equal(abs(output[which(output$stratum == "forest"), "Nestedness"] - 0.1421628) < 1E-6, TRUE)
  expect_equal(abs(output[which(output$stratum == "forest"), "stdErrNestedness"] - 0.02596182) < 1E-8, TRUE)
  expect_equal(abs(output[which(output$stratum == "forest"), "Alpha"] - 20.42857) < 1E-5, TRUE)
  expect_equal(abs(output[which(output$stratum == "forest"), "stdErrAlpha"] - 4.110713) < 1E-5, TRUE)
  expect_equal(abs(output[which(output$stratum == "forest"), "Gamma"] - 138.8225) < 1E-4, TRUE)
  expect_equal(abs(output[which(output$stratum == "forest"), "stdErrGamma"] - 27.33917) < 1E-4, TRUE)

  expect_equal(abs(output[which(output$stratum == "parking"), "Simpson"] - 0.2996533) < 1E-6, TRUE)
  expect_equal(abs(output[which(output$stratum == "parking"), "stdErrSimpson"] - 0.03163888) < 1E-8, TRUE)
  expect_equal(abs(output[which(output$stratum == "parking"), "Sorensen"] - 0.4675060) < 1E-6, TRUE)
  expect_equal(abs(output[which(output$stratum == "parking"), "stdErrSorensen"] - 0.01100832) < 1E-8, TRUE)
  expect_equal(abs(output[which(output$stratum == "parking"), "Nestedness"] - 0.1678526) < 1E-6, TRUE)
  expect_equal(abs(output[which(output$stratum == "parking"), "stdErrNestedness"] - 0.04033552) < 1E-8, TRUE)
  expect_equal(abs(output[which(output$stratum == "parking"), "Alpha"] - 13.09091) < 1E-5, TRUE)
  expect_equal(abs(output[which(output$stratum == "parking"), "stdErrAlpha"] - 2.986194) < 1E-5, TRUE)
  expect_equal(abs(output[which(output$stratum == "parking"), "Gamma"] - 181.8853) < 1E-4, TRUE)
  expect_equal(abs(output[which(output$stratum == "parking"), "stdErrGamma"] - 33.16738) < 1E-4, TRUE)
})

J4R::shutdownJava()

