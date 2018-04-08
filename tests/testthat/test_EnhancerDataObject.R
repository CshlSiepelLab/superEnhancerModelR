context("superEnhancerDataObject")

## Create a test design data.frame
design=expand.grid(E1=c(0,1),E2=c(0,1),E3=c(0,1))

## Create fake expression data
expression=c(0.001,0.2,0.3,0.7,0.4,0.8,0.6,1)*100

## Create activity function
actFun=formula(~E1+E2+E3+E1:E2)

## Test that object is created correctly
testthat::test_that("enhancerDataObject builds correctly",{
  testthat::expect_s4_class(enhancerDataObject(expression,design,actFun,linkFunction = "logistic"),"enhancerDataObject")
})

