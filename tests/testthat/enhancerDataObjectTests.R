context("superEnhancerDataObject")

## Create a test design matrix
design=matrix(c(0,1,0,1,0,0,1,1),nrow=4)
colnames(design)=c("E1","E2")
design=as.data.frame(design)

## Create fake expression data
expression=c(0.001,0.2,0.3,0.9)*100

## Create activity function
actFun=formula(~E1+E2+E1*E2)

## Test that object is created correctly
test_that("enhancerDataObject builds correctly",{
  expect_s4_class(enhancerDataObject(expression,design,actFun,linkFunction = "logistic"),"enhancerDataObject")
})
