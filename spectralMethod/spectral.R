rm(list = ls())
setwd("/Users/yang.2677/Desktop/SBM/spectralMethod")

## Select your own path in which spectral function is located
source("spectralFun.r")
require(grid)

## Step1: generate a random graph and membership assignment
nN = 100	# node number
k = 2	# community number
P = matrix(c(.8,.3,.3,.8), ncol = 2) # Group structure matrix
Amem = DCSBM(nN, k = 2, P,random.community.assignment=TRUE) 
A = Amem$Adjacency	# Adjacency matrix
EA = Amem$Expected	# Expectation of adjacency matrix
mem = Amem$Membership	# Membership assignment

## Step2: calculate laplacian of adjacency matrix
L = spectral.Laplacian(EA, 0)

## Step3: run kmeans and print output
out = spectral.kmeans(L)
X = out$EigenVecMat
mem.opt = out$mem.opt
printOut(mem, mem.opt)
