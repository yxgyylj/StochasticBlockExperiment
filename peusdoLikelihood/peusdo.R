rm(list = ls())
setwd("/Users/yang.2677/Desktop/SBM/peusdoLikelihood")

## Select your own path in which spectral function is located
source("peusdoFun.r")
library("e1071")
require(grid)

## Step1: Generate a random graph and membership assignment
N = 100	# node number
K = 2	# community number
P = matrix(c(.8,.3,.3,.8), ncol = 2) # Group structure matrix
Amem = DCSBM(N, k = 2, P,random.community.assignment=FALSE) 
A = Amem$Adjacency	# Adjacency matrix
EA = Amem$Expected	# Expectation of adjacency matrix
mem = Amem$Membership	# True membership assignment

## Step2: Initialize data
step = 10 # total iterations
memini = c(rep(1,floor(N/2)), rep(2,N-floor(N/2))) # initial guess
b = matrix(0, N, K)
IniData = IniLize(A, memini)
O = IniData$Omat
Pmat = IniData$Pmat
Lam = IniData$LamMat
nMat = IniData$nMat
pai = IniData$pai
PEst = matrix(0, K, K)

## Step3: Run algorithm
for(i in 1:step){
	# substep1: update b
	for(k in 1:K){
		b[,k] = as.numeric(memini == k)%*%A
	}
	Data = ParaUpdt(A, O, Pmat, Lam, nMat, pai, b, mem, memini)
	pai = Data$paiUpdt
	Lam = Data$Lamupdt
}
paiMat = Data$paiMat
memopt = Data$memopt

## Step4: Finally update P
for(l in 1:K){
	for(k in 1:K){
		PEst[l,k] = (t(paiMat)%*%A%*%paiMat)/nMat[l,k]
	}
}
printOut(mem, memopt, P, PEst)
image(A, col = paste("gray", 99:1, seq = ""))
