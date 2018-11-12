## Spectral method for SBM, written by Xige Yang

# generate adjacency matrix
DCSBM = function(n, k = 2, P, sizes = c(round(n / 2), n - round(n / 2)), random.community.assignment = c(FALSE, TRUE), delta = rep(0, k), 
edge.list = c(FALSE, TRUE), theta = rep(1, n)){

  #Check that entries are appropriate
  if(sum(sizes) != n){
    stop("argument sizes must sum to n")
  }
  if(length(sizes) != k){
    stop("argument sizes must be of length k")
  }

  #Create community labels
  Y = matrix(rep(0, n*k), ncol = k)
  index = list()
  possible = 1:n
  Membership = rep(1, n)
  random.community.assignment = random.community.assignment[1]

  #assign vertex labels randomly if random.assignment == TRUE
  if(random.community.assignment == TRUE){
    for(i in 1:k){
      index[[i]] = sample(possible, sizes[i])
      Membership[index[[i]]] = i
      possible = setdiff(possible, index[[i]])
    }
  }

  #assign vertex labels in order if random.assignment == FALSE
  if(random.community.assignment == FALSE){
    for(i in 1:k){
      index[[i]] = possible[1:sizes[i]]
      Membership[index[[i]]] = i
      possible = setdiff(possible, index[[i]])
    }
  }
  for(i in 1:k){
    Y[index[[i]], i] = 1
  }

  edge.list = edge.list[1]

  #Generate expected value of the adjacency matrix
  Theta = diag(theta)
  expected.A = Theta%*%Y%*%P%*%t(Y)%*%Theta

  #Generate a realization of the adjacency matrix via Poisson draws
  foo = matrix(rbinom(n^2, 1, matrix(expected.A, ncol = 1)), ncol = n)
  Adj = matrix(0, ncol = n, nrow = n)
  Adj[upper.tri(Adj)] = foo[upper.tri(foo)]
  Adj = Adj + t(Adj)
  diag(Adj) = diag(foo)

  if(edge.list == TRUE){
    return(list(Adjacency = as.vector(Adj), Membership = Membership))
  }

  if(edge.list == FALSE){
    return(list(Adjacency = Adj, Expected = expected.A, Membership = Membership))
  }
}

## Initialize
IniLize = function(A, memini, K = 2){
	N = ncol(A)
	
	# local variables
	nVec = rep(0, K)
	O = matrix(0, K, K)
	P = matrix(0, K, K)
	Lam = matrix(0, K, K)
	
	# get n-Matrix
	for(k in 1:K){
		nVec[k] = sum(as.numeric(memini == k))
	}
	nVec = as.matrix(nVec)
	nMat = nVec%*%t(nVec)
	diag(nMat) = nVec*(nVec-1)
	
	# get O-Matrix
	for(k in 1:K){
		fook = as.numeric(mem == k)
		for(l in 1:K){
			fool = as.numeric(mem == l)
			O[k,l] = fool%*%A%*%fook
		}
	}
	
	# get initial evaluations
	pai = nVec/N
	R = diag(as.numeric(pai))
	P = O/nMat
	for(l in 1:K){
		for(k in 1:K){
			if(nMat[l,k] == 0){
				stop(c("nMat[", l, ',', k, "] = 0!"))
			}
			Lam[l,k] = N*R[k,]%*%P[,l]
		}
	}
	return(list(Omat = O, Pmat = P, LamMat = Lam, nMat = nMat, pai = pai))
}

## Parameter update
ParaUpdt = function(A, O, P, Lam, nMat, pai, b, mem, memini, K = 2, iterMax = 100, eps = .01){
	# initialize
	N = ncol(A)
	paiMat = matrix(0, N, K)
	paiUpdt = pai
	LamUpdt = Lam
	paiOri = pai + 1
	LamOri = Lam + 1
	memopt = memini
	iterCt = 0
	
	# substep2-4: update paiMat, pai and Lam
	while(iterCt <= iterMax && norm(as.matrix(LamOri-LamUpdt,'F')) > eps && norm(as.matrix(paiOri-paiUpdt),'F') > eps){
		paiOri = paiUpdt
		LamOri = LamUpdt
		iterCt = iterCt + 1
		for(i in 1:N){
			for(l in 1:K){
				prod = 1
				for(m in 1:K){
					prod = prod*exp(b[i,m]*log(Lam[l,m]) - Lam[l,m])
				}
				paiMat[i,l] = pai[l]*prod
			}
		}
		for(i in 1:N){
			foo = sum(paiMat[i,])
			for(l in 1:K){
				paiMat[i,l] = paiMat[i,l]/foo
			}
		}
	
		for(l in 1:K){
			paiUpdt[l] = sum(paiMat[,l])/N
			if(paiUpdt[l] == 0){
				stop(c("paiUpdt[", l, "] = 0!"))
			}
		}
		for(l in 1:K){
			foo = sum(paiMat[,l])
			for(k in 1:K){
				LamUpdt[l,k] = sum(paiMat[,l]*b[,k])/foo
			}
		}
	}
	if(iterCt >= iterMax){
		warning("Iteration count exceeds upper bound.")
	}
	
	# substep5: update mem assignment
	for(i in 1:N){
		memopt[i] = which(paiMat[i,] == max(paiMat[i,]))[1]
	}
	
	return(list(pai = paiUpdt, paiMat = paiMat, Lam = LamUpdt, memopt = memopt))
}

# Print outcomes
printOut = function(mem, memopt, P, PEst){
	cat("True membership assignment at mem = ", mem, '\n')
	cat("Optimum membership assignment at mem = ", memopt, '\n')
	cat("Hamming distance between true and estimated assignment = ", hamming.distance(mem, memopt), '\n')
	cat("True P value = \n", P, '\n')
	cat("Estimated P = \n", PEst)
}