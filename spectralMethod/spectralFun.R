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
  diag(Adj) = 0

  if(edge.list == TRUE){
    return(list(Adjacency = as.vector(Adj), Membership = Membership))
  }

  if(edge.list == FALSE){
    return(list(Adjacency = Adj, Expected = expected.A, Membership = Membership))
  }
}

## Spectral matrix
# Laplacian of adjacency matrix
spectral.Laplacian = function(A,tau){
	dimA = dim(A)
	D = diag(rep(tau,dimA[1]))
	L = D
	
	for(i in 1:dimA[1]){
		D[i,i] = D[i,i] + sum(A[i,])
		if(D[i,i] == 0){
			stop("D is singular!")
		}
	}
	D = sqrt(D)
	
	for(i in 1:dimA[1]){
		for(j in 1:dimA[2]){
			L[i,j] = A[i,j]/D[i,i]/D[j,j]
		}
	}
	return(L)
}

# Extract matrix X w.r.t k-th largest eigenvalue of L and run k-mean
spectral.kmeans = function(L, k = 2){
	prin = prcomp(L)
	# extract matrix X for k-mean algorithm
	X = prin$rotation[,1:k]
	return(list(EigenVecMat = X, mem.opt = kmeans(X, k)$cluster))
}

# Print outcomes
printOut = function(A, X, mem, mem.opt){
	cat("Original membership assignment at mem = ", mem, '\n')
	cat("Optical membership assignment at mem = ", mem.opt, '\n')
	#plot(X, col = mem.opt, main = "Clusters of X")
	return(list(Adjacency = A, EigenVecMat = X))
}

