# Matrix operations ----
	is.pd = function(X, spd = F, tol = 0, warn = T) {
		# Check if X matrix is positive definite
		# X:    a symmetric matrix
		# spd:  check if the matrix is *semi*-positive definide?
		# tol:  a tolerance for the autovalues
		# warn: print a warnning if the matrix is not symmetric
		
		X = unname(as.matrix(X))
		
		# Check simmetry
		if (!isSymmetric.matrix(X)) {
			if (warn)
				warning("[is.pd]: The matrix is not symmetric")
			return(FALSE)
		}
		
		# Get autovalues
		autoV = eigen(X, symmetric = T, only.values = T)$values[nrow(X)]  
		
		if (is.numeric(autoV)) { # Check if the autovalues aren't complex
			if (spd) {
				return(autoV[length(autoV)] >= -tol)
			} else {
				return(autoV[length(autoV)] > 0)
			}
		} else {
			return(FALSE)
		}
	}
