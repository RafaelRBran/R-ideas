# univariate case
	# Gumpbell distribution
		dgumpbell = function(x) {
			# probability density of gumpbell in x
			return(exp(-exp(-x) - x))
		}
		
		pgumpbell = function(q) {
			# probability function of gumpbell in q
			return(exp(-exp(-q)))
		}
		
		qgumpbell = function(p) {
			# p-esim quantile of gumpbell
			return(-log(-log(p)))
		}
		
		rgumbell = function(n) {
			# gumbell generator by inversion method
			u = runif(n, 0, 1)
			return(qgumpbell(u))
		}
		
	# Weibull distribution
		dweibull = function(x, c, alpha = 1, xi = 0, log = FALSE) {
			# probability density function of weibull in x, according Johnson, Kotz and Balakrishnan in Continuous Univariate Distributions
			if (c <= 0 | alpha <= 0 | x <= xi)
				stop("[dweibull]: inconsistent parameters")
			
			pdf = (c/alpha)*((x - xi)/alpha)^(c - 1) * exp(-((x - xi)/alpha)^c)
			
			if (log) 
				pdf = log(pdf)
			
			return(pdf)
		}
		
		dweibull2 = function(D, b, a = 1, log = FALSE) {
			# probability density function of weibull in D, according rugarch package
			if (b <= 0 | a <= 0 | D <= 0)
				stop("[dweibull2]: inconsistent parameters")
			
			pdf = b*log(a) + log(b) + (b - 1)*log(D) - (a * D)^b
			
			if (log) 
				pdf = log(pdf)
			
			return(pdf)
		}
		
		pweibull = function(q, c, alpha = 1, xi = 0, survival = FALSE) {
			# probability function of weibull in x, according Johnson, Kotz and Balakrishnan in Continuous Univariate Distributions
			if (c <= 0 | alpha <= 0 | q <= xi)
				stop("[pweibull]: inconsistent parameters")
			
			cdf = 1 - exp(-((q - xi)/alpha)^c)
			
			if (survival)
				cdf = 1 - cdf
			
			return(cdf)
		}
		
		pweibull2 = function(q, b, a = 1, survival = FALSE) {
			# probability function of weibull in D, according rugarch package
			if (b <= 0 | a <= 0 | q <= 0)
				stop("[pweibull2]: inconsistent parameters")
			
			cdf = 1 - exp(-(a * q)^b)
			
			if (survival) 
				cdf = 1 - cdf
			
			return(cdf)
		}
		
		qweibull = function(p, c, alpha = 1, xi = 0) {
			# p-esim quantile of weibull (JKB)
			if (c <= 0 | alpha <= 0 | p <= xi)
				stop("[qweibull]: inconsistent parameters")
			
			return(alpha * (-log(1 - p)) ^ (1/c) + c)
		}
		
		rweibull = function(n, c, alpha = 1, xi = 0) {
			# weibull generator by inversion method (JKB)
			if (c <= 0 | alpha <= 0 | p <= xi)
				stop("[rweibull]: inconsistent parameters")
			
			u = runif(n, 0, 1)
			return(qweibull(u, c, alpha, xi))
		}

# multivariate case
	# Multivariate normal based on Falk (1999)
		rmnorm = function(n, mean = 0, S = 1) {
			# Generate a multivariate normal
			# n:    number of observations
			# mean: vector of expectations
				# if mean is an sigle value, the mean will be repeated
			# S:    Variance and covariance matrix
				# if S is a vector, the diag(S) will be used
			
			p = max(NROW(S), length(mean))
			
			# S options (if S isn't pd, chol(S) will return an error)
			if (is.vector(S)) {
				S = diag(S, ncol = p)
			} else {
				S = unname(as.matrix(S)) # check if is a matrix
			}
			
			# Arrumando a média
			if (length(mean) == 1) {
				mean = rep(mean, p)
			} else if (length(mean) != p) {
				stop("[rmnorm]: vetor mean ou matriz S de dimensões erradas")
			}

			return(do.call(rbind, apply(simplify = F,
				matrix(rnorm(n*p, mean = 0, sd = 1), ncol = p) %*% chol(S), # var-cov
				1, '+', mean # add mean
			)))
		}
		
	# Multivariate normal based on Falk (1999)
		rmunif = function(n, min = 0, max = 1, S = ((max - min)^2)/12, exato = F, tol = 1e-3) {
			# Generate a multivarate uniform
			# n: number of observations
			# min: vector of minimus
				# if is an sigle value, it will be repeated
			# max: vector of maximus
				# if is an sigle value, it will be repeated
			# S: var-cov matrix.
				# if S is a vector, the diag(S) will be used
			# exact: Code should try corect the S to the exact matrix?
			# tol: a tolerance
			
			p = max(NROW(S), length(min), length(max))
			
			if (is.vector(S)) {
				S = diag(S, nrow = p)
			} else {
				S = unname(as.matrix(S))
			}
			
			# Check dimensions of min and max
			if (p %% length(min) != 0 ||  p %% length(max) != 0) {
				stop("[rmunif]: min or max dimension are not multiple of the dimension")
			}
			
			# Getting amplitude
			amplitude = max - min
			if (any(amplitude <= 0))
				stop("[rmunif]: max must be bigger than min")
			amplitude_diag     = diag(amplitude, nrow = p)
			amplitude_inv_diag = diag(1/amplitude, nrow = p)
			
			# Getting S of multivariate U(0, 1)
			S01 = amplitude_inv_diag %*% S %*% amplitude_inv_diag 
			
			# Testing if the variance are equal to the teorical one
			if (any(abs(diag(S01) - 1/12) > tol))
				warning("[rmunif]: Na distribuição uniforme, a variância vale ((max - min)^2)/12")
			
			# Exact matrix
			if (exact) {
				S01 = 2*sin((pi/6)*S01)
				if (!is.pd(S, warn = T, spd = F, tol = tol)) # in this case, tol is the same of diag(S01) - 1/12) == 0. Don't need to 
					stop("[rmunif]: S01 corected is not PD")
			}
			
			# Generating multivarate U(0, 1)
			normal = rmnorm(n, mean = 0, S = S01)
			unif01 = pnorm(normal, mean = 0, sd = sqrt(1/12))
			
			# Back to multivariate uniform
			return(do.call(rbind, apply(simplify = F,
				unif01 %*% amplitude_diag,
				1, '+', min
			))) # var-cov is aredy corected
		}
