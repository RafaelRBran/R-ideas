# Time Series related ----
	# TODO: Improve all returns adding test statistics, confidence intervals etc.
	# Tests based on artigo de Christofersen and Peletier (2003)
	markov_uc.test = function(x, VaR, tau = 0.05) {
		# P-value of Kupiec's test
		# x:   a vector representing a time series
		# VaR: tau level value-at-risk
		# tau: the level of VaR
		
		if (tau < 0 | tau > 1)
			stop("[markov_uc.test]: tau should be in the interval [0, 1]")
		
		n = length(x)
		n1 = sum(x < VaR) # number of VaR violations
		n0 = n - n1 # number of non-violations
		
		if (n1 == 0) {
			warning("[markov_uc.test]: there are no violations to the VaR")
			log.LR = n0*log(1 - tau)
		} else if (n0 == 0) {
			warning("[markov_uc.test]: every value violates the VaR")
			log.LR = n1 * log(tau)
		} else {
			log.LR = n1*log(tau) + n0*log(1 - tau) - n1*log(n1/n) - n0*log(n0/n)
		}
		
		stat.LR = -2*log.LR
		valor.p = 1 - pchisq(stat.LR, 1)
		return(valor.p)
	}
	
	markov_cc.test = function(x, VaR, tau = 0.05) {
		# Markovian test of conditional coverage
		# x:   a vector representing a time series
		# VaR: tau level value-at-risk
		# tau: The level of VaR (expected proportion of violations)
		if (tau < 0 | tau > 1)
			stop("[markov_cc.test]: tau should be in the interval [0, 1]")
		
		viol = x < VaR
		viol_ind = which(viol)
		
		n = length(x)
		
		n1 = sum(viol)
		n0 = n - n1
		
		n11 = sum(viol[viol_ind + 1], na.rm = T)
		n01 = sum(!viol[viol_ind - 1])
		
		if (n1 == 0) {
			warning("[markov_cc.test]: n1 = 0")
			log.LR = n*log(1 - tau)
		} else if (n0 == 0) {
			warning("[markov_cc.test]: n1 = 0")
			log.LR = n*log(tau)
		} else if (n11 == 0) {
			warning("[markov_cc.test]: n11 = 0")
			log.LR = n1*log(tau) + n0*log(1 - tau) - 
				(n0 - n1)*log(1 - n1/n0) - n1*log(n1/n0)
		} else if (n01 == 0) {
			warning("[markov_cc.test]: n01 = 0")
			log.LR = n1*log(tau) + n0*log(1 - tau)
		} else {
			pi01_hat = n01/n0
			pi11_hat = n11/n1
			
			log.LR = n1*log(tau) + n0*log(1 - tau) - 
				(n0 - n01)*log(1 - pi01_hat) - n01*log(pi01_hat) - 
				(n1 - n11)*log(1 - pi11_hat) - n11*log(pi11_hat)
		}
		
		stat.LR = -2*log.LR
		valor.p = 1 - pchisq(stat.LR, 2)
		return(valor.p)
	}
	
	markov_ind.test = function(x, VaR) {
		# Markovian test of independence
		# x:   a vector representing a time series
		# VaR: a vector representing the value-at-risk time series
		
		viol = x < VaR
		viol_ind = which(viol)
		
		n = length(x)
		n1 = sum(viol)
		n0 = n - n1
		
		if (n1 == 0 | n0 == 0)
			stop("[markov_ind.test]: there are no violations or only every vale on x are a violation")
		
		n11 = sum(viol[viol_ind + 1], na.rm = T)
		n01 = sum(!viol[viol_ind - 1])
		
		pi01_hat = n01/n0
		pi11_hat = n11/n1
		
		pi1_hat = n1/n
		
		if (n11 == 0) {
			warning("[markov_ind.test]: n11 = 0")
			log.LR = n1*log(pi1_hat) + n0*log(1 - pi1_hat) - 
				(n0 - n01)*log(1 - pi01_hat) - n01*log(pi01_hat) 
		} else if (n01 == 0) {
			warning("[markov_ind.test]: n01 = 0")
			log.LR = n1*log(pi1_hat) + n0*log(1 - pi1_hat)
		} else {
			log.LR = n1*log(pi1_hat) + n0*log(1 - pi1_hat) - 
				(n0 - n01)*log(1 - pi01_hat) - n01*log(pi01_hat) - 
				(n1 - n11)*log(1 - pi11_hat) - n11*log(pi11_hat)
		}
		
		stat.LR = -2*log.LR
		valor.p = 1 - pchisq(stat.LR, 1)
		return(valor.p)
	}
	
# Refs ----
	# Christoffersen and Pelletier (2003) Backtesting Value-at-Risk: A Duration-Based Approach, Working paper
