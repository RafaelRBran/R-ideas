# plot related codes ----
	# Guides ----
		plot_guide = function() {
			# A quick guide for colors and types of points and lines
			
			mfrow = par('mfrow')
			on.exit(par(mfrow = mfrow))
			par(mfrow = c(1, 2))
			
			# Point
			x = rep(1:5, 5)
			y = rep(5:1, each = 5)
			plot(y ~ x, pch = 1:25, xlim = c(.7, 5.6), ylim = c(.7, 5.2), main = "pch")
			text(x + .3, y, 1:25)
			
			# lines
			x = rep(c(4, 5), 4)
			y = rep(4:1, each = 2)
			plot(y ~ x, pch = 19, col = 1:25, xlim = c(0.5, 5.5), ylim = c(0.5, 4.5), main = "lty, col")
			text(x + .3, y, 1:8)
			points(c(1, 3), c(3.5, 3.5), lty = 1, type = "l"); text(3.3, 3.5, "1")
			points(c(1, 3), c(2.5, 2.5), lty = 2, type = "l"); text(3.3, 2.5, "2")
			points(c(1, 3), c(1.5, 1.5), lty = 3, type = "l"); text(3.3, 1.5, "3")
		}
		
	# Quantilic curve ----
		quantilic.plot = function(x, ..., from = 0, to = 1) {
			# plot the quantilic curve
			# x:    a numeric vector of observations
			# ...:  additional parameters to plot()
			# from: the first quantile level to be plotted
			# to:   the last quantile level to be plotted
			tau    = seq(from = from, to = to, length.out = min(1000, length(x[x >= quantile(x, from) & x <= quantile(x, to)])))
			quants = quantile(x, tau, na.rm = T)
			plot(quants ~ tau, type = "l", ...)
			
			invisible(quants)
			
			# Bloat stuff
			#if (show.boxp.stats || show.mean) {
			#	half.lin = (to - from)*.05
			#	acfun    = ecdf(x)
			#}
			#if (show.boxp.stats) {
			#	boxp     = boxplot.stats(x, do.conf = FALSE)
			#	for (stat in boxp$stats) {
			#		x.axis = rep(acfun(stat), 2) + c(-half.lin, half.lin); y.axis = rep(stat, 2)
			#		points(y.axis ~ x.axis, type = "l", lty = 2, col = "red")
			#		if (show.text)
			#			text(x = ifelse(stat < boxp$stats[4], x.axis[2] + half.lin, x.axis[1] - half.lin), y = stat, round(stat, 2))
			#	}
			#		points(x = rep(from, length(boxp$out)), boxp$out, pch = "+")
			#}
			#if (show.mean) {
			#	stat   = mean(x, na.rm = T)
			#	dist.y = (max(quants) - min(quants))*.05
			#	points(x = rep(acfun(stat), 2) + half.lin*c(-2, 2), y = rep(stat, 2), type = "l", col = "blue")
			#	if (show.text)
			#		text(x = acfun(stat) - half.lin, y = stat + dist.y, round(stat, 2), col = "blue")
			#}
		}
			
		bands = function (x, side = 1, boxplot_division = 0, ...) {
			# Inspired on Tukey (1977), Chap. 2
			# TODO: check if an active plot window exists
			# TODO: Add confidence intervals based on Variations of Box Plots by McGill, Tukey and Larsen (1978)
			# Having a quantilic plot, draws lines on y-axis corresponding to observations
			# x:                a numeric vector of observations
			# side:             the axis where the bands are drawn
			# boxplot_division: a single number or logical.
				# If F or 0, all points are represented in '_'.
				# If T or 1, central values are represented as "-", non central values are represented as '+' and extreme values ar represented as '*'.
				# Else, very extreme values are represented as '.'.
				
			# Generating pch
			if (!boxplot_division)
				pch = ifelse(side %% 2 == 0, "_", "|")
			else if (boxplot_division) {
				quants = quantile(x, c(.25, .75))
				step = 1.5*(quants[2] - quants[1])
				pch =
					ifelse(!(boxplot_division %in% c(0, 1)) & (x < quants[1] - 2*step | x > quants[2] + 2*step), 1, # far out
					ifelse(                                    x < quants[1] - step   | x > quants[2] + step,    8, # outside
					ifelse(                                    x < quants[1]          | x > quants[2],           3, # inside outer fences
					                                                                                            20  # inside inner fences
					)))
			}
			
			usr = par('usr')
			dist_x = usr[2] - usr[1]
			dist_y = usr[4] - usr[3]
			
			# ploting points
			if      (side == 1) points(x = x, y = rep(usr[3] + (dist_y*.05), length(x)), pch = pch, ...)
			else if (side == 2) points(y = x, x = rep(usr[1] + (dist_x*.05), length(x)), pch = pch, ...)
			else if (side == 3) points(x = x, y = rep(usr[4] - (dist_y*.05), length(x)), pch = pch, ...)
			else if (side == 4) points(y = x, x = rep(usr[2] - (dist_x*.05), length(x)), pch = pch, ...)
			else stop('[bands]: side bust be 1, 2, 3 or 4')
		}

# Refs ----
	# John W. Tukey (1977) Descriptive data analysis, Addison-Wesley Publishing Company
