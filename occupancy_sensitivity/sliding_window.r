# x - x variable
# y - y variable
# f - function to apply to each window
# w - width of each window in units of x
# n - number of windows
sliding_window <- function(x, y, f, w, n=NULL){
	# find default n (max non-overlap)
	if(is.null(n)){
		n <- floor((max(x) - min(x)) / w)
	}
	# find midpoints of far left and far right windows
	w2 <- w/2
	l_mp <- min(x) + w2
	r_mp <- max(x) - w2
	# get all midpoints
	midps <- seq(from=l_mp, to=r_mp, length=n)
	# midpoint index wrapper function
	fi <- function(i){
		mp <- midps[i]; lb=mp-w2; ub=mp+w2
		if(i != n){
			yw <- y[x > lb & x < ub]
		}else{
			yw <- y[x > lb & x <= ub]
		}
		return(data.frame(
			f_y=f(yw),
			x_mean=mean(x[x > lb & x <= ub]),
			x_lo=lb,
			x_hi=ub
		))
	}
	output <- do.call("rbind", lapply(X=1:n, FUN=fi))
	colnames(output)[colnames(output) == "f_y"] <- deparse(substitute(f))
	return(output)
}


# up - uniform proportion
# sdev - standard deviation
# ideal - mean of distribution
# env - environmental vector
# ... - arguments for plot(), xlab, ylab, etc
plot_combined_dist_vec <- function(up, sdev, ideal, env, ...){
	fakeenv <- seq(from=min(env), to=max(env), length=150)
	unif00 <- dnorm(x=fakeenv, mean=ideal, sd=sdev)
	pmax00 <- max(unif00)
	# downscale unif part
	unif_scaled <- unif00 * (1-up)
	unif_scaled <- unif_scaled + sum(unif00 - unif_scaled)/length(unif_scaled)
	# add bottom corners for plotting polygon
	fakeenv <- c(fakeenv[1], fakeenv, fakeenv[length(fakeenv)])
	unif_scaled <- c(0, unif_scaled, 0)
	plot(NULL, xlim=range(fakeenv), ylim=c(0, pmax00), xaxs="i", yaxs="i", ...)
	polygon(x=fakeenv, y=unif_scaled, col="lightblue")
	segments(x0=env, y0=0, x1=env, y1=0.03 * pmax00, col="blue")

}



# up - uniform proportion
# sdev - standard deviation
# ideal - ideal host
# tree - phylogeny of hosts
# hosts - hosts in data set
# ... - arguments for plot(), xlab, ylab, etc
plot_combined_dist_phy <- function(up, sdev, ideal, tree, hosts, ...){
	ns <- make_nested_set(tree)
	dists2ideal <- sapply(X=unique(hosts), FUN=bl_distance_ns, tipb=ideal, tree=tree, ns=ns)
	fakedis <- seq(from=0, to=max(dists2ideal), length=150)
	unif00 <- dnorm(x=fakedis, mean=0, sd=sdev)
	pmax00 <- max(unif00)
	# downscale unif part
	unif_scaled <- unif00 * (1-up)
	unif_scaled <- unif_scaled + sum(unif00 - unif_scaled)/length(unif_scaled)
	# add bottom corners for plotting polygon
	fakedis <- c(fakedis[1], fakedis, fakedis[length(fakedis)])
	unif_scaled <- c(0, unif_scaled, 0)
	plot(NULL, xlim=range(dists2ideal), ylim=c(0, pmax00), xaxs="i", yaxs="i", ...)
	polygon(x=fakedis, y=unif_scaled, col="lightblue")
	segments(x0=dists2ideal, y0=0, x1=dists2ideal, y1=0.03 * pmax00, col="blue")
}