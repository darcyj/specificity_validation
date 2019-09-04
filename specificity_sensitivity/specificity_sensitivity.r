
## load package
	library(specificity)

## universal parameters for all sims (g=global)
	g_nsim <- 300
	g_nspec <- 500
	g_nsamp <- 500
	g_nobs <- 1400
	g_ncores <- 3
	g_ups <- seq(from=0, to=1, length.out=g_nspec)

## generate fake data
	envvar <- seq(from=0, to=1000, length.out=g_nsamp)
	# artificial tree
	set.seed(12345)
	hostre <- rtree(200)
	hosvar <- rep(hostre$tip.label, length=g_nsamp)
	#hosidl <- sample(hostre$tip.label, 1)
	hosidl <- "t85"
	# grid centered on 0, just for ease of simulation
	sgrid <- randomgrid(n_samp=g_nsamp)


## Vector input (elevation)

	## Sensitivity of SES to specificity, with specificity done via expanding sd of p(s)
		# simulate matrix
		vec_sdevs <- seq(from=100, to=450, length.out=g_nspec)
		vec_sim_sd <- env_spec_sim( sdev=vec_sdevs, ideal=600, env=envvar, n_obs=g_nobs, 
			n_cores=g_ncores, up=0)
		# calculate specificity
		vec_spec_sd <-  phy_or_env_spec(vec_sim_sd$matrix, env=envvar, 
			n_cores=g_ncores, n_sim=g_nsim)

	## Sensitivity of SES to specificity, with specificity done via expanding up 
		# simulate matrix
		vec_sim_up <- env_spec_sim( sdev=min(vec_sdevs), ideal=600, env=envvar, n_obs=g_nobs,
			n_cores=g_ncores, up=g_ups)
		# calculate specificity
		vec_spec_up <-  phy_or_env_spec(vec_sim_up$matrix, env=envvar, n_cores=g_ncores,
			n_sim=g_nsim)

	## save (just in case)
		save(list=ls(), file="spec_sens.rdata")
	## plots (if running line-by-line)
		par(mfrow=c(2,1))
		plot(vec_spec_sd$SES ~ vec_sim_sd$params$sd)
		plot(vec_spec_up$SES ~ vec_sim_up$params$up)


## Geographic input (sgrid)

	## Sensitivity of SES to specificity, with specificity done via expanding sd of p(s)
		# simulate matrix
		geo_sdevs <- seq(from=20, to=100, length.out=g_nspec)
		geo_sim_sd <- geo_spec_sim( sdev=geo_sdevs, n_obs=g_nobs, grid=sgrid,
			ideal_x=0, ideal_y=0, n_ideal=1, up=0, n_cores=g_ncores)
		# calculate specificity
		geo_spec_sd <- phy_or_env_spec(geo_sim_sd$matrix, env=dist(sgrid), 
			n_cores=g_ncores, n_sim=g_nsim)

	## Sensitivity of SES to specificity, with specificity done via expanding up 
		# (uniform proportion) of p(s)
		# simulate matrix
		geo_sim_up <- geo_spec_sim( min(geo_sdevs), n_obs=g_nobs, grid=sgrid,
			ideal_x=0, ideal_y=0, n_ideal=1, up=g_ups, n_cores=g_ncores)
		# calculate specificity
		geo_spec_up <- phy_or_env_spec(geo_sim_up$matrix, env=dist(sgrid), 
			n_cores=g_ncores, n_sim=g_nsim)

	## save (just in case)
		save(list=ls(), file="spec_sens.rdata")
	## plots (if running line-by-line)
		par(mfrow=c(2,1))
		plot(geo_spec_sd$SES ~ geo_sim_sd$params$sd)
		plot(geo_spec_up$SES ~ geo_sim_up$params$up)


## Phylogenetic input

	## Sensitivity of SES to specificity, with specificity done via expanding sd of p(s)
		phy_sdevs <- seq(from=2, to=5, length.out=500)
		# simulate matrix
		phy_sim_sd <- phy_spec_sim( sdev=phy_sdevs, ideal=hosidl, n_ideal=1, hosts=hosvar, 
			hosts_phylo=hostre,	n_obs=g_nobs, up=0, n_cores=g_ncores)
		# calculate specificity
		phy_spec_sd <- phy_or_env_spec(phy_sim_sd$matrix, hosts=hosvar, hosts_phylo=hostre,
			n_cores=g_ncores, n_sim=g_nsim)

	## Sensitivity of SES to specificity, with specificity done via expanding up 
		# (uniform proportion) of p(s)
		# simulate matrix
		phy_sim_up <- phy_spec_sim( sdev=min(phy_sdevs), ideal=hosidl, n_ideal=1, hosts=hosvar,
			hosts_phylo=hostre, n_obs=g_nobs, up=g_ups, n_cores=g_ncores)
		# calculate specificity
		phy_spec_up <- phy_or_env_spec(phy_sim_up$matrix, hosts=hosvar, hosts_phylo=hostre,
			n_cores=g_ncores, n_sim=g_nsim)

	## save (just in case)
		save(list=ls(), file="spec_sens.rdata")
	## plots (if running line-by-line)
		par(mfrow=c(2,1))
		plot(phy_spec_sd$SES ~ phy_sim_sd$params$sd)
		plot(phy_spec_up$SES ~ phy_sim_up$params$up)


## make plot
	pdf("spec_sens.pdf", useDingbats = F)
	par(mfrow=c(3,2))
	# first row - vector
	plot(vec_spec_sd$SES ~ vec_sim_sd$params$sd, pch=20, ylab="SES", xlab="SD of P(s|elev), UP=0")
	plot(vec_spec_up$SES ~ vec_sim_up$params$up, pch=20, ylab="SES", xlab="UP of P(s|elev), SD=100")
	# second row - geography
	plot(geo_spec_sd$SES ~ geo_sim_sd$params$sd, pch=20, ylab="SES", xlab="SD of P(s|geodist), UP=0")
	plot(geo_spec_up$SES ~ geo_sim_up$params$up, pch=20, ylab="SES", xlab="UP of P(s|geodist), SD=20")
	# third row - phylogeny
	plot(phy_spec_sd$SES ~ phy_sim_sd$params$sd, pch=20, ylab="SES", xlab="SD of P(s|phylodist), UP=0")
	plot(phy_spec_up$SES ~ phy_sim_up$params$up, pch=20, ylab="SES", xlab="UP of P(s|phylodist), SD=2")
	# writeout pdf
	dev.off()


	