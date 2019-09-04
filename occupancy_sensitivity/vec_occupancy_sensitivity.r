
## prepare environment for testing
	library(specificity)
	#load data
	data(endophyte)
	envvar <- endophyte$metadata$Elevation
	source("sliding_window.r")

## Sensitivity of SES to occupancy
	# using default "global_unif" denominator setting
	# simulate matrices three ways:
		# 1. no uniform proportion (all normal)
		# 2. up = 0.20
		# 3. up = 0.40
	# each will be analyzed both weighted and unweighted.
	# define factors common to all 3 sims:
	prm <- list( sdev = 200, ideal = 1200, env = envvar, 
		n_obs = round(seq(from=5, to=1250, length=500)),
		n_cores=3, seed = 1234567, n_sim = 300, verbose = T)
	
	# simulate matrices
	sim_up00 <- env_spec_sim( sdev=prm$sdev, ideal=prm$ideal, env=prm$env, 
		n_obs=prm$n_obs, n_cores=prm$n_cores, seed=prm$seed, up=0.00)
	sim_up30 <- env_spec_sim( sdev=prm$sdev, ideal=prm$ideal, env=prm$env, 
		n_obs=prm$n_obs, n_cores=prm$n_cores, seed=prm$seed, up=0.30)
	sim_up60 <- env_spec_sim( sdev=prm$sdev, ideal=prm$ideal, env=prm$env, 
		n_obs=prm$n_obs, n_cores=prm$n_cores, seed=prm$seed, up=0.60)

	# calculate specificity of each
	spec_up00 <- phy_or_env_spec(sim_up00$matrix, env=prm$env, n_cores=prm$n_cores, n_sim=prm$n_sim)
	spec_up30 <- phy_or_env_spec(sim_up30$matrix, env=prm$env, n_cores=prm$n_cores, n_sim=prm$n_sim)
	spec_up60 <- phy_or_env_spec(sim_up60$matrix, env=prm$env, n_cores=prm$n_cores, n_sim=prm$n_sim)
	save(list=ls(), file="vec_occ.rdata")

	# make plot - left to right, top to bottom
	ylimits <- range(c(spec_up00$SES, spec_up30$SES, spec_up60$SES))
	shannon <- function(x){ x <- x[x > 0]; x <- x/sum(x); -1 * sum(x * log(x))}
	plotsrow <- function(up, sim, spec, ylm=ylimits){
		# 1. probability distribution for up
		plot_combined_dist_vec(up, prm$sdev, prm$ideal, envvar, ylab=paste("P(s|env). up =", up), xlab="Elevation")
		# 2. effect of occupancy
		occ <- colSums(sim$matrix > 0)
		lm_occ <- lm(spec$SES ~ occ)
		plot(spec$SES ~ occ, pch=20, xlab="occupancy", ylab="SES", ylim=ylm)
		#abline(lm_occ, col="red", lwd=2)
		# 3. SES sensitivity to n_obs
		plot(spec$SES ~ colSums(sim$matrix), pch=20, xlab="# observations", ylab="SES", ylim=ylm)
		# 4 effect of shannon
		#shan <- apply(X=sim$matrix, MAR=2, FUN=shannon)
		#plot(spec$SES ~ shan, pch=20, xlab="Shannon", ylab="SES", ylim=ylm)
		# 4. variance
		#vdat <- sliding_window(colSums(sim$matrix), lm_occ$residuals, var, 30)
		#plot(var ~ x_mean, data=vdat, type="l", col="red", lwd=2, xlab="# observations", ylab="var(specificity)")
	}
	pdf("vec_occ.pdf", useDingbats = F)
	par(mfrow=c(3,3))
	plotsrow(0, sim_up00, spec_up00)
	plotsrow(0.3, sim_up30, spec_up30)
	plotsrow(0.6, sim_up60, spec_up60)
	dev.off()

# using default "species_sim" denominator setting
	# just to show how it doesn't work
	spec_bad <- phy_or_env_spec(sim_up00$matrix, env=prm$env, n_cores=prm$n_cores,
		n_sim=prm$n_sim, denom_type="species_sim")
	occs <- colSums(sim_up00$matrix > 0)
	shannons <- apply(X=sim_up00$matrix, MAR=2, FUN=shannon)
	pdf("vec_occ_speciessim.pdf", useDingbats = F)
	plot(spec_bad$SES ~ occs, pch=20, xlab="occupancy", ylab="SES (\"species sim\")")
	dev.off()