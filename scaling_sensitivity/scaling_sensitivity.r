
## load package
	library(specificity)

## generate fake data
	envvar <- seq(from=0, to=1000, length.out=500)
	sim <- env_spec_sim( sdev=200, ideal=600, env=envvar, n_obs=10000, 
			n_cores=3, up=0, oceanp=0.25)
plot(sim$matrix[,1] ~ envvar)

## scale envvar by 1 to 100 and calculate SES
	scales <- 1:100
	pval <- ses <- rep(0, length(scales))
	for(i in 1:length(scales)){
		envvar_s <- envvar * scales[i]
		res_s <-  phy_or_env_spec(sim$matrix, env=envvar_s, n_cores=3, n_sim=1000)
		ses[i] <- res_s$SES[1]
		pval[i] <- res_s$Pval[1]
	}

## found tiny tiny difference due to rounding errors, so round everything to 4 decimal places
	pval <- round(pval, 4)
	ses <- round(ses, 4)

## make a plot
	pdf("scaling_sens.pdf")
	par(mfrow=c(2,1))
	plot(ses ~ scales, pch=20, xlab="Scaling multiplier", ylab="specificity (SES)")
	plot(pval ~ scales, pch=20, xlab="Scaling multiplier", ylab="P-value")
	dev.off()

	