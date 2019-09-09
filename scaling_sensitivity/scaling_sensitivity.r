
## load package
	library(specificity)

## generate fake data
	envvar <- seq(from=0, to=1000, length.out=500)
	sim <- env_spec_sim( sdev=200, ideal=600, env=envvar, n_obs=10000, 
			n_cores=3, up=0, oceanp=0.25)

## scale envvar by 1 to 100 and calculate SES
	scales <- (1:25) * 4
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


## add intercept to envvar by 1 to 100 and calculate SES
	ints <- (1:25) * 4
	pval2 <- ses2 <- rep(0, length(ints))
	for(i in 1:length(ints)){
		envvar_s <- envvar * ints[i]
		res_s <-  phy_or_env_spec(sim$matrix, env=envvar_s, n_cores=3, n_sim=1000)
		ses2[i] <- res_s$SES[1]
		pval2[i] <- res_s$Pval[1]
	}

## found tiny tiny difference due to rounding errors, so round everything to 4 decimal places
	pval2 <- round(pval2, 4)
	ses2 <- round(ses2, 4)

## make a plot
	pdf("intercept_sens.pdf")
	par(mfrow=c(2,1))
	plot(ses2 ~ ints, pch=20, xlab="Intercept Adder", ylab="specificity (SES)")
	plot(pval2 ~ ints, pch=20, xlab="Intercept Adder", ylab="P-value")
	dev.off()

## save workspace
	save(list=ls(), file="stuff.rdata")