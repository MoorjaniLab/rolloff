######### Fit exponential distribution to ROLLOFF output ############
# Author: Priya Moorjani
# rexpfit.r is a program that can be used to fit an exponential distribution to the output of ROLLOFF. This program uses the nls function in R to determine the nonlinear (weighted) least-squares estimates of the parameters of a nonlinear model. It also uses the DEoptim package to pick the initial values of all the parameters (Katharine Mullen, David Ardia, David Gil, Donald Windover, James Cline (2011). DEoptim: An R Package for Global Optimization by Differential Evolution. Journal of Statistical Software, 40(6), 1-26. URL http://www.jstatsoft.org/v40/i06/.).

#More information on these functions and packages can be found at-
#http://stat.ethz.ch/R-manual/R-patched/library/stats/html/nls.html
#http://cran.r-project.org/web/packages/DEoptim/index.html

#####################################################################

expfit <- function(input, 			# output from rolloff
					output,  		# output of expfit			
                	col=4,   		 # data column containing weighted correlation
         	       	lval=0.5,		 # lower value of dist to use
            	    hval=100,		 # higher value of dist to use                	
                	affine=TRUE	 	# affine = TRUE/FALSE
)
{

nexp <- 1  		 # number of exponentials = 1
# Read input file
data <- read.table(input, header = F)

# set dist and wcorr
dist <- data[,1]
wcorr <- data[,col]
ndist <- length(dist)  ## number of rows in dataset

# check x lower value and y lower value
data.sub <- data
if ((lval > dist[1]) || (hval < dist[ndist])) {
	data.sub <- subset(data, ((dist <= hval) & (dist >= lval)))
} 
dist <- data.sub[,1]		# updated x values
wcorr <- data.sub[,col]		# updated y values

# fit exponential distribution to data
#install.packages('DEoptim')
library('DEoptim')

if(affine){
	# Pick optimized parameters for model 1: A*exp(m1 * -x) + C
	fm1 <- function(x) x[1] + x[2]*exp(-x[3] * dist/100)
	fm2 <- function(x) sum((wcorr-fm1(x))^2)
	fm3 <- DEoptim(fm2, lower=rep(0,3), upper=c(10e7, 10, 300), control=list(trace=FALSE))
	par1 <- fm3$optim$bestmem

	# parameters for y ~ Ae-mt + C
	names(par1) <- c("C", "A", "m")

	# fit the model of Single exponential
	fit1 <- nls(wcorr ~ (C + A*exp(-m * dist/100)), start=par1, control=list(maxiter=1000, warnOnly=TRUE))
	
	# estimated parameters
	aff <- coef(fit1)[1]
	xco <- coef(fit1)[2]
    xexp <- coef(fit1)[3]	# rate of decay of exponential
	} else { 
	# Pick optimized parameters for model 1: A*exp(m1 * -x) + C
	fm1 <- function(x) x[1]*exp(-x[2] * dist/100)
	fm2 <- function(x) sum((wcorr-fm1(x))^2)
	fm3 <- DEoptim(fm2, lower=rep(0,2), upper=c(10, 300), control=list(trace=FALSE))
	par1 <- fm3$optim$bestmem

	# parameters for y ~ Ae-mt
	names(par1) <- c("A", "m")

	# fit the model of Single exponential
	fit1 <- nls(wcorr ~ (A*exp(-m * dist/100)), start=par1, control=list(maxiter=1000, warnOnly=TRUE))
	
	#estimated parameters
	xco <- coef(fit1)[1]
	xexp <- coef(fit1)[2]  	# rate of decay of exponential
	}
	
#print output
outlog <- paste(output, "log", sep='.')
cat("Summary of fit:\n", file=outlog, append = FALSE)
capture.output(summary(fit1), file = outlog, append = TRUE)

# Convert distance in Morgans
xexp <- format(xexp, digits=4, width=0, justify="right")
cat(paste("Estimated number of generations since mixture: ", xexp, "\n"), file=outlog, append=TRUE)

# print fit
fitout <- paste("fit_", output, sep='')
fitdata <- data.frame(dist,fitted(fit1))
fitdata <- format(fitdata, digits=5, justify="right")
write.table(fitdata, fitout, sep="\t", col.names = c("genetic_distance(cM)", "Fitted_values"), row.names=F)

# plot results
outpdf <- paste(output, "pdf", sep='.')
pdf(outpdf)
plot(dist, predict(fit1), type="l", col= "red", xlab = "Genetic Distance (cM)", ylab = "Weighted Correlation", main = "ROLLOFF results")
points(dist,wcorr, col= "blue")
legend("topright",cex=0.75, pch=c(1,NA), lty=c(NA,1), col=c("blue", "red"), legend=c("Admixture LD", "fit"))
dev.off()


} # end of expfit
