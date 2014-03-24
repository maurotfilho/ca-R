
# Correlation of CAs


###############################################################################
# Entropy of n
###############################################################################
###############################################################################
correlation_first_last <- function (y){

	N_iters <- nrow(y)

    # Correlation between the input and output of the CA
	cat("---------------------------------------------\n")
	cat("Correlation (1st, last):", cor(y[1,],y[N_iters,]), "\n")
	cat("---------------------------------------------\n")
	
}


###############################################################################
###############################################################################
###############################################################################
correlation_first_iter <- function (y){

	N_iters <- nrow(y)
	cor_v <- numeric(N_iters)

	for (i in 1:N_iters)
	{
		cor_v[i] <- abs(cor(y[1,],y[i,]))
	}

	plot(cor_v, type="o", main="Correlação 1st-iter", ylim=c(0,1))

}



###############################################################################
###############################################################################
###############################################################################
correlation_iter_iter <- function(y) {

	N_iters <- nrow(y)
	cor_v <- numeric(N_iters)

	for (i in 2:N_iters)
	{
		cor_v[i-1] <- abs(cor(y[i,],y[i-1,]))
	}

	plot(cor_v, type="o", main="Correlação iter-iter", ylim=c(0,1))

}



###############################################################################
###############################################################################
###############################################################################
correlation_iter <- function(y, iter) {

	N_iters <- nrow(y)
	cor_v <- numeric(N_iters-iter)

	for (i in iter:(N_iters-1))
	{
		cor_v[i-iter+1] <- abs(cor(y[iter,],y[i+1,]))
	}

	plot(cor_v, type="o", main=sprintf("Correlação iter %d",iter), ylim=c(0,1))

}

run_correlation <- function() {

	source("ca.R")
	layout(mat=matrix(c(1,2,3,4), nrow=2, byrow=TRUE))
	
	y1 <- runCA(rule=30, seed=initializeCells(0))

	correlation_first_iter(y1[1:100,])
	correlation_iter_iter(y1[1:100,])

	y2 <- runCA(rule=30, seed=initializeCells(1))

	correlation_first_iter(y2[1:100,])
	correlation_iter_iter(y2[1:100,])



}
