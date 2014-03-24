source("ca.R")


# Unidimensional

# caso simples: 
#		horizontal 	=> seleciona os dados a cada T iterações
#					se T for igual a N_iters, seleciona apenas a última linha da evolução
#
#		vertical 	=> seleciona os dados da coluna C
					

extract1d <- function(y, type, T=1, C=1, i_begin=1)
{
	n_cells <- ncol(y)
	n_iters <- nrow(y)
	
	if (type == "vertical")
	{
		dados <- y[i_begin:n_iters,C]
	}
	else if (type == "horizontal")
	{
		# caso mais simples que seleciona apenas a ultima linha da evolução
		T <- n_iters
		dados <- y[T,]
	}
	dados

}

runExtract1d <- function(type, output.filename)
{
	write(x="", file=output.filename)

	N_cells <- 1024
	N_iters <- 500000

	y <- runCA(rule=c(30), seed=initializeCells(1), n_iters=N_iters)

	# horizontal ignores C, i_begin
	# vertical uses them
	dados <- extract1d(y, type, C=N_cells/2, i_begin=1)	
	write(paste(dados, collapse=""), file=output.filename, append=TRUE)

}


output.filename <- "data.ca30h"
output.filename <- "data.ca30v"

runExtract1d("vertical", output.filename)

