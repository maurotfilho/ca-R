

library(digest)


get_position <- function(input, begin=1, end=2)
{
	hex <- substring(input,begin,end)
	
	hex <- paste("0x",hex, sep="")
	
	return(as.numeric(hex))
}

bloom <- function(y)
{
	# bits tabela
	bits <- 20

	# cria tabela
	tables <- numeric(2**bits)

	# contador de ciclos
	count <- 0

	N_iters <- nrow(y)

	for (i in 1:N_iters)
	{
		hash <- digest(y[i,], algo="sha1")
		pos <- get_position(hash, begin=1, end=5) + 1
		#cat("position ", pos, "\n")
		#cat("count ", count, "\n")

		if (tables[pos] == 1)
		{
			count <- count + 1
			if (count == 8)
			{
				cat("Ciclo! Iter: ", i,"\n")
				count <- 0
			}
		}
		else
		{
			count <- 0
		}

		tables[pos] <- 1
	}

	#cat(tables)
}


run_bloom <- function()
{

	source('ca.R')

	y1 <- runCA(rule=30, seed=initializeCells(0))
	bloom(y1)

	y2 <- runCA(rule=30, seed=initializeCells(1))
	bloom(y2)

}

