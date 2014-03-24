
source('ca.R')
#library(igraph)

bin2dec <- function(x)
{
	sum(x * 2^(rev(seq_along(x)) - 1))
}

dec2bin <- function(x, k, n)
{

	b <- numeric(n)
	j <- 1
	
	while(x > 0)
	{
		b[j] <- x %% k	# modular division (remainder)
		j <- j + 1
		x <- x %/% k	# integer division (quocient)
	}
	rev(b)
}

N_cells <- 7

main <- function(N_cells, k=2)
{
	n <- k**N_cells
	check <- logical(n)
	nstate.next <- numeric(N_cells)

	adj.mat <- matrix(logical(n**2), ncol=n, nrow=n)

	cycles <- vector()

	# Percorrer todas as possibilidades
	for (i in 0:(n-1))
	{
		#cat('Seed: ', i, '\n')
		if (check[i+1] == FALSE)
		{
			#i_cycle <- 1
		
			check[i+1] <- TRUE

			seed <- dec2bin(i, k, N_cells)
			cycle <- i

			nstate.next <- nextStateGen(seed)
			nstate.prev <- nstate.next

			cycle <- c(cycle, bin2dec(nstate.next) )

			adj.mat[bin2dec(seed)+1,bin2dec(nstate.next)+1] <- 1

			#cat(bin2dec(seed), ' ', bin2dec(nstate.next), ' ')
			while (check[bin2dec(nstate.next) + 1] == FALSE)
			{
				check[bin2dec(nstate.prev) + 1] <- TRUE
				nstate.next <- nextStateGen(nstate.prev)
				adj.mat[bin2dec(nstate.prev)+1,bin2dec(nstate.next)+1] <- 1
				#cat(bin2dec(nstate.next), ' ')
				nstate.prev <- nstate.next
				cycle <- c(cycle, bin2dec(nstate.next) )
			}

			cycle_check <- which(cycle == tail(cycle,1))
			if (length(cycle_check) > 1)
			{
				cycles <- c(cycles, cycle_check[2] - cycle_check[1])
				#cat('\n')
				#cat('ciclo tamanho: ', cycle_check[2] - cycle_check[1], '\n')
			}

			#cat('\n')

		}

	}

	#ratio <- max(cycles)/n
	#cat("total, min, 1st, median, 3rd, max, ratio\n")
	#cat(length(cycles), quantile(cycles), ratio)
	#cat("\n")

	#if (ratio >= 0.5)
	#{cat("Regra boa: ",rules,"\n")}

	g <- graph.adjacency(adj.mat)
	#summary(g)
	#L <- layout.fruchterman.reingold(g)
	#plot(g,layout=L)
	#plot(g)
	tkplot(g)

	cycles
}

go <- function()
{
	N_cells <- 12
	d <- data.frame()

	filename <- sprintf('quartis-n%d-r%d-t%d.csv', N_cells, r, totalistic)

	for(i_rule in 1:(k**(k**(2*r+1)) / 4 ) )
	#for(i_rule in 1:100)
	{
		cat('Regra', i_rule, '\n')
		rules <- i_rule
		cycles <- main(N_cells)
		d <- rbind(d, data.frame(t(c(i_rule, length(cycles), t(quantile(cycles)), max(cycles)/(k**N_cells)))))

	}
	write.csv(d, file=filename)
}