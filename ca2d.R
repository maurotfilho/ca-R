###############################################################################
# Le packages
###############################################################################
require(plotrix)
#require(bitops)
#require(rgl)
#require(mgcv)
#require(scatterplot3d)
#require(Rcmdr)

###############################################################################
#  Some definitions
###############################################################################
random.initializer 	<- 0
N_layouts <- 9  	# number of layouts on the grid
N_cells <- 25 		# number of cells
N_iter 	<- 50 		# number of iterarions
r 		<- 1 		# TODO: neighborhood radius -> JUST WORKS WITH r = 1, for now
k		<- 2 		# numbers of states
dim		<- 2 		# number of dimensions
totalistic <- 2 	# TODO: 0 -> non-totalistic; 
					# 1 -> totalistic-furreca; 
					#		sum (binary): ... 9 8 7 6 5 4 3 2 1 0
					# 2 -> totalistic-plus
					#		center:			1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 
					#		sum (binary): 	 8   7   6   5   4   3   2   1   0
					#
					#		center: 		...	2 1 0 2 1 0 2 1 0 2 1 0 2 1 0 2 1 0
					#		sum (ternary):	...	  5     4     3     2     1     0
rules <- c(666)

neigh.type <- "vonneumann"		# vonneumann 	<- 23123
								# moore			<- 32194
								# diagonal1 	<- 32194
								# diagonal2		<- 
								# x-unnamed		<- 

PRINT_CELLS = 1
EXPORT_CELLS = 1
SUBSEQ_OVERLAP = 0
PRINT_ENTROPY = 0
PRINT_SCATTERPLOT = 0
PRINT_CONTOURPLOT = 0

#
#	vonneumann			moore 				cole
#
#	#	X	#			X	X	X 			#	X	#
#	X	X	X 			X	X	X 			#	X	#
#	#	X	# 			X	X	X 			X	#	X
#
#
#	diagonal1 			diagonal2 			smith
#
#	#	X	X 			X	X	# 			#	X	#
#	X	X	X 			X	X	X 			#	X	X
#	X	X	# 			#	X	X  			#	#	#
#

###############################################################################
# Rule Code
###############################################################################
#
# Generic rule decode from decimal to binary
#
###############################################################################
decodeRule <- function(rule.v.decimal) 
{

	# outer totalistic takes double bits on rule
	# outer totalistic is represented by 2
	# normal totalistic is represented by 1

	if (neigh.type == "vonneumann")
	{
		rule.v.binary <- matrix(numeric(totalistic * length(rules) * k^(2*r*dim + 1)), # zero fill initializer
							nrow=length(rules), 
							ncol=totalistic * k^(2*r*dim + 1))	
	}
	else if (neigh.type == "moore")
	{
		rule.v.binary <- matrix(numeric(totalistic * length(rules) * k^((2*r+1)**dim)), # zero fill initializer
							nrow=length(rules), 
							ncol=totalistic * k^((2*r+1)**dim))	
	}
	else if (neigh.type == "diagonal1" || neigh.type == "diagonal2")
	{
		rule.v.binary <- matrix(numeric(totalistic * length(rules) * k^7), # zero fill initializer
							nrow=length(rules), 
							ncol=totalistic * k^7)
	}
	else 
	{
		rule.v.binary <- matrix(numeric(totalistic * length(rules) * k^(2*r + 1)), # zero fill initializer
							nrow=length(rules), 
							ncol=totalistic * 7)		
	}
	
	for (i in 1:length(rules)) 
	{
		j <- 1
		rule.decimal <- rule.v.decimal[i]
		while(rule.decimal > 0)
		{
			rule.v.binary[i,j] <- rule.decimal %% k	# modular division (remainder)
			j <- j + 1
			rule.decimal <- rule.decimal %/% k	# integer division (quocient)
		}
	}
#	cat("rule binary: ", rule.v.binary, "\n")

	rule.v.binary
}

###############################################################################
# Neighborhood get
###############################################################################
#
# obtain the neighborhood of xy on index
#
###############################################################################
getNeighborhood <- function(xy, index_i, index_j)
{
	# obtain the indexes of north, east, south, west
	north <- index_j - 1
	if (north <= 0){ north <- north + N_cells}
	south <- (index_j + 1) %% N_cells
	if (south == 0){ south <- N_cells }
	west <- index_i - 1
	if (west <= 0){ west <- west + N_cells}
	east <- (index_i + 1) %% N_cells
	if (east == 0){ east <- N_cells }

	#cat("center:",index_i,',',index_j,'\n')
	#cat(north,',',south,',',west,',',east,'\n')

	if (neigh.type == "vonneumann")
	{
		if (r == 1)
		{
			neighborhood <- numeric(4 * r + 1)		# zero fill initializer
			neighborhood[1] <- xy[index_i,index_j]	# center
			neighborhood[2] <- xy[index_i,north]  	# north
			neighborhood[3] <- xy[east,index_j]		# east
			neighborhood[4] <- xy[index_i,south]	# south
			neighborhood[5] <- xy[west,index_j]		# west
		}
	}

	else if (neigh.type == "moore")
	{
		if (r == 1)
		{
			neighborhood <- numeric((2*r+1)**dim)	# zero fill initializer
			neighborhood[1] <- xy[index_i,index_j]	# center
			neighborhood[2] <- xy[index_i,north]  	# north
			neighborhood[3] <- xy[east,index_j]		# east
			neighborhood[4] <- xy[index_i,south]	# south
			neighborhood[5] <- xy[west,index_j]		# west
			neighborhood[6] <- xy[east,north]  		# northeast
			neighborhood[7] <- xy[east,south]		# southeast
			neighborhood[8] <- xy[west,south]		# southwest
			neighborhood[9] <- xy[west,north]		# northwest
		}
	}

	else if (neigh.type == "diagonal1")
	{
		if (r == 1)
		{
			neighborhood <- numeric(7)
			neighborhood[1] <- xy[index_i,index_j]	# center
			neighborhood[2] <- xy[index_i,north]  	# north
			neighborhood[3] <- xy[east,index_j]		# east
			neighborhood[4] <- xy[index_i,south]	# south
			neighborhood[5] <- xy[west,index_j]		# west
			neighborhood[6] <- xy[east,north]  		# northeast
			neighborhood[7] <- xy[west,south]		# southwest
		}
	}

	else if (neigh.type == "diagonal2")
	{
		if (r == 1)
		{
			neighborhood <- numeric(7)
			neighborhood[1] <- xy[index_i,index_j]	# center
			neighborhood[2] <- xy[index_i,north]  	# north
			neighborhood[3] <- xy[east,index_j]		# east
			neighborhood[4] <- xy[index_i,south]	# south
			neighborhood[5] <- xy[west,index_j]		# west
			neighborhood[6] <- xy[east,south]		# southeast
			neighborhood[7] <- xy[west,north]		# northwest
		}
	}
	else
	{
		cat("Neighborhood Type not defined\n")
	}

    neighborhood
}

###############################################################################
# NextStateGen
###############################################################################
#
# state^{t+1} = rule(state^{t})
#
###############################################################################
nextStateGen <- function(x) 
{
	x_next <- matrix(numeric(N_cells*N_cells), nrow=N_cells, ncol=N_cells)

	rule.v <- decodeRule(rules)

	for (index_i in 1:N_cells) 
	{
		for (index_j in 1:N_cells)
		{
			# obtain the neighborhood of x on index i
			neighborhood <- getNeighborhood(x, index_i, index_j)

			# totalistic
			if(totalistic == 1)
			{
				index.neigh <- sum(neighborhood) + 1
			}
			# outer totalistic
			else if(totalistic == 2)
			{
				index.neigh <- k*sum(neighborhood[2:length(neighborhood)]) + neighborhood[1] + 1
			}

			# NOT-IMPLEMENTED YET!
			# non-totalistic
			else
			{
				# obtain the index on decimal of the neighborhood
				index.neigh <- 1
				for (aux in 0:(2*r)) 
				{
					index.neigh <- index.neigh + k^aux * neighborhood[length(neighborhood) - aux]
				}
			}

			# obtain which rule is going to use
			# in case of uniform CA, it will be only one rule
			# in case of hybrid CA, it will be more than one
			index.rules <- ((index_i - 1) %% length(rules)) + 1
			
			# obtain the value of next state by the rule
			x_next[index_i,index_j] <- rule.v[index.rules,index.neigh]
		}
	}
	x_next
}

###############################################################################
# Evolve
###############################################################################
# 
# evolution of CA for N_iter iterations
# y[iter, cell] <- holds every step on each line
#
###############################################################################
evolve <- function(xy)
{
	# y[N_iter,N_cells] <- 0,...,0
	#					   0,...,0
	#					     ...
	#					   0,...,0
	z <- array(numeric(N_cells*N_cells*N_iter),dim=c(N_cells,N_cells,N_iter))

	xy_aux <- xy
	for (i in 1:N_iter) {
		#cat("Iter:", i, "\n")
		z[,,i] <- xy_aux
		xy_next <- nextStateGen(xy_aux)
		xy_aux <- xy_next
		# cat(y[i,],'\n')
	}
	z
}
###############################################################################
# Entropy of subsequence size 1
###############################################################################
# 
# Basic entropy just to see if the ratio (of ones and zeros) is fine
#
###############################################################################
entropy.1 <- function(x, kbase=k)
{ 
  	# obtain p and log(p) for each element
  	# if its probability p is zero, its log will be -infinite
  	# so it must verify that no one of this appear
  	# because if it's appear, will give us an NaN (and it's not good)
	p <- numeric(kbase)
	log_p <- numeric(kbase)

	for(i in 1:kbase)
	{
	  	# probability of each element to appear in 'x'
	 	p[i] <- sum(x == (i-1))/(length(x))
		if (p[i] <= 0)
		{
			log_p[i] <- 0
		}
		else
		{
			log_p[i] <- log(p[i], base=kbase)
		}
	}

	-sum(log_p*p)
}

###############################################################################
# All Combinations of N
###############################################################################
# TODO:
###############################################################################
allcombn <- function(n, k)
{
	cmbs <- matrix(numeric(n*(k**n)),nrow=k**n,ncol=n)
	
	for(i in 0:((k**n)-1))
	{
		divisor <- i

		for(j in (n-1):0)
		{
			cmbs[i+1,n-j] <- divisor %/% (k**j)
			divisor <- divisor %% (k**j)
		}
	}

	cmbs
}

###############################################################################
# All Subsequences
###############################################################################
#
# obtain all subseqs of x with size n
# with or without overlap
# TODO: try not to truncate de subsequences
#
###############################################################################
allsubseq <- function(x, n, overlap=SUBSEQ_OVERLAP)
{
	if(overlap)
	{
		# WITH OVERLAP
		subseq <- matrix(numeric((length(x)-n+1)*n), nrow=(length(x)-n+1), ncol=n)
		for (i in 1:(length(x) - n + 1))
		{
			#cat(i, ":", i+n-1,  ":")
		    #cat(x[i:(i+n-1)], "\n")
		    subseq[i,] <- x[i:(i+n-1)]
		}
	}
	else
	{
		# WITHOUT OVERLAP
		subseq <- matrix(numeric((length(x)%/%n)*n), nrow=(length(x)%/%n), ncol=n)
    	i <- 1
		for (j in seq(1,(nrow(subseq))))
		{
		    # cat(i, ":", i+n-1,  ":")
		    # cat(x[i:(i+n-1)], "\n")
			subseq[j,] <- x[i:(i+n-1)]
      		i <- i + n
		}
	}

	subseq
}

###############################################################################
# Entropy of n
###############################################################################
#
# x <- CA
# n <- size of subseq 
#      the size must be divisor of the CA length or it will be truncated
#
###############################################################################
entropy <- function(x, n=1, kbase=k) 
{
	# generate all subsequences of x using size n
	sbsq <- allsubseq(x, n, SUBSEQ_OVERLAP)

	# generate all combinations that is possible with size n and k base
	cmbs <- allcombn(n, kbase)
  
  	# count all subseqs
	count <- numeric(nrow(cmbs))
  	for (i in 1:nrow(sbsq))
  	{
    	for (j in 1:nrow(cmbs))
    	{
      		count[j] <- count[j] + all(sbsq[i,] == cmbs[j,])
    	}
  	}

	# probability of each subseq to appear in 'x'
	p <- numeric(nrow(cmbs))
	p <- count/(nrow(sbsq))

  	# log of probability p
	log_p <- numeric(nrow(cmbs))

	for(i in 1:nrow(cmbs))
	{
		if (p[i] <= 0)
		{
			log_p[i] <- 0
		}
		else
		{
			log_p[i] <- log(p[i], base=kbase)
		}
	}

	#  cat("count", count, "\n")
	#  cat("prob", p, "\n")

	-sum(log_p*p)/n
}

###############################################################################
# Print cells
###############################################################################
#
#	Print the cells evolution using plotrix package
# 	  ex: mat = matrix(elemento, ncol=N_COLS, byrow=TRUE)
#    	  color2D.matplot(mat, main=title)
#
###############################################################################
printCells <- function(y)
{
	if (dim == 2)
	{
		color2D.matplot(abs(y-1))
	}
		
}

###############################################################################
# Print cells evolution
###############################################################################
#
# 	Print cells evolution using plotrix package
#
###############################################################################
printCellsEvolution <- function(z)
{
	layout(matrix(1:N_layouts,sqrt(N_layouts),sqrt(N_layouts),byrow=TRUE))
	color2D.matplot((abs(z[,,1]-k+1))/k, main="Iter 1")
	
	# goes from 2 to N_iter-1 equally spaced 
	# (dividing in N_lay spaces)
	for(i in floor(seq(2, (N_iter-1), length=(N_layouts-2)))) 
	{
		#iter <- (i - 1) * (N_iter%/%N_layouts) + 1
		iter <- i
		color2D.matplot((abs(z[,,iter]-k+1)/k), main=iter)
		
	}
	color2D.matplot((abs(z[,,N_iter]-k+1)/k), main=N_iter)		
}

###############################################################################
# Export cells evolution
###############################################################################
#
# Export cells evolution using plotrix package to a PNG file
#
###############################################################################
exportCellsEvolution <- function(z)
{

	for(i in 1:N_iter) 
	{
		#iter <- (i - 1) * (N_iter%/%N_layouts) + 1
		iter <- i

		filename <- sprintf("images/ca2dplot/plotca2d-%04d.png", i)
		png(filename)
		color2D.matplot((abs(z[,,iter]-k+1)/k), main=iter)
		dev.off()
		
	}
	system("convert -delay 40 \"images/ca2dplot/*.png\" images/ca2dplot/animated.gif")
	system("rm images/ca2dplot/*.png")
		
}

###############################################################################
# Initialize Cells
###############################################################################
#
# random = 0 -> 1 cell with 1, and the others are 0
# random = 1 -> random bits to initialize the cells
#
###############################################################################
initializeCells <- function(random = 0) 
{
	if (random)
	{
		x <- matrix(round(runif(N_cells^dim)), nrow=N_cells, ncol=N_cells)
	}
	else
	{
		x <- array(0, dim=rep(N_cells,dim)) # array(data, dim= (N_cells, N_cells)) <=> NxN
											# the size of rep(N_cells,dim) indicates the number of dimensions
											# each one with size N_cells
		if (dim == 2)
			x[N_cells%/%2 + 1,N_cells%/%2 + 1] <- 1
	}
	x
}

###############################################################################
# Detect Loop
###############################################################################
#
# Given all states of a Cellular Automata (z)
# Checks if any state repeats more than once
#
###############################################################################
detectLoop <- function(z) 
{
	for (i_cur in 1:(N_iter-1)) 
	{
		#x_cur <- y[i_cur]
		for (i_compared in (i_cur+1):N_iter) 
		{
			if (all(z[,,i_cur] == z[,,i_compared]))
			{
				cat("Wait! I've found a LOOP on:", i_cur, i_compared, "\n")

			}
		}
	}
	
}

###############################################################################
# show entropy
###############################################################################
#
#
###############################################################################
show.entropy <- function(H, n=1, plot.main="", xlab="", ylab="Entropy", info="", type=NULL)
{
	cat("---------------------------------------------\n")

    H.min <- min(H)
    H.mean <- mean(H)
    H.max <- max(H)		

    if (type == "cell")
    {
    	# contour
    	if (PRINT_CONTOURPLOT)
    		contour(H)
    	# scattered plot 3d
    	if (PRINT_SCATTERPLOT)
    	{
			library(scatterplot3d)
			layout(1,1,1)
			s3d.x=rep(1:N_cells, N_cells)
			s3d.y=numeric(N_cells**dim)
			for(i in 1:N_cells)
			{
				s3d.y[(N_cells*i-N_cells+1):(N_cells*i)] <- rep(i,N_cells)
			}
			s3d.z=H

			s3d <- scatterplot3d(s3d.x,s3d.y,s3d.z, highlight.3d=TRUE, type="p", pch=18)
    	}
    }
    else if (type == "state")
    {
		plot(H, xlab=xlab, ylab=ylab, main=plot.main, ylim=c(0,1))
    	abline(h=H.mean)
    }

    # text(0.1*length(H), 0.4, sprintf("Max Entropy: %.4f", H.max), pos=4)
    # text(0.1*length(H), 0.3, sprintf("Mean Entropy: %.4f", H.mean), pos=4)
    # text(0.1*length(H), 0.2, sprintf("Min Entropy: %.4f", H.min), pos=4)

    cat(info)
	cat("\n")

	cat("Min Entropy:", H.min,'\n')
	cat("at:")
	if(type == "state")
	{
		for (i in 1:length(H)){if(H[i]==H.min) cat(i,' ')}
	}
	else if (type == "cell")
	{

		for (i in 1:N_cells)
		{
			for (j in 1:N_cells)
			{	
				if(H[i,j]==H.min) 
					cat('[',i,',',j,'] ')
			}
		}	
	}

	cat("\n")

	cat("Mean Entropy:", H.mean,'\n')

	cat("Max Entropy:", H.max,'\n')
	cat("at:")
	if(type == "state")
	{
		for (i in 1:length(H)){if(H[i]==H.max) cat(i,' ')}
	}
	else if (type == "cell")
	{
		for (i in 1:N_cells)
		{
			for (j in 1:N_cells)
			{	
				if(H[i,j]==H.max) 
					cat('[',i,',',j,'] ')
			}
		}	
	}
	cat("\n")

	cat("---------------------------------------------\n")
}

###############################################################################
# Run Entropy
###############################################################################
#
# Obtain the entropy and other statistical measures for analysis
# about the randomness of the cellular automata (z[x,y,t]) given
#
###############################################################################
runEntropy <- function(z, n=1, type=NULL)
{
    if (length(type)==0)
    {
        cat("Error: Type not defined.\n")
    }
    else
    {
    	if (type == "cell")
    	{

			# Plot: entropy xy each cells during evolution
			H <- matrix(numeric(N_cells*N_cells), nrow=N_cells, ncol=N_cells)
			for (i in 1:N_cells) 
			{
				for (j in 1:N_cells)
				{
					#cat(i," ",j," ",z[i,j,],"\n")
					H[i,j] <- entropy(z[i,j,],n)	
				}
				
			}

			if (PRINT_ENTROPY)
            {
                plot.main=sprintf("Entropy of subseq size %d",n)
                info.main=sprintf("Entropy of subseq size %d on Cells",n)
                axis.x.lab <- "cell"
                show.entropy(H=H,n=n, plot.main=plot.main, xlab=axis.x.lab, ylab="Entropy",
                    info=info.main, type)
            }
		}
		# perspective 3d
		#persp(z=z)

		#identify(x, y)
		#coords <- locator(type="l") # add lines
		#coords # display list

		else if (type == "state")
		{
			# Plot: entropy x each cells during evolution
			H <- numeric(N_iter)
			for (i in 1:N_iter) 
			{
				H[i] <- entropy(c(z[,,i])) # converting matrix to vector to obtain entropy of xy
			}			
            if (PRINT_ENTROPY)
            {
                plot.main=sprintf("Entropy of subseq size %d",n)
                info.main=sprintf("Entropy of subseq size %d on States",n)
                axis.x.lab <- "iters"
                show.entropy(H=H,n=n, plot.main=plot.main, xlab=axis.x.lab, ylab="Entropy",
                    info=info.main, type)
            }
        }
    }
}


###############################################################################
# Run Stats
###############################################################################
#
# Obtain the entropy and other statistical measures for analysis
# about the randomness of the cellular automata (y[t,x]) given
#
###############################################################################
runStats <- function(z)
{
	layout(1,1,1)
    # Entropy of each state (line) during a CA evolution
    runEntropy(z, n=1, type="state")

    # Correlation between the input and output of the CA
	cat("---------------------------------------------\n")
	cat("Correlation (1st, last):", cor(c(z[,,1]),c(z[,,N_iter])), "\n")
	cat("---------------------------------------------\n")

	# Entropy of each cell (column) during a CA evolution
    runEntropy(z, n=1, type="cell")
    runEntropy(z, n=2, type="cell")
    runEntropy(z, n=3, type="cell")
    runEntropy(z, n=4, type="cell")

	# layout(1,1,1)
	# H <- matrix(numeric(N_iters*N_cells), nrow=N_iters, ncol=N_cells)
	# for (i in 1:N_iters) 
	# {
	# 	for (j in 1:N_cells) 
	# 	{
	# 		H[i,j] <- entropy(y[1:i,j])	
	# 	}
		
	# 	#cat("H[",i,"] = ", H[i],"\n")
	# }
	# exportEntropyPlots(H)

    #autoCorrelation(y, type="state")

}


###############################################################################
# Run CA
###############################################################################
#
# initialize
# evolve
# ... and print pretty-good-CA in two dimensions for us :)
#
###############################################################################
runCA <- function(rule=rules, seed=NULL, n_iters=N_iters) 
{
	# initialize
	cat("Initializing...\n")

    if (length(seed)==0)
        xy <- initializeCells(random.initializer)
    else
        xy <- seed

	# evolve
	cat("Evolving...\n")
	z <- evolve(xy)

	# show some stuff
	cat("Printing...\n")
	if (PRINT_CELLS)
		printCellsEvolution(z)
	if (EXPORT_CELLS)
		exportCellsEvolution(z)

	z
}

###############################################################################
#   Main function
###############################################################################
main <- function() 
{
	z <- runCA()

	cat("Some stats...\n")

	runStats(z)

	#cat("Detecting loop...\n")	
	#detectLoop(z)

	cat("Done!\n")
}

###############################################################################
###############################################################################

main()

# Making animated GIF of .PNG
# convert   -delay 20   -loop 0   files*.png   animated.gif