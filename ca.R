###############################################################################
# Le packages
###############################################################################
#require(plotrix)
#require(bitops)
source('entropyCA.R')

###############################################################################
#  Some definitions
###############################################################################
N_cells <- 128		# number of cells
N_iters	<- 1280	    # number of iterarions
r    	<- 1 		# neighborhood radius
k		<- 2 		# numbers of states
dim		<- 1 		# number of dimensions
totalistic <- 0 	# 0 -> non-totalistic; 1 -> totalistic
hybrid  <- "CELL" 	# CELL -> each cell evolves by one different rule
					# ITER -> each iteration the rule will be different

wolframrule <- c(30)
garule <- c(1436194405, 1436965290, 1704302169, 1721325161, 1705400746)

rules 	<- c(45)


#c(51,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,
#    30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,
#    30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30)
		# wolfram rule: 30
		# rule: 45
		# fake-entropy rule: 51; 127
        # fake r2 rule: 252645135
		# nandi rules: 51, 153, 195
		# nandi rules: 102, 150, 90, 150
        #  90, 105, 150, 165
		# tomassini and perrenoud rules: 90, 105, 150, 165
		# szaban old rules: 86, 90, 101, 105, 150, 154, 165, 1436194405
		# szaban new rules: 1436194405, 1436965290, 1704302169, 1721325161, 1705400746

random.initializer  = 1
SUBSEQ_OVERLAP = 0
PLOT_AUTOCORRELATION = 0
PRINT_ENTROPY = 0
PRINT_CA = 0

BOUNDARY = 2    # BOUNDARY CONDITIONS:      0 -> null boundary
                #                           1 -> 1 boundary
                #                           2 -> cyclic boundary


###############################################################################
# Rule Code
###############################################################################
#
# Generic rule decode from decimal to binary
#
###############################################################################
decodeRule <- function(rule.v.decimal, kbase=k) 
{
	# zero fill initializer
    if(totalistic)
    {
       rule.v.binary <- matrix(numeric(length(rule.v.decimal) * (1 + (kbase-1)*(2*r + 1))), 
                            nrow=length(rule.v.decimal), 
                            ncol=(1+(kbase-1)*(2*r + 1))
                            ) 
    }
    else
    {
	   rule.v.binary <- matrix(numeric(length(rule.v.decimal) * kbase^(2*r + 1)), 
							nrow=length(rule.v.decimal), 
							ncol=kbase^(2*r + 1)
							)
	}

	for (i in 1:length(rule.v.decimal)) 
	{
		j <- 1
		rule.decimal <- rule.v.decimal[i]
		while(rule.decimal > 0)
		{
			rule.v.binary[i,j] <- rule.decimal %% kbase	# modular division (remainder)
			j <- j + 1
			rule.decimal <- rule.decimal %/% kbase	# integer division (quocient)
		}
	}

	rule.v.binary
}

###############################################################################
# Neighborhood get
###############################################################################
#
# obtain the neighborhood of x on index
###############################################################################
getNeighborhood <- function(x, index)
{
	neighborhood <- numeric(2 * r + 1)	# zero fill initializer
	n_cells <- length(x)
    # cyclic boundary
    if (BOUNDARY == 2)
    {
        pos <- index - r 	# from left to right on the neighborhood
        if (pos <= 0)		# 
            pos <- n_cells + pos
        for (i in 1:length(neighborhood))
        {
            neighborhood[i] <- x[pos]
            pos <- (pos + 1) %% n_cells
            if (pos == 0)
                pos <- n_cells
        }
    }
    # null or 1 boundary
    else
    {
        pos <- index - r    # from left to right on the neighborhood
        for (i in 1:length(neighborhood))
        {
            if (pos <= 0 || pos > n_cells)
            {
                neighborhood[i] <- BOUNDARY
            }
            else
            {
                neighborhood[i] <- x[pos]
            }
            pos <- pos + 1
        }
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
nextStateGen <- function(x_cur, rule.v=decodeRule(rules)) 
{
	x_next <- numeric(length(x_cur))

	#cat('x_cur ', length(x_cur), '\n')
  	#cat("rule.v",  rule.v, "\n")

  	#	cat("nrow(rule.v)",  nrow(rule.v), "\n")
  	#	cat("ncol(rule.v)",  ncol(rule.v), "\n")

	for (i in 1:length(x_next)) 
	{
		# obtain the neighborhood of x on index i
		neighborhood <- getNeighborhood(x_cur, i)
		#cat('neighborhood',neighborhood,'\n')
		if (totalistic)
		{
			index.neigh <- sum(neighborhood) + 1
		}
        else
        {
            # obtain the index on decimal of the neighborhood
            index.neigh <- 1
            for (j in 0:(2*r)) 
            {
                index.neigh <- index.neigh + k^j * neighborhood[length(neighborhood) - j]
            }
        }

		# obtain which rule is going to use
		# in case of uniform CA, it will be only one rule
		# in case of hybrid CA, it will be more than one
		index.rules <- ((i - 1) %% nrow(rule.v)) + 1

    #	cat("index.rules", index.rules, "\n")
    #	cat("index.neigh", index.neigh, "\n")
    
    #	cat('x_next[i]', rule.v[index.rules,index.neigh], '\n')
		x_next[i] <- rule.v[index.rules,index.neigh]
	}

	x_next
}

###############################################################################
# Evolve
###############################################################################
# 
# evolution of CA for N_iters iterations
# y[iter, cell] <- holds every step on each line
#
###############################################################################
evolve <- function(x_initial, rule=rules, n_iters=N_iters)
{
	# y[N_iters,N_cells] <- 0,...,0
	#				        0,...,0
	#					      ...
	#					    0,...,0
	y <- matrix(rep(0,n_iters*N_cells), nrow=n_iters, ncol=N_cells, byrow=TRUE)
	
	x_aux <- x_initial

	rule.v <- decodeRule(rule)

	for (i in 1:n_iters) {
		y[i,] <- x_aux

		# this is another type of hybrid-CA
		# instead of alternating each cell
		# it alternates on each iteration for using a specified rule
		if(hybrid == "ITER")
		{
			index.rules <- ((i - 1) %% nrow(rule.v)) + 1
			x_next <- nextStateGen(x_aux, matrix(rule.v[index.rules,],nrow=1) )
		}
		else
		{
			x_next <- nextStateGen(x_aux, rule.v)
		}

		x_aux <- x_next
		#cat(y[i,],'\n')
	}

	y
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
			cmbs[i+1,n-j] <- divisor %/% (k**j) # integer division (quocient)
			divisor <- divisor %% (k**j) # modular division (remainder)
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
# Print the cells evolution using plotrix package
###############################################################################
#
# 		mat = matrix(elemento, ncol=N_COLS, byrow=TRUE)
#    	color2D.matplot(mat, main=title)
#
###############################################################################
printCells <- function(y, rule=rules)
{
    title <- c("rule ", paste(rule[1:(1+length(rule)%/%3)], collapse=", "))
    
	color2D.matplot(abs(y-k+1), main=title)
}

###############################################################################
# Initialize Cells
###############################################################################
#
# random = 0 -> 1 cell with 1, and the others are 0
# random = 1 -> random bits to initialize the cells
#
###############################################################################
initializeCells <- function(random = 0, n=N_cells) 
{
	if (random)
	{
		x <- round(runif(n, min=0, max=(k-1)))
	}
	else
	{
		x <- array(0, dim=rep(n,dim))
		if (dim == 1)
			x[n%/%2 + 1] <- 1  		# %/% -> integer division		
		else if (dim == 2)
			x[n%/%2 + 1,n%/%2 + 1] <- 1
		else
			error("dimensions exceeded\n")
	}

	x
}

###############################################################################
# Detect Loop
###############################################################################
#
# Given all states of a Cellular Automata (y)
# Checks if any state repeats more than once
#
###############################################################################
detectLoop <- function(y) 
{
	n_iters <- nrow(y)
	for (i_cur in 1:(n_iters-1)) 
	{
		#x_cur <- y[i_cur]
		for (i_compared in (i_cur+1):n_iters) 
		{
			if (all(y[i_cur,] == y[i_compared,]))
			{
				cat("Wait! I've found a LOOP on:", i_cur, i_compared, "\n")

			}
		}
	}
}

###############################################################################
# Auto Correlation 
###############################################################################
#
# TODO: 
#
###############################################################################
autoCorrelation <- function(y, type)
{
    if (type == "state")
    {
        # auto correlation on columns
        ac <- matrix(numeric(length(y)), ncol=ncol(y), nrow=nrow(y)) 
        for (i in 1:ncol(y))
        {

            # acf function of R
            ac[,i] <- (acf(y[,i], lag.max=nrow(y),plot=FALSE, na.action=na.omit))$acf

            if (PLOT_AUTOCORRELATION)
            {
                filename <- sprintf("plots/autocorrelation-ca1d-evol-cell-%04d.png", i)
                png(filename)
                plot(abs(ac[1:(10*(log(nrow(y),base=10))),i]), 
                    type = "h", xlab = "Lag", ylab = NULL,
                    ylim = c(0,1), main = sprintf("Abs Autocorrelation Cell %d", i)
                )
                dev.off()
            }
        }

        if (PLOT_AUTOCORRELATION)
        {
            system("convert -delay 50 \"plots/*.png\" plots/ca1d-evol-autocorrelation.gif")
            system("rm plots/autocorrelation-ca1d-evol-cell-*.png")
        }
    }
    else if (type == "cell")
    {
        # auto correlation on rows 
        for (i in 1:nrow(y))
        {
        }
    }
    else 
    {
        cat("Error: type not defined.\n")
    }
    ac
}

###############################################################################
# Run Stats
###############################################################################
#
# Obtain the entropy and other statistical measures for analysis
# about the randomness of the cellular automata (y[t,x]) given
#
###############################################################################
runStats <- function(y)
{


    # Entropy of each state (line) during a CA evolution
#    runEntropy(y, n=1, type="state")
#    runEntropy(y, n=2, type="state")
#    runEntropy(y, n=3, type="state")
#    runEntropy(y, n=4, type="state")

    # Correlation between the input and output of the CA
#	cat("---------------------------------------------\n")
#	cat("Correlation (1st, last):", cor(y[1,],y[N_iters,]), "\n")
#	cat("---------------------------------------------\n")

	# Entropy of each cell (column) during a CA evolution
    runEntropy(y, n=1, type="cell")
    runEntropy(y, n=2, type="cell")
    runEntropy(y, n=3, type="cell")
    runEntropy(y, n=4, type="cell")
    runEntropy(y, n=5, type="cell")

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

#    autoCorrelation(y, type="state")

}

###############################################################################
# Run CA
###############################################################################
#
# initialize
# evolve
# ... and print pretty-good-CA for us :)
#
###############################################################################
runCA <- function(rule=rules, seed=NULL, n_iters=N_iters) 
{
	# initialize
	cat("Initializing...\n")
    if (length(seed)==0)
        x <- initializeCells(random.initializer)
    else
        x <- seed

	# evolve
	cat("Evolving...\n")
	y <- evolve(x, rule, n_iters)

    # show some stuff or not
    if (PRINT_CA)
	{
    	cat("Printing...\n")
        printCells(y, rule)
    }

	y
}

###############################################################################
#   Main function
###############################################################################
main <- function() 
{

    #seed <- initializeCells(random.initializer)
    #r <- 1

#    layout(matrix(c(1,1,2,5,3,6,4,7),nrow=2,ncol=4,byrow=FALSE))
    layout(c(1,1,1))
	y <- runCA()

	#runStats(y)

	#cat("Detecting loop...\n")	
	#detectLoop(y)

	cat("Done!\n")
}
###############################################################################
###############################################################################



# main()







# Rule 30 <=> 1Eh
#
# 111 110 101 100 011 010 001 000
#  0   0   0   1   1   1   1   0   => 30
#
# x_center^{t+1} = x_left^t XOR (x_center^t OR x_right^t)
#rule30 <- function(left, center, right) 
#{
#	bitXor(left, bitOr(center,right))
#}




# for (i in 254:255)
# {
#     filename <- sprintf("plots/entropy/ca1d-%04d.png", i)
#     png(filename, bg="white", width=1000, height=700)
#     layout(matrix(c(1,5,2,6,3,7,4,8),nrow=2,ncol=4,byrow=FALSE))
#     y <- runCA(rule=i)
#     runStats(y)
#     dev.off()
# }
# system("convert -delay 50 \"plots/entropy/*.png\" plots/entropy/ca1d-entropy.gif")


