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
        # TODO:
        # método está ineficiente, pois ele varre o total de subseqs para fazer 1 soma
        # deveria ser O(1) ao invés de O(k^n)

        # método O(1)
        index_aux <- sbsq[i,]
        index <- 0
        for(j in 1:n)
        {
            index <- index + index_aux[j] * (k**(n-j))
        }
        index <- index + 1
        count[index] <- count[index] + 1

        #for (j in 1:nrow(cmbs))
        #{
        #   count[j] <- count[j] + all(sbsq[i,] == cmbs[j,])
        #}
    }

# Contar e eliminar aqueles valores que estão aparecendo pouco
# Montar uma variação da entropia
# 
#    cat('\n')
#    cat(count)
#    cat('\n')

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

    (-sum(log_p*p))/n
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
# Variation of entropy
###############################################################################
#
# x <- CA
# n <- size of subseq 
#      the size must be divisor of the CA length or it will be truncated
#
###############################################################################
entropy.lawful <- function(x, n=1, kbase=k) 
{
    # generate all subsequences of x using size n
    sbsq <- allsubseq(x, n, SUBSEQ_OVERLAP)

    # generate all combinations that is possible with size n and k base
    cmbs <- allcombn(n, kbase)
  
    # count all subseqs
    count <- numeric(nrow(cmbs))
    for (i in 1:nrow(sbsq))
    {
        # TODO:
        # método está ineficiente, pois ele varre o total de subseqs para fazer 1 soma
        # deveria ser O(1) ao invés de O(k^n)

        # método O(1)
        index_aux <- sbsq[i,]
        index <- 0
        for(j in 1:n)
        {
            index <- index + index_aux[j] * (k**(n-j))
        }
        index <- index + 1
        count[index] <- count[index] + 1
    }

    # Remove from count those cmbs that never got achieved
    cmbs.removed <- which(count==0)
    if (length(cmbs.removed) > 0)
        count <- count[-cmbs.removed]

    # probability of each subseq to appear in 'x'
    p <- numeric(nrow(cmbs) - length(cmbs.removed))
    p <- count/(nrow(sbsq))

    # log of probability p
    log_p <- numeric(nrow(cmbs) - length(cmbs.removed))

    for(i in 1:length(log_p))
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

    (-sum(log_p*p))/n
}

###############################################################################
# Cover Degree
###############################################################################
# 
###############################################################################
coverdegree <- function(x, n=1, kbase=k)
{
        # generate all subsequences of x using size n
    sbsq <- allsubseq(x, n, SUBSEQ_OVERLAP)

    # generate all combinations that is possible with size n and k base
    cmbs <- allcombn(n, kbase)
  
    # count all subseqs
    count <- numeric(nrow(cmbs))
    for (i in 1:nrow(sbsq))
    {
        index_aux <- sbsq[i,]
        index <- 0
        for(j in 1:n)
        {
            index <- index + index_aux[j] * (k**(n-j))
        }
        index <- index + 1
        count[index] <- count[index] + 1
    }

    # Remove from count those cmbs that never got achieved
    cmbs.removed <- which(count==0)
    if (length(cmbs.removed) > 0)
        count <- count[-cmbs.removed]

    # cover degree 1-(removed/total)
    #cat("Cover Degree", 1-length(cmbs.removed)/nrow(cmbs), "\n")

    1-length(cmbs.removed)/nrow(cmbs)
}

###############################################################################
# Export Entropy Plots
###############################################################################
#
# exports multiple plots of the evolution of entropy of each cell
#
###############################################################################
exportEntropyPlots <- function(H) 
{
    for (i_cell in 1:N_cells) 
    {
        filename <- sprintf("plots/entropy-ca1d-evol-cell-%04d.png", i_cell)
        png(filename)
        plot(H[,i_cell], main=sprintf("Entropy of cell %d",i_cell))
        dev.off()
    }

    system("convert -delay 50 \"plots/*.png\" plots/ca1d-evol-entropy.gif")
    #system("rm plots/*.png")
}

###############################################################################
# show entropy
###############################################################################
#
#
###############################################################################
show.entropy <- function(H, n=1, plot.main="", xlab="", ylab="Entropy", info="")
{
    cat("---------------------------------------------\n")

    H.min <- min(H)
    H.mean <- mean(H)
    H.max <- max(H)

    plot(H, xlab=xlab, ylab=ylab, main=plot.main, ylim=c(0,1))
    abline(h=1)
    abline(h=H.mean)

    text(0.1*length(H), 0.4, sprintf("Max Entropy: %.4f", H.max), pos=4)
    text(0.1*length(H), 0.3, sprintf("Mean Entropy: %.4f", H.mean), pos=4)
    text(0.1*length(H), 0.2, sprintf("Min Entropy: %.4f", H.min), pos=4)

    cat(info)
    cat("\n")

    cat("Min Entropy:", H.min,'\n')
    cat("at:")
    for (i in 1:length(H)){if(H[i]==H.min) cat(i,' ')}
    cat("\n")

    cat("Mean Entropy:", H.mean,'\n')

    cat("Max Entropy:", H.max,'\n')
    cat("at:")
    for (i in 1:length(H)){if(H[i]==H.max) cat(i,' ')}
    cat("\n")

    cat("---------------------------------------------\n")
}


###############################################################################
# Run Entropy
###############################################################################
#
# Obtain the entropy and other statistical measures for analysis
# about the randomness of the cellular automata (y[t,x]) given
#
###############################################################################
runEntropy <- function(y, n=1, type, lawful=0)
{

    if (length(type)==0)
    {
        cat("Error: Type not defined.\n")
    }

    else
    {
        n_iters <- nrow(y)
        n_cells <- ncol(y)

        # Entropy of each state (line) during a CA evolution
        if (type == "state")
        {
            H <- numeric(n_iters)
            for (i in 1:n_iters) 
            {
                H[i] <- entropy(y[i,], n)
            }

            if (PRINT_ENTROPY)
            {
                plot.main=sprintf("Entropy of subseq size %d",n)
                info.main=sprintf("Entropy of subseq size %d on States",n)
                axis.x.lab <- "iters"
                show.entropy(H=H,n=n, plot.main=plot.main, xlab=axis.x.lab, ylab="Entropy",
                    info=info.main)
            }
        }

        # Entropy of each cell (column) during a CA evolution
        else if (type == "cell")
        {
            H <- numeric(n_cells)
            if(lawful)
            {
                for (i in 1:n_cells) 
                {
                    H[i] <- entropy.lawful(y[,i], n)
                }
            }
            else
            {
                for (i in 1:n_cells) 
                {
                    H[i] <- entropy(y[,i], n)
                }
            }

            if (PRINT_ENTROPY)
            {
                plot.main=sprintf("Entropy of subseq size %d",n)
                info.main=sprintf("Entropy of subseq size %d on Cells",n)
                axis.x.lab <- "cell"
                show.entropy(H=H,n=n, plot.main=plot.main, xlab=axis.x.lab, ylab="Entropy",
                    info=info.main)
            }
        }
    }

    H
}


###############################################################################
# Run Cover Degree
###############################################################################
#
# Obtain the cover degree for analysis
# of the cellular automata (y[t,x]) given
#
###############################################################################
runCover <- function(y, n=1, type)
{

    if (length(type)==0)
    {
        cat("Error: Type not defined.\n")
    }

    else
    {
        n_iters <- nrow(y)
        n_cells <- ncol(y)

        # Cover degree of each state (line) during a CA evolution
        if (type == "state")
        {
            CD <- numeric(n_iters)
            for (i in 1:n_iters) 
            {
                CD[i] <- coverdegree(y[i,], n, k)
            }
        }
        # Cover degree of each cell (column) during a CA evolution
        else if (type == "cell")
        {
            CD <- numeric(n_cells)
            for (i in 1:n_cells) 
            {
                CD[i] <- coverdegree(y[,i], n, k)
            }
        }
    }

    CD.min <- min(CD)
    CD.mean <- mean(CD)
    CD.max <- max(CD)

    #cat("---------------------------------------------\n")
    #cat("Min Cover:", CD.min,'\n')
    #cat("Mean Cover:", CD.mean,'\n')
    #cat("Max Cover:", CD.max,'\n')
    #cat("---------------------------------------------\n")

    #plot(CD, main=sprintf("cover degree, hmax=%d",n), ylim=c(0,1))

    CD
}

