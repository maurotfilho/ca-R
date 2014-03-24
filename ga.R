source("ca.R")

#goodrules <- c(30 , 45)
#goodrules <- c(196 , 532 , 682 , 920 , 1263 , 1265 , 1600 , 1747 , 1770)
goodrules <- NULL

###############################################################################
# Definitions
###############################################################################
FITNESS_MIN = 0.9
sink_filename = "cover_t_%d_k_%d_r_%d.txt"
random.test = 1

###############################################################################


###############################################################################
# do test
###############################################################################
#
# run the tests
#
###############################################################################
doTest <- function(rule.v=goodrules, seed=NULL, type="entropy")
{
    cat("Type:", type, "\n")
    good.rules <- list()

    if (totalistic)
        n_rules <- k**(1 + (k-1)*(2*r+1))
    else
        n_rules <- k**(k**(2*r+1))
    
    fitness <- numeric(n_rules)

    # if no seed is given, it loads the seed with random entry
    if(length(seed)==0)
    {
        seed <- initializeCells(1)
    }

    # if no rules is given, it loads with every rule! bruteforce OVER 9000!
    if(length(rule.v)==0)
    {
        rule.v <- 1:(n_rules-1)
        cat("Rules:", n_rules,"\n")
    }

    stats.v <- matrix(7*numeric(n_rules), nrow=n_rules, ncol=7)

    # do the tests for every rule in rule vector (rule.v)
    for (index_rule in rule.v)
    {
        cat("---------------------------------------------\n")
        cat("Rule", index_rule, '\n')

        # runs the CA for rule given and seed
        # using random seed ! ! !
        y <- runCA(rule=index_rule, seed=initializeCells(1))

        cat("Evaluating...\n")

        #ac <- autoCorrelation(y, type="state")

        # vector that stores the statistical values obtained        
        #stats <- numeric(6)
        #stats[1] <- 1-max((abs(ac[2:(10*(log(length(ac),base=10)))])))
        #if (is.nan(stats[1])){ stats[1] <- 0 }
        if (type=="entropy")
        {
            stats <- numeric(6)
            stats[1] <- min(runEntropy(y, n=1, type="cell", lawful=1))
            stats[2] <- min(runEntropy(y, n=2, type="cell", lawful=1))
            stats[3] <- min(runEntropy(y, n=3, type="cell", lawful=1))
            stats[4] <- min(runEntropy(y, n=4, type="cell", lawful=1))
            stats[5] <- min(runEntropy(y, n=5, type="cell", lawful=1))
            stats[6] <- min(runEntropy(y, n=6, type="cell", lawful=1))
            #stats[1] <- mean(runEntropy(y, n=1, type="state"))
        }
        else if (type == "coverdegree")
        {
            stats <- numeric(6)
            stats[1] <- min(runCover(y, n=1, type="cell"))
            stats[2] <- min(runCover(y, n=2, type="cell"))
            stats[3] <- min(runCover(y, n=3, type="cell"))
            stats[4] <- min(runCover(y, n=4, type="cell"))
            stats[5] <- min(runCover(y, n=5, type="cell"))
            stats[6] <- min(runCover(y, n=6, type="cell"))
        }
        # fitness value
        fitness.aux <- 0
        for(i in 1:length(stats)){fitness.aux <- fitness.aux + i*stats[i]}
        fitness[index_rule] <- fitness.aux/(sum(1:length(stats)))

        cat('Stats', stats, '\n')
        cat('Fitness', fitness[index_rule],'\n')

        # R U GUD ENOUGH?! .o) _o.
        if (fitness[index_rule] >= FITNESS_MIN)
        {
            good.rules <- append(good.rules, index_rule)
        }

        stats.v[index_rule,] <- c(fitness[index_rule],stats)

        cat("---------------------------------------------\n")
    }

    # No rules that are good enough
    if (length(good.rules) == 0)
    {
        cat("No rules!\n")
    }
    # AWWWWWWWWWWWYEEAAAAAAAAH
    else 
    {
        for (i in 1:length(good.rules))
        {
            cat(good.rules[[i]], ', ')
        }
    }
    cat('max fitness', max(fitness),'\n')

    stats.v
}

###############################################################################
# Run Tests
###############################################################################
#
# Here i can put so much tests
# and the output from ones can be the input for others!
#
###############################################################################
runTests <- function()
{
    sink(sprintf(sink_filename,totalistic,k,r), append=FALSE, split=FALSE)
    doTest(seed=initializeCells(random.test), type="coverdegree")
    sink()
}


###############################################################################
#   Main function
###############################################################################
main <- function() 
{
    runTests()
}

###############################################################################
###############################################################################

main()




# goodrules <- NULL

# CONVENTIONAL
# rules.v <- c(15 , 30 , 42 , 45 , 75 , 85 , 86 , 89 , 101 , 102 , 105 , 106 , 112 , 120 , 135 , 
# 138 , 149 , 150 , 153 , 154 , 165 , 166 , 169 , 170 , 174 , 180 , 208 , 210 , 225 , 240 , 241 
# )
# goodrules <- c(30 , 45 , 75 , 86 , 89 , 101 , 135 , 149)

# TOTALISTIC
#  fitness 6
# goodrules <- c(46 , 66 , 92 , 96 , 102 , 138 , 142 , 146 , 147 , 150 , 
#    177 , 194 , 196 , 204 , 223 , 277 , 289 , 307 , 309 , 311 , 329 , 
#    335 , 339 , 340 , 345 , 347 , 385 , 388 , 389 , 416 , 424 , 439 , 
#   469 , 470 , 471 , 473 , 520 , 532 , 538 , 550 , 552 , 554 , 572 , 
#   578 , 581 , 582 , 583 , 586 , 587 , 592 , 624 , 626 , 628 , 632 , 
#   633 , 635 , 636 , 640 , 642 , 659 , 662 , 663 , 682 , 690 , 709 , 
#   750 , 763 , 775 , 776 , 777 , 788 , 793 , 794 , 795 , 797 , 825 , 
#   831 , 844 , 866 , 867 , 869 , 876 , 893 , 896 , 902 , 905 , 906 , 
#   915 , 920 , 923 , 933 , 938 , 939 , 942 , 960 , 983 , 1018 , 1031 , 1033 , 1037 , 1038 , 
#   1074 , 1077 , 1078 , 1136 , 1139 , 1149 , 1150 , 1155 , 1158 , 1226 , 1230 , 1236 , 1262 , 
#   1263 , 1265 , 1266 , 1267 , 1274 , 1279 , 1280 , 1281 , 1283 , 1301 , 1310 , 1311 , 1315 , 
#   1316 , 1317 , 1319 , 1320 , 1321 , 1379 , 1382 , 1392 , 1492 , 1504 , 1522 , 1524 , 1527 , 
#   1544 , 1546 , 1550 , 1553 , 1554 , 1558 , 1560 , 1562 , 1563 , 1568 , 1573 , 1594 , 1595 , 
#   1596 , 1600 , 1603 , 1604 , 1605 , 1607 , 1608 , 1609 , 1612 , 1614 , 1617 , 1635 , 1663 , 
#   1712 , 1716 , 1735 , 1747 , 1753 , 1765 , 1767 , 1769 , 1770 , 1771 , 1774 , 1780 , 1787 , 
#   1793 , 1797 , 1798 , 1801 , 1802 , 1803 , 1805 , 1806 , 1807 , 1874 , 1879 , 1959 , 1978 , 
#   1990 , 1996 , 2010 , 2032 , 2040 , 2041)

# fitness 7
# goodrules <- c(102 , 146 , 177 , 196 , 204 , 289 , 307 , 309 , 339 , 345 , 347 , 388 , 389 , 
#     416 , 424 , 469 , 532 , 538 , 582 , 583 , 586 , 592 , 626 , 633 , 635 , 642 , 659 , 663 , 
#     682 , 763 , 775 , 776 , 788 , 793 , 794 , 795 , 825 , 831 , 844 , 869 , 902 , 906 , 920 , 
#     923 , 933 , 938 , 939 , 942 , 1031 , 1074 , 1149 , 1262 , 1263 , 1265 , 1266 , 1267 , 1274 , 
#     1279 , 1280 , 1281 , 1301 , 1315 , 1317 , 1319 , 1320 , 1321 , 1379 , 1492 , 1524 , 1554 , 
#     1560 , 1563 , 1573 , 1595 , 1596 , 1600 , 1604 , 1605 , 1608 , 1609 , 1614 , 1712 , 1716 , 
#     1735 , 1747 , 1753 , 1765 , 1767 , 1769 , 1770 , 1771 , 1774 , 1803 , 1805 , 1959 , 1978 , 
#     1996 , 2010 , 2032 , 2040)

# fitness 8
#goodrules <- c(196 , 532 , 682 , 920 , 1263 , 1265 , 1600 , 1747 , 1770)


# totalistic, r=3, k=2, 0-255
# max fitness 0.9663076
#41 , 45 , 54 , 90 , 105 , 113 , 142 , 150 , 165 , 210

# totalistic, r=1, k=3, 0-2186
