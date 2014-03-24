library(stringr)

source('ca.R')

double.re <- "([[:digit:]]+)[[:punct:]]?([[:digit:]]*)"

#FILENAME="data-teste.txt"
#FILENAME="data-teste2.txt"
FILENAME="data_t_0_k_2_r_1.txt"
#FILENAME="data_t_1_k_3_r_1.txt"
#FILENAME="data_t_1_k_2_r_3.txt"

getValues <- function(data.v)
{
    regexp <- double.re
    data.v <- str_extract_all(data.v, regexp)
    data.v
}

getDataFromFile <- function(filename)
{
    data.file <- scan(filename, sep='\n', what="")
}

getFitness <- function(stats.m, h.max)
{
    fitness.v <- numeric(nrow(stats.m))
    for(i in 1:h.max){fitness.v <- fitness.v + i*stats.m[,i]}
    fitness.v <- fitness.v / (sum(1:h.max))
    fitness.v
}

getFitnessGeom <- function(stats.m, h.max)
{
    fitness.v <- numeric(nrow(stats.m))
    fitness.v <- stats.m[,1]
    if (h.max > 1)
    {
        for(i in 2:h.max){fitness.v <- fitness.v * stats.m[,i]}
    }
    fitness.v
}

data.raw <- getDataFromFile(FILENAME)

cat("----\n")
data.v <- data.raw[grep("Fitness", data.raw)]
fitness.v <- as.double(getValues(data.v ))

cat("----\n")
data.v <- data.raw[grep("Stats", data.raw)]
stats.v <- getValues(data.v)
stats.m <- matrix(numeric(5*length(stats.v)), ncol=5)
for (i in 1:length(stats.v))
{
    stats.m[i,] <- matrix(as.double(stats.v[[i]]), ncol=5)
}

cat("----\n")
fit.m <- matrix(numeric(5*nrow(stats.m)), ncol=5)
for (i in 1:5)
{
    fit.m[,i] <- getFitness(stats.m, i)
}

cat("----\n")
fitgeo.m <- matrix(numeric(5*nrow(stats.m)), ncol=5)
for (i in 1:5)
{
    fitgeo.m[,i] <- getFitnessGeom(stats.m, i)
}

#
# presents the behavior of curves of metric M and E_hmax for hmax = {1:5}
#
layout(matrix(c(1:5),nrow=5,byrow=TRUE))
for(i in 1:5)
{ 
    filename<-sprintf("plots/ca1d-t-%d-k-%d-r-%d-%d.png",totalistic,k,r,i)

    png(filename)
    
    plot(sort(fit.m[,i]), pch=19, main=sprintf("hmax = %d",i), ylab="Valor da métrica", xlab="Regras ordenadas pelo valor das métricas", ylim=c(0,1))
    points(sort(stats.m[,i]), pch=22, col="green")
    points(sort(fitgeo.m[,i]), pch=24, col="red")
    legend("topleft", inset=.05, title="Tipo de curva",c("E(hmax)","Ma(hmax)", "Mg(hmax)"), fill=c("green", "black", "red"), horiz=FALSE)
    dev.off()
}

n_rules = nrow(stats.m)

#
# presents the behavior of TOP10 rules on metric M and E_hmax, for hmax = {1:5}
#
layout(matrix(c(1:5),nrow=5,byrow=TRUE))
for(i in 1:5)
{ 
    filename<-sprintf("plots/ca1d-zoom-t-%d-k-%d-r-%d-%d.png",totalistic,k,r,i)

    png(filename)
    
    y<-sort(fit.m[,i])[(n_rules-9):n_rules]
    x<-list()
    for (j in 1:10)
    {
        x <- append(x, paste(which(y[j]==fit.m[,i]), collapse=" "))
    }

    plot(y=y, x=1:10, main=sprintf("hmax = %d",i), ylab="M", xlab="Top 10 Regras", ylim=c(0,1), xaxt="n",
        axes=FALSE)

    # Make x axis tick marks without labels
    axis(1, at=1:10, format(x), las=2, cex.axis=0.5)

    # Create box around plot
    box()

    y<-sort(stats.m[,i])[(n_rules-9):n_rules]
    x<-list()
    for (j in 1:10)
    {
        x <- append(x, which(y[j]==stats.m[,i])[1])
    }
    points(y=y, x=1:10, pch=16, col="green")
    legend("bottomright", inset=.05, title="Type of curve",c("E_hmax","M"), fill=c("green", "black"), horiz=FALSE)
    dev.off()
}



#system("convert -delay 50 \"plots/ca1d-t-*.png\" plots/ca1d-evol-autocorrelation.gif")
#system("rm plots/*.png")

cat('\n')
cat("Top 10: M Arit")
cat('\n')
for (i in sort(fit.m[,5], decreasing=TRUE)[1:10])
{
    cat(which(i==fit.m[,5], arr.ind=TRUE)," ", i)
    cat('\n')
}

cat('\n')
cat("Top 10: M Geom")
cat('\n')
for (i in sort(fitgeo.m[,5], decreasing=TRUE)[1:10])
{
    cat(which(i==fitgeo.m[,5], arr.ind=TRUE)," ", i)
    cat('\n')
}



for (h in 1:5)
{
    cat('\n')
    cat("Top 10: h_max = ", h)
    cat('\n')
    for (i in sort(stats.m[,h], decreasing=TRUE)[1:10])
    {
        cat(which(i==stats.m[,h], arr.ind=TRUE)," ", i)
        cat('\n')
    }
}

