
#fn <- "script_cenapad_30_141516171819.Rout"
#fn <- "script_cenapad_30_202122.Rout"
#fn <- "script_cenapad_30_23.Rout"
fn <- "script_cenapad_30_24.Rout"

#fn <- "script_cenapad20212223.Rout"
#fn <- "script_cenapad151719.Rout"
#fn <- "script_cenapad141618.Rout"
fn <- "script_cenapad_45_24.Rout"

data <- scan(fn, what=character())

i_20 <- which(data=="main(20)")
i_21 <- which(data=="main(21)")
i_22 <- which(data=="main(22)")
i_23 <- which(data=="main(23)")
i_24 <- which(data=="main(24)")


i_14 <- which(data=="main(14)")
i_16 <- which(data=="main(16)")
i_18 <- which(data=="main(18)")

i_15 <- which(data=="main(15)")
i_17 <- which(data=="main(17)")
i_19 <- which(data=="main(19)")

i_ciclos <- which(data=="ciclo")+2
#i_ciclos <- i_ciclos[i_ciclos>i_17 & i_ciclos<i_19]
i_ciclos <- i_ciclos[i_ciclos>i_24]
ciclos <- as.integer(data[i_ciclos])
summary(ciclos)
length(ciclos)




> i_ciclos <- i_ciclos[i_ciclos>i_14 & i_ciclos<i_16]
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    2.0    84.0   168.0   443.8   534.0  2198.0 

> i_ciclos <- i_ciclos[i_ciclos>i_15 & i_ciclos<i_17]
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    1.0    32.0    60.0   606.8    60.0  6820.0 

> i_ciclos <- i_ciclos[i_ciclos>i_16 & i_ciclos<i_18]
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    2.0    48.0    48.0   241.8   296.0  2816.0 

> i_ciclos <- i_ciclos[i_ciclos>i_17 & i_ciclos<i_19]
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    2.0   255.0   408.0  5958.0   854.2 78810.0 

> i_ciclos <- i_ciclos[i_ciclos>i_18]
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1.00   18.00   36.00   31.79   36.00  504.00 

> i_ciclos <- i_ciclos[i_ciclos>i_19]
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      2     114     456   14560     912  183900 

> i_ciclos <- i_ciclos[i_ciclos>i_20 & i_ciclos<i_21]
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    2.0    60.0    60.0   274.8   120.0  9112.0 

> i_ciclos <- i_ciclos[i_ciclos>i_21 & i_ciclos<i_22]
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      1      42      84    8004     168  352900 

> i_ciclos <- i_ciclos[i_ciclos>i_22 & i_ciclos<i_23]
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      2     132     264    1199     264  122900 

> i_ciclos <- i_ciclos[i_ciclos>i_23]
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      2     276     276   37450     552 3460000 



