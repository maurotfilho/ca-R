
source('ca.R')
source('cycle-detection.R')

k<-2
r<-1
totalistic<-0
N_cells <- 25
rules<-30

filename <- sprintf('qua-n%d-rule%d.csv', N_cells, rules)
cycles <- main(N_cells)

d <- data.frame(cycles)

write.csv(d, file=filename)

