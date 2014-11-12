k <- read.table("k562.chromhmmclass")
g <- read.table("gm12878.chromhmmclass")

sum(k[,6] == 1 & k[,7] == 0)/NROW(k) ## Overlaps *just* promoter.
sum(k[,6] == 0 & k[,7] == 1)/NROW(k) ## Overlaps enahncer.
sum(k[,6] == 1 & k[,7] == 1)/NROW(k) ## Overlaps both promoter and enhancer.
sum(k[,6] == 0 & k[,7] == 0)/NROW(k) ## Neither/ other.
sum(k[,6] == 0 & k[,7] == 0 &  k[,9] == 1)/NROW(k) ## Trans.
sum(k[,6] == 0 & k[,7] == 0 &  k[,9] == 0 &  k[,10] == 1)/NROW(k) ## Het.

sum(g[,6] == 1 & g[,7] == 0)/NROW(g) ## Overlaps *just* promoter.
sum(g[,6] == 0 & g[,7] == 1)/NROW(g) ## Overlaps *just* enhancer.
sum(g[,6] == 1 & g[,7] == 1)/NROW(g) ## Overlaps both promoter and enhancer.
sum(g[,6] == 0 & g[,7] == 0)/NROW(g) ## Neither/ other.
sum(g[,6] == 0 & g[,7] == 0 &  g[,9] == 1)/NROW(g) ## Trans.
sum(g[,6] == 0 & g[,7] == 0 &  g[,9] == 0 &  g[,10] == 1)/NROW(g) ## Het.

ka <- c("Promoter"= sum(k[,6])/ NROW(k),
  "Enhancer_(excluding_promoter)"= sum(k[,6] == 0 & k[,7] == 1)/ NROW(k),
  "Trans"= sum(k[,6] == 0 & k[,7] == 0 &  k[,9] == 1)/NROW(g),
  "Other"= sum(k[,6] == 0 & k[,7] == 0 &  k[,9] == 0)/ NROW(k))
			
ga <- c("Promoter"= sum(g[,6])/ NROW(g),
  "Enhancer_(excluding_promoter)"= sum(g[,6] == 0 & g[,7] == 1)/ NROW(g),
  "Trans"= sum(g[,6] == 0 & g[,7] == 0 &  g[,9] == 1)/NROW(g),
  "Other"= sum(g[,6] == 0 & g[,7] == 0 &  g[,9] == 0)/ NROW(g))
			
df <- data.frame(k562= ka, gm12878= ga)
df



## < 0.5% are *JUST* insulators.
sum(k[,8] == 1 & k[,6] == 0 & k[,7] == 0)/NROW(k)
sum(g[,8] == 1 & g[,6] == 0 & g[,7] == 0)/NROW(g)