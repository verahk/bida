


library(bida)
library(dplyr)

wd <- "/mn/sarpanitu/ansatte-u6/verahk/Documents/git/bida/"
wd <- setwd(wd)
wd <- "/mn/sarpanitu/ansatte-u6/verahk/Documents/git/bida/"
dir_sim  <- paste0(wd, "inst/simulations/sim_bnlearn/")
dir_data <- paste0(wd, "inst/data/")
dir_simres <- paste0(wd, "../sim_bnlearn/simres/")
#print(paste0(wd, "inst/R/tab_print.R"))
source(paste0(wd, "/inst/R/tab_print.R"))
source(paste0(wd, "/inst/R/labeller.R"))


df <- read.csv("./inst/data/bnlearn_network_params.csv")
names(df) = tolower(names(df))

df$network <- tolower(df$network)
#df <- df[df$nodes < 100, ]
df$desc <- NA
df$descr <- NA
df$card <- NA
df$maxcard <- NA
#df$doprobs <- NA

for (r in 1:nrow(df)){
  bn  <- try(readRDS(paste0(dir_data, df$network[r], ".rds")))
  if (all(class(bn) == "try-error")) next
  n <- length(bn)
  dag  <- bnlearn::amat(bn)
  desc <- bida:::descendants(dag)
  nlevels <- sapply(bn, function(x) dim(x$prob)[1])
  df$desc[r] <- sum(desc)
  df$descr[r] <- df$desc[r]/(n*(n-1))
  df$card[r] <- mean(nlevels)
  df$maxcard[r] <- max(nlevels)
  #df$doprobs[r] <- sum(outer(nlevels, nlevels))

}

caption <- paste("Characteristics of bnlearn networks.",
                   "mbsize = Average Markov blanket size.",
                   "d = Average degree.",
                   "maxin = Maximum indegree.",
                   "card = Average cardinality.",
                   "maxcard = Max cardinality.",
                   "desc = Number of descendants.",
                   "desc = Rate of desc (desc/(n*(n-1)).",
                   "*Water2 correspond to the bnlearn network, but the distribution of all root nodes are positive.")
df <- df %>%
  filter(nodes < 100, maxcard < 15, !network == "water") %>%
  rename(n = nodes, maxin = maxindegree) %>%
  arrange(n, network) %>%
  mutate(across(!c("mbsize", "d", "card", "descr"), ~format(.x, nsmall=0, big.mark=" ") ))

# save cvv
write.csv(df, "./inst/data/network_summary.csv", row.names = F)

# save tex table
df %>%
  group_by(size)  %>%  #-> grouped_df
  tab_to_file(filename = "./inst/data/network_summary",
              type = "latex",
              digits = 2,
              align  = paste0(c("l", "l", "l", rep("r", ncol(df)-2)), collapse = ""),
              caption = caption)


