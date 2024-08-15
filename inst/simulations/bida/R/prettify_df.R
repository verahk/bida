prettify_df <- function(df){


   # lookup tables that give pretty names and order variables
  lookup <- list()

    lookup$adjset <- matrix(c(#"exact", "pa", # bidaexact
      "parents","pa",
      "pa", "pa",
      "parentsmin","pa-min",
      "pamin", "pa-min",
      "o", "o",
      "o*", "o*",
      "osetmin","o-min",
      "omin", "o-min",
      "oset","o",
      "cond","cond.",
      "marg","marg.",
      "ancprob", "ancp",
      "-", "-"),
      ncol = 2, byrow = T)

    lookup$network <- matrix(c("asia", "Asia (n=8)",
                               "sachs", "Sachs (n=11)",
                               "child", "Child (n=20)",
                               "insurance", "Insurance (n=27)",
                               "water2", "Water2 (n=32)",
                               "alarm", "Alarm (n=37)",
                               "hailfinder", "Hailfinder (n=56)",
                               "hepar2", "Hepar2 (n=70)",
                               "win95pts", "Win95pts (n=76)"),
                                ncol = 2, byrow = T)

  if ("N" %in% names(df)) {
    df <- mutate(df,
                 N = as.numeric(gsub("N", "", N)),
                 N = factor(N, sort(unique(N))))
  }


  for (var in c("adjset", "network")) {
    if (var %in% names(df)) {
      x <- gsub("_", "", df[[var]]) # variable with any "_" removed
      df[[var]] <- factor(x, lookup[[var]][, 1], lookup[[var]][, 2])
    }
  }

  return(df)
}




