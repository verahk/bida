
parent_support_exact <- function(data, type, max_parent_size, temp_id, filepath_apssolver, th = .999, ess = 1, nlev){
  n <- ncol(data)
  if (n > 20){
    stop("Don't use on systems with more than 20 variables.")
  }

  # compute family scores and write to file
  if (tolower(type) == "categorical") {
    filepath_score   <- paste("temp",temp_id,".bdeu.score", sep="")
    filepath_support <- paste("temp",temp_id,".bdeu.support", sep="")
    famscore_to_file_bdeu(data,  max_parent_size, file_out=filepath_score, ess = ess, nlev = nlev)
  } else if (tolower(type) == "gaussian") {
    filepath_score   <- paste("temp", temp_id,".fml.score", sep="")
    filepath_support <- paste("temp", temp_id,".fml.support", sep="")
    famscore_to_file_fml(data, max_parent_size, file_out=filepath_score)
  } else {
    stop("Data type must be either Cateogorical or Gaussian!")
  }

  on.exit(file.remove(filepath_score, filepath_support))

  # compute parent support
  aps_type <- "modular"
  system(paste(paste0(filepath_apssolver, "/aps"), aps_type, filepath_score, filepath_support, sep=" "))
  ps <- read_in_aps_parent_support(filepath_support, max_parent_size)

  # store sets drop parent sets with insuff support
  for (x in seq_len(n)) {
    if (th < 1 && any(ps$support[[x]] < 1-th)) {
      sorted <- sort.int(ps$support[[x]], decreasing = TRUE, index.return = TRUE)
      indx  <- sorted$ix[cumsum(sorted$x) <= th]
      ps$sets[[x]] <-  ps$sets[[x]][indx,, drop = FALSE]
      ps$support[[x]] <- ps$support[[x]][indx]/sum(ps$support[[x]][indx])
    }
  }
  return(ps)
}
