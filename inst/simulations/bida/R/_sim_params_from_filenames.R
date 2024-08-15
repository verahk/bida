

sim_params_from_filenames <- function(files, names, pattern = "_|\\.", include_params = T){
  df <- data.frame(fileid = as.character(seq_along(files)),
                   do.call("rbind", stringr::str_split(files, pattern)))
  df[ncol(df)] <- NULL
  names(df)[-1] <- names

  if (include_params) {
    stopifnot("network" %in% names)
    params <- read.csv("inst/data/network_summary.csv")
    params <- dplyr::mutate(params,
                            across(where(is.character), ~gsub(" ", "", .x)))
    df <- dplyr::left_join(df, params, by  = "network", multiple = "all")
  }
  return(df)
}
