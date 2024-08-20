

tag_from_sim_settings <- function(bnname, N, r, prefix = NULL, postfix = NULL) {
  tag <- sprintf("%s_N%s_r%02.0f", bnname, N, r)
  paste(prefix, tag, postfix, sep = "")
}

