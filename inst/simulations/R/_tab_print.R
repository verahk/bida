

tab_to_file <- function(grouped_df,
                        filename,
                        type = "latex",
                        digits = 1,
                        caption = "",
                        fbold = NULL,
                        byrow = T,
                        align = NULL)
{

  if (!is.null(fbold)) {
    tabs <- tab_print_bold(grouped_df, byrow = byrow, fbold = fbold, digits = digits, type = type)
  } else {
    tabs <- setNames(rep(list(grouped_df), length(type)), type)
  }

  for (t in type){
    xtab <- xtable::xtable(tabs[[t]],
                           digits = digits,
                           align = align,
                           caption = caption)
    if (t == "latex"){
      addtorow <- tab_add_space_by_row(dplyr::group_indices(grouped_df))
    } else {
      addtorow <- NULL
    }


    print(xtab,
          sanitize.text.function=identity,
          width = "\\textwidth",
          type = t,
          file = paste0(filename, ".", gsub("^la", "", t)),
          tabular.environment = "tabularx",
          rotate.colnames = t == "latex" && ncol(xtab) > 15,
          include.rownames = FALSE,
          booktabs = TRUE,
          size = "\\footnotesize",
          html.table.attributes = "border = 0",
          add.to.row = addtorow)
  }
}


tab_print_bold <- function(df, fbold, byrow = T, indx_values = sapply(df[1, ], class) == "numeric", digits = 1,  type = c("latex", "html"))
{
  stopifnot(all(sapply(df[indx_values], class) == "numeric"))
  if (length(digits) == 1) digits <- rep(digits, ncol(df)+1)

  mat  <- as.matrix(df[, indx_values])
  bold <- mat*0

  if (byrow) {
    for (i in 1:nrow(mat)){
      pos <- fbold(mat[i, ])
      bold[i, pos] <- 1
    }
  } else {
    for (i in 1:ncol(mat)){
      pos <- fbold(mat[, i])
      bold[pos, i] <- 1
    }
  }



  for (i in 1:ncol(mat)) mat[, i] <- round(mat[, i], digits[1+which(indx_values)][i])
  mode(mat) <- "character"

  tabs <- list()
  indx <- which(bold == 1)
  if ("latex" %in% type)
  {
    tabs$latex <- df
    tabs$latex[indx_values] <- replace(mat, indx, paste("\\textbf{", mat[indx], "}", sep = ""))
  }
  if ("html" %in% type)
  {
    tabs$html  <- cbind(df[, !indx_values],
                        replace(mat, indx, paste("<strong>", mat[indx], "</strong>", sep = "")))
  }
  return(tabs)
}

tab_add_space_by_row <- function(groups, startpos = 0) {
  addtorow <- list()
  if (length(unique(groups)) > 0){
    uniq  <- unique(groups)
    grpos <- startpos + cumsum(tabulate(groups)[uniq[-length(uniq)]])

    addtorow$pos <- as.list(grpos)
    addtorow$command <- as.vector(rep("\\addlinespace[0.5em]", length(grpos)), mode = "character")
  }
  return(addtorow)
}


#groups <- lapply(c("^pc", "MCMCPC", "MCMCPC", "MCMCCP"), function(x) which(grepl(x, names(tab))))

# tab_print_bold <- function(df, groups, fbold = which.max, text = "bf")
# {
#
#   if (length(fbold) == 1) rep(list(fbold), length(groups))
#   if (length(text) == 1)  rep(text, length(groups))
#
#   mat  <- as.matrix()
#   bold <- matrix(nrow = nrow(df), ncol = ncol(df))
#
#   for (g in seq_along(groups)) {
#       group <- groups[[g]]
#       for (i in 1:nrow(df)) {
#         pos <- fbold[[g]](df[i, group])
#         bold[i, group[pos]] <- text[[g]]
#     }
#   }
#
#   tmp <- matrix(c("\\textbf{", "}",
#                    "\\textit{", "}"), 2, 2,
#                  byrow = T,
#                  dimnames = list(c("bf", "it"), NULL))
#
#   tab <- mutate(df, across(where(is.numeric), ~ as.character(round(.x, 2))))
#   tab <-
#   for (i in 1:nrow(df)){
#     for (j in 1:nrow(cold)) {
#       if (!is.na(bold[i, j])){
#
#       }
#     }
#   }
#
#   if (byrow) {
#     for (i in 1:nrow(mat)){
#       pos <- fbold(mat[i, ])
#       bold[i, pos] <- 1
#     }
#   } else {
#     for (i in 1:ncol(mat)){
#       pos <- fbold(mat[, i])
#       bold[pos, i] <- 1
#     }
#   }


#
#   for (i in 1:ncol(mat)) mat[, i] <- round(mat[, i], digits[1+which(indx_values)][i])
#   mode(mat) <- "character"
#
#   tabs <- list()
#   indx <- which(bold == 1)
#   if ("latex" %in% type)
#   {
#     tabs$latex <- df
#     tabs$latex[indx_values] <- replace(mat, indx, paste("\\textbf{", mat[indx], "}", sep = ""))
#   }
#   if ("html" %in% type)
#   {
#     tabs$html  <- cbind(df[, !indx_values],
#                         replace(mat, indx, paste("<strong>", mat[indx], "</strong>", sep = "")))
#   }
#   return(tabs)
# }

tab_add_space_by_row <- function(groups, startpos = 0) {
  addtorow <- list()
  if (length(unique(groups)) > 0){
    uniq  <- unique(groups)
    grpos <- startpos + cumsum(tabulate(groups)[uniq[-length(uniq)]])

    addtorow$pos <- as.list(grpos)
    addtorow$command <- as.vector(rep("\\addlinespace[0.5em]", length(grpos)), mode = "character")
  }
  return(addtorow)
}


## helper-function for xtable:
# ## add space to rows based on a variable group
# ## return nested list with
# addtorow_linespace_by_group <- functions(pos, x = NULL) {
#   if (is.null(x)) x <- list(pos = vector("list", 0L),
#                             command = vector("character", 0L))
#   if (length(unique(groups)) > 0){
#     x$pos <- c(x$pos, as.list(grpos))
#     x$command <- c(x$command, cmd)
#     #x$command <- as.vector(rep("\\addlinespace[0.5em]", length(grpos)), mode = "character")
#   }
#   return(x)
# }


addtorow_linespace_by_group <- function(groups,startpos = 0) {
   x <- list(pos = vector("list", 0L),
             command = vector("character", 0L))
  if (length(unique(groups)) > 0){
    uniq  <- unique(groups)
    grpos <- startpos + cumsum(tabulate(groups)[uniq[-length(uniq)]])
    cmd   <- rep("\\addlinespace[0.5em]", length(grpos))

    x$pos <- c(x$pos, as.list(grpos))
    x$command <- c(x$command, cmd)
    #x$command <- as.vector(rep("\\addlinespace[0.5em]", length(grpos)), mode = "character")
  }
  return(x)
}

## helper-function for xtable:
## add multiheader by pasting unique combs of columns to a character vector
## input:
## - keys: character vector with key variables (id-variables in wide format)
## - cols: a data frame with 2 columns, with all combinations of the two variables

tab_multi_header <- function(df, to_cols = NULL, sep_by = NULL, caption = NULL, label = NULL, filename = NULL, sideways = F, standalone = T) {
  df <- arrange(df, across(all_of(to_cols)))
  tab <- tidyr::pivot_wider(df, names_from = to_cols)
  add.to.row <- list(pos = vector("list", 0L),
                     command = vector("character", 0L))

  if (length(to_cols) > 0) {
    keys   <- names(select(df, -all_of(to_cols), -value))
    pos <- list(-1, nrow(tab))
    command <- c(addtorow_multi_header(keys, distinct(select(ungroup(df), all_of(to_cols))), sideways = ncol(tab) > 15),
                    "\\bottomrule\n")

    add.to.row$pos <- c(add.to.row$pos, pos)
    add.to.row$command <- c(add.to.row$command, command)
  }

  if (!is.null(sep_by)) {
    tab <- select(tab, all_of(sep_by), everything())
    groups <- group_indices(group_by(tab, across(all_of(sep_by))))
    tmp    <- addtorow_linespace_by_group(groups, 0)

    add.to.row$pos <- c(add.to.row$pos, tmp$pos)
    add.to.row$command <- c(add.to.row$command, tmp$command)
  }

   xtab <- xtable::xtable(tab, caption = caption, label = label, digits = 0)



  if (!is.null(filename)) {


    print(xtab,
          #sanitize.text.function=identity,
          width = "\\textwidth",
          type = "latex",
          table.placement = "h",
          file = filename,
          append = F,
          tabular.environment = "tabularx",
          #rotate.colnames = ncol(tab) > 10,
          #include.rownames = FALSE,
          #header
          include.rownames = FALSE, #Don't print rownames
          include.colnames = is.null(to_cols), #We create them ourselves
          caption.placement = "top", #"top", NULL
          hline.after=NULL, #We don't need hline; we use booktabs
          floating=TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
          sanitize.text.function = force, # Important to treat content of first column as latex function
          booktabs = TRUE,
          size = "\\footnotesize",
          html.table.attributes = "border = 0",
          add.to.row = add.to.row)

    filename_preview <- gsub("\\.tex", "_preview.tex", filename)
    file.create(filename_preview)
    cat(c("\\documentclass[]{article}",
          "\\usepackage[margin = 2cm]{geometry}",
          "\\usepackage[utf8]{inputenc}",
          "\\usepackage{tabularx}",
          "\\usepackage{booktabs}",
          "\\usepackage{caption}",
          "\\usepackage{rotating}",
          "\\usepackage{lscape}",
          "\\begin{document}",
          paste0("\\input{", gsub(".*\\/", "", filename), "}"),
          "\\end{document}"),
        sep = "\n",
        file = filename_preview,
        append = F)

  } else {
    out <- list(xtab = xtab,
                add.to.row = add.to.row)
    return(out)
  }
}

addtorow_multi_header <- function(keys, cols, sideways = T) {
  cols <- lapply(cols, as.character)
  nkeyvars <- length(keys)
  addrow <- list("\\toprule \n",
                 labs1 = list(),
                 midrule = list())
  if (sideways) cols[[2]] <- paste("\\begin{sideways}", cols[[2]], "\\end{sideways}")
  addrow$labs1[[1]] <- sprintf("\\multicolumn{%i}{c}{}", nkeyvars)
  cumcount <- nkeyvars
  for (x in unique(cols[[1]])) {
    count <- sum(cols[[1]] == x)
    addrow$labs1[[x]] <- paste0(sprintf("\\multicolumn{%i}{c}", count), "{", x, "}")
    addrow$midrule[[x]] <- sprintf("\\cmidrule(l{2pt}r{2pt}){%i-%i}", cumcount+1, cumcount+count)

    cumcount <- cumcount + count
  }
  addrow$labs1 <- paste(addrow$labs1, collapse = " & ")
  addrow$midrule <- paste(addrow$midrule, collapse = " ")
  addrow$labs2 <- paste(c(keys, cols[[2]]), collapse = "&")
  addrow$rule <- "\\midrule \n"

  addrow$labs1 <- paste0(addrow$labs1, "\\\\\n")
  addrow$labs2 <- paste0(addrow$labs2, "\\\\\n")

  paste(unlist(addrow), collapse = " ")
}



print_df_as_xtab <- function(df, caption, label, filename, add.to.row = NULL, include.colnames = T) {
  if (!"xtable" %in% class(df)){
    xtab <- xtable::xtable(df, caption = caption, label = label, digits = 0)
  }
  print(xtab,
        #sanitize.text.function=identity,
        width = "\\textwidth",
        type = "latex",
        table.placement = "h",
        file = filename,
        append = F,
        NA.string = ".",
        tabular.environment = "tabular",
        #rotate.colnames = ncol(tab) > 10,
        #include.rownames = FALSE,
        #header
        include.rownames = FALSE, #Don't print rownames
        include.colnames = include.colnames, #We create them ourselves
        caption.placement = "top", #"top", NULL
        #hline.after=NULL, #We don't need hline; we use booktabs
        floating=TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
        sanitize.text.function = force, # Important to treat content of first column as latex function
        booktabs = F,
        size = "\\footnotesize",
        html.table.attributes = "border = 0",
        add.to.row = add.to.row)

  filename_preview <- gsub("\\.tex", "_preview.tex", filename)
  file.create(filename_preview)
  cat(c("\\documentclass[]{article}",
        "\\usepackage[margin = 2cm]{geometry}",
        "\\usepackage[utf8]{inputenc}",
        "\\usepackage{tabularx}",
        "\\usepackage{booktabs}",
        "\\usepackage{caption}",
        "\\usepackage{rotating}",
        "\\usepackage{lscape}",
        "\\begin{document}",
        paste0("\\input{", gsub(".*\\/", "", filename), "}"),
        "\\end{document}"),
      sep = "\n",
      file = filename_preview,
      append = F)
}
