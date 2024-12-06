#' Title
#'
#' @param df
#' @param values_from
#' @param names_from
#' @param file
#' @param caption
#' @param label
#' @param digits
#' @param group_keys_in_rows
#' @param footnote
#' @param addTimeStamp
#' @param bold_fun
#' @param sanitize.text.function
#'
#' @return
#' @export
#'
#' @examples
#'
#' df <- data.frame(group = rep(letters[1:3], each = 6),
#'                  key   = rep(rep(1:3, each = 2), 3),
#'                  name  = rep(rep(LETTERS[1:2], 3), 3),
#'                  value1 = 1:18,
#'                  value2 = 18:1)
#'
#' values_from <- c("value1", "value2")
#' names_from  <- "name"
#' group_by    <- "group"
#' df_to_tex(df, values_from, names_from, group_by)
#'
#' # put group keys in separating row
#' df_to_tex(df, values_from, names_from, group_by, group_keys_in_rows = TRUE)
#'
#' # make maximum of each row bold
#' df_to_tex(df, values_from, names_from, group_by, group_keys_in_rows = TRUE, bold_fun = "max")
#'
df_to_tex <- function(df,
                      values_from = "value",
                      names_from = character(0),
                      group_by = character(0),
                      file = "",
                      caption = NULL,
                      label = NULL,
                      digits = 2,
                      group_keys_in_rows = FALSE,
                      footnote = "",
                      addTimeStamp = FALSE,
                      bold_fun = NULL,
                      sanitize.text.function = NULL)
{

  if (nrow(df) == 0) {
    warning("Empty df, table", file, "is not generated.")
    return(df)
  }

  # init
  align <- NULL
  add_to_row <- list(pos = list(),
                     command = c())

  keys  <- names(df)[!names(df) %in% c(values_from, names_from)]
  if (length(group_by) > 0) {
    df <- dplyr::group_by(df, across(dplyr::all_of(group_by)))
    if (group_keys_in_rows) {
      keys <- keys[!keys %in% group_by]
    }
  }

  # from long to wide ---
  # make multi-header column if more than one value column and a names_from
  if (length(names_from) > 1) {
    stop("names_from can only refer to a single column.")
  } else if (length(names_from) == 0) {
    include.colnames <- TRUE
  } else if (length(values_from) == 1 && length(names_from) == 1) {
    include.colnames <- TRUE
    df <- tidyr::pivot_wider(df, values_from = values_from, names_from = names_from, names_expand = FALSE)
  } else {  # make multi-column header
    include.colnames <- FALSE
    df <- tidyr::pivot_wider(df,
                             values_from = values_from,
                             names_from = dplyr::all_of(names_from),
                             names_expand = FALSE,
                             names_sep = "..")


    # remove columns with all NA
    indx <- vapply(df, function(x) all(is.na(x)), logical(1))
    df   <- df[, !indx]

    colnames  <- names(df)[!names(df) %in% c(keys, group_by)]
    colnames_top <- sub("\\.\\..*$", "", colnames)
    colnames_bottom  <- sub("^.*\\.\\.", "", colnames)

    pos <- list(0)
    tmp <- mapply(function(x, n) sprintf("\\multicolumn{%s}{c}{%s}", n, x),
                  unique(colnames_top), table(colnames_top)[unique(colnames_top)])
    command1 <- paste0(c(rep(" ", length(keys)), tmp), collapse = "&")
    command2 <- sprintf("\\cline{%s-%s}", length(keys)+1, length(keys)+length(colnames))
    command3 <- paste0(c(keys, colnames_bottom), collapse = " & ")
    command  <- sprintf("%s   \\\\ \n %s  \\\\ \n %s  \\\\ \n",
                        command1, command2, command3)

    add_to_row$pos <- c(add_to_row$pos, pos)
    add_to_row$command <- c(add_to_row$command, command)

    align <- paste0(c(paste0(rep("c", length(keys)+1), collapse = ""),  # add one for rownames - even when not printed
                      sapply(table(colnames_top)[unique(colnames_top)], function(n) paste0(rep("c", n), collapse = ""))),
                    collapse = "|")
  }

  if (!is.null(bold_fun)) {
    # bold extremes in each row -----
    tmp <- round(df[, !names(df) %in% c(group_by, keys)], digits)
    pos_extreme <- which(tmp == apply(tmp, 1, match.fun(bold_fun)))
    tmp <- apply(tmp, 2, as.character)
    tmp[pos_extreme] <- paste0("\\textbf{", tmp[pos_extreme], "}")
    df <- data.frame(df[, unique(c(group_by, keys))], tmp)
    sanitize.text.function <- identity
  }

  if (length(group_by) > 0) {
    # add empty row between groups -----

    gr_keys    <- dplyr::group_keys(df)
    gr_indices <- dplyr::group_indices(df)

    if (group_keys_in_rows) {
      # ungroup and remove key-columns
      df  <- dplyr::select(dplyr::ungroup(df), -all_of(group_by))

      # write multicolumn-lines separating each group
      tmp <- apply(gr_keys, 1, paste, collapse = ". ")
      command <- sapply(tmp,
                        function(x) sprintf("\\\\ \n %s \\multicolumn{%s}{c}{%s} \\\\ \n ",
                                            paste(rep("&", length(keys), collapse = "")),
                                            ncol(df)-length(keys),
                                            x))
      command[1] <- paste0("\n\\midrule", command[1])
      pos  <- as.list(c(0, cumsum(tabulate(gr_indices)[-max(gr_indices)])))
    } else {
      pos  <- as.list(cumsum(tabulate(gr_indices)[-max(gr_indices)]))
      command <- rep("\\\\ \n", length(pos))
    }
    add_to_row$pos <- c(add_to_row$pos, pos)
    add_to_row$command <- c(add_to_row$command, command)
  }

  if (length(footnote) > 0 || addTimeStamp) {
    # add footnote -----
    if (addTimeStamp) {
      footnote <- paste0(footnote,
                         "Table printed at ", format(Sys.time(), "%d.%m.%y %H:%M"))
    }
    command <- sprintf("\\bottomrule \n \\multicolumn{%s}{l}{\\scriptsize %s} \n ", ncol(df), footnote)
    pos <- nrow(df)
    add_to_row$pos <- c(add_to_row$pos, pos)
    add_to_row$command <- c(add_to_row$command, command)
  }

  # debug
  if (!is.null(align) && nchar(gsub("\\|", "", align)) != ncol(df) +1) {
    print(list(colnames = colnames,
               values_from = values_from,
               align = align,
               head(df)))
    stop()
  }

  hline_after <- c(-1,
                   switch(group_keys_in_rows+1, 0, integer(0)),
                   switch((length(footnote) > 0 || addTimeStamp) +1, nrow(xtab), integer(0)))

  # replace underscore  "_" with "."
  add_to_row$command <- gsub("_", ".", add_to_row$command)
  colnames(df) <- gsub("_", ".", colnames(df))


  cat("Save file", file, "\n")
  print(names(df))
  if (file.exists(file)) file.remove(file)
  print(#xtab,
    xtable::xtable(df, align = align, caption = caption, label = label, digits = digits),
    file = file,
    include.rownames = FALSE,
    booktabs = TRUE,
    add.to.row = add_to_row,
    caption.placement = "top",
    hline.after = hline_after,
    include.colnames = include.colnames,
    table.placement = "ht",
    sanitize.text.function = sanitize.text.function)

}
