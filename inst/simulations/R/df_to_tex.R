df_to_tex <- function(df,
                      values_from = "value",
                      names_from = NULL,
                      file = "",
                      caption = NULL,
                      label = NULL,
                      digits = 2,
                      group_keys_in_rows = FALSE,
                      note = "",
                      addTimeStamp = FALSE)
{


  keys  <- names(df)[!names(df) %in% c(values_from, names_from)]
  gr_vars <- group_vars(df)

  if (length(gr_vars) > 0 && group_keys_in_rows) {
    keys <- keys[!keys %in% gr_vars]
  }
  align <- NULL
  add_to_row <- list(pos = list(),
                     command = c())

  if (length(names_from) == 0 || length(values_from) == 1) {
    df <- tidyr::pivot_wider(df, values_from = values_from, names_from = names_from)
  } else {

    colnames <- levels(interaction(df[names_from]))
    pos <- list(0)

    tmp <- sapply(colnames,
                  function(x) sprintf("\\multicolumn{%s}{c}{%s}", length(values_from), x))
    command1 <- paste0(c(rep(" ", length(keys)), tmp), collapse = "&")
    command2 <- sprintf("\\cline{%s-%s}", length(keys)+1, length(keys)+length(values_from)*length(colnames))
    command3 <- paste0(c(keys, rep(values_from, length(colnames))), collapse = " & ")
    command  <- sprintf("%s   \\\\ \n %s  \\\\ \n %s  \\\\ \n",
                        command1, command2, command3)

    add_to_row$pos <- c(add_to_row$pos, pos)
    add_to_row$command <- c(add_to_row$command, command)


    df <- df %>%
      tidyr::pivot_wider(values_from = all_of(values_from),
                         names_from = all_of(names_from))

    order_cols <- c(keys, sapply(colnames, function(x) names(df)[endsWith(names(df), x)]))
    df <- df %>%
      select(all_of(c(gr_vars, keys, order_cols))) %>%
      arrange(across(all_of(c(gr_vars, keys))), .by_group = TRUE)


    align <- paste0(c(paste0(rep("r", length(keys)+1), collapse = ""),  # add one for rownames - even when not printed
                      rep(paste0(rep("c", length(values_from)), collapse = ""), length(colnames))),
                      collapse = "|")
  }

  if (length(gr_vars) > 0) {
    gr_keys <- group_keys(df)
    gr_indices <- group_indices(df)

    if (group_keys_in_rows) {
      # ungroup and remove key-columns
      df  <- df %>%
        ungroup() %>%
        select(-all_of(gr_vars))

      # write multicolumn-lines separating each group
      tmp <- apply(gr_keys, 1, paste, collapse = ". ")
      command <- sapply(tmp,
                        function(x) sprintf("\\\\ \n %s \\multicolumn{%s}{c}{%s} \\\\ \n ",
                                            paste(rep("&", length(keys), collapse = "")),
                                            ncol(df)-length(keys),
                                            x))
      command[1] <- paste0("\n\\midrule", command[1])
      pos  <- as.list(c(0, cumsum(tabulate(gr_indices)[-1])))
    } else {
      pos  <- as.list(cumsum(tabulate(gr_indices)[-1]))
      command <- rep("\\\\ \n", length(pos))
    }
    add_to_row$pos <- c(add_to_row$pos, pos)
    add_to_row$command <- c(add_to_row$command, command)
  }

  if (length(note) > 0 || addTimeStamp) {

    if (addTimeStamp) {
      note <- paste0(note,
                    "Table printed at ", format(Sys.time(), "%d.%m.%y %H:%M"))
    }
    command <- sprintf("\\bottomrule \n \\multicolumn{%s}{l}{\\scriptsize %s} \n ", ncol(df), note)
    pos <- nrow(df)
    add_to_row$pos <- c(add_to_row$pos, pos)
    add_to_row$command <- c(add_to_row$command, command)
  }

  xtab <- xtable::xtable(df, align = align, caption = caption, label = label, digits = digits)
  if (file.exists(file)) file.remove(file)
  hline_after <- c(-1,
                   switch(group_keys_in_rows+1, 0, integer(0)),
                   switch((length(note) > 0 || addTimeStamp) +1, nrow(xtab), integer(0)))
  print(xtab,
        file = file,
        include.rownames = FALSE,
        booktabs = TRUE,
        add.to.row = add_to_row,
        caption.placement = "top",
        hline.after = hline_after,
        include.colnames = is.null(names_from) || length(values_from) == 1,
        table.placement = "ht")
}
