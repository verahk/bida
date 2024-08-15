

partition_from_cpt_array <- function(cpt) {

    if (length(dim(cpt)) < 2) return(NULL)

    r <- dim(cpt)[1]
    q <- prod(dim(cpt))/r

    # enumerate unique probs for each outcome
    tmp <- apply(matrix(cpt, q, r, byrow = T), 2,
                 function(x) match(x, unique(x)))

    # enumerate distributions
    parts <- tmp%*%nrow(tmp)**seq.int(0, r-1)
    # split(data.frame(tmp), parts) check

    # return list with  members of each part
    unname(split(seq_len(q)-1, parts))
}
