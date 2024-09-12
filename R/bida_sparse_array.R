


#' class: cida sparse array
#'
#' Sparse representation of multi-dimensional arrays.
#'
#' @name bida_sparse_array
#' @param x (numeric vector) values of non-zero elements
#' @param index (integer vector) position of non-zero elements
#' @param dims (integer vector) dimension of array
#'
#' @return an object of class `bida_sparse_array`
#' @keywords internal
#'
#' @examples
#'
#' vals <- seq(3, 27, by = 3)
#' x <- new_bida_sparse_array(vals, vals-1, c(3, 3, 3))
#' x
#' as.array(x)
#'
#' # addition
#' x+1
#' y <- new_bida_sparse_array(0:2, 0:2, c(3, 3, 3), default = 1)
#' x+y
#' stopifnot(all(as.array(x+y) == as.array(x) + as.array(y)))
#'
#' # rep
#'
#' dims <- 2:4
#' x <- new_bida_sparse_array(rep(1, prod(dims)), seq_len(prod(dims))-1, dims)
#' y <- aperm(x, c(3:1))
#' stopifnot(all.equal(as.array(y), aperm(as.array(x), c(3:1))))
#'
#'
#' # compare size of sparse and non-sparse arrays
#' obj_size <- matrix(NA, 24, 2)
#' sparsity <- .25 # share of non-zero elements
#' for (i in seq_len(nrow(obj_size))) {
#'   dims <- rep(2, i)
#'   index  <- sample(prod(dims), round(prod(dims)*sparsity))
#'   counts <- sample(100, size = length(index), replace = T)
#'
#'   sparse <- new_bida_sparse_array(counts, index, dims)
#'   arr    <- as.array(sparse)
#'   obj_size[i, 1] <- pryr::object_size(sparse)
#'   obj_size[i, 2] <- pryr::object_size(arr)
#' }
#'
#' matplot(2**seq_len(nrow(obj_size)), obj_size, type = "l")
#' plot(2**seq_len(nrow(obj_size)), obj_size[, 1]/obj_size[, 2], type = "l")
#'


#' @rdname bida_sparse_array
#' @export
bida_sparse_array <- function(value, index, dim, dimnames = NULL, default = 0) {
  tmp <- sort.int(index, index.return = TRUE)
  new_bida_sparse_array(value[tmp$ix], tmp$x, dim, dimnames, default)
}

new_bida_sparse_array <- function(value, index, dim, dimnames = NULL, default = 0) {
  structure(list(value = value,
                 index = index,
                 dim = dim,
                 dimnames = dimnames,
                 default = default),
            class = "bida_sparse_array")
}





# Generics ----


#' @rdname bida_sparse_array
#' @export
as.array.bida_sparse_array <- function(x) {
  y <- array(x$default, x$dim, x$dimnames)
  y[x$index+1] <- x$value
  return(y)
}


#' @rdname bida_sparse_array
#' @export
rep.bida_sparse_array <- function(x, times = 1, each = 1) {
  if (times == 1 && each == 1) return(x)
  len <- prod(x$dim)

  rep_index   <- rep(x$index, each = each, times = times)
  each_index  <- rep(seq_len(each)-1, times = length(x$index)*times)
  times_index <- rep(each*len*(seq_len(times)-1), each = length(x$index)*each)

  new_index   <- each_index + each*rep_index + times_index
  new_value   <- rep(x$value, each = each, times = times)
  new_dim     <- times*each*prod(x$dim)
  new_bida_sparse_array(new_value, new_index, dim = new_dim, default = x$default)
}


#' @rdname bida_sparse_array
#' @export
aperm.bida_sparse_array <- function(x, perm) {
  indx <- x$index
  dims <- x$dim
  n <- length(dims)

  stride <- c(1, cumprod(dims[-n]))
  conf <- mapply(function(s, k) (indx%/%s)%%k, stride, dims)

  new_dims <- dims[perm]
  new_indx <- c(conf[, perm, drop = F]%*%c(1, cumprod(new_dims[-n])))

  x$dim    <- new_dims
  x$index  <- new_indx

  return(x)
}


## Arithmetics ----
#' @rdname bida_sparse_array
#' @export
`+.bida_sparse_array` <- function(x, y){
  aritmethics_bida_sparse_array("+", x, y)
}

`*.bida_sparse_array` <- function(x, y){
  aritmethics_bida_sparse_array("*", x, y)
}

#' @rdname bida_sparse_array
#' @export
`-.bida_sparse_array` <- function(x, y){
  aritmethics_bida_sparse_array("-", x, y)
}
#' @rdname bida_sparse_array
#' @export
`/.bida_sparse_array` <- function(x, y) {
  aritmethics_bida_sparse_array("/", x, y)
}

aritmethics_bida_sparse_array <- function(fun, x, y) {
  stopifnot(inherits(x, "bida_sparse_array"))

  f <- match.fun(fun)
  if (is.numeric(y) && length(y) == 1) {
    stopifnot(length(y) == 1)
    x$value <- f(x$value, y)
    x$default <-  f(x$default, y)
    return(x)
  } else if (inherits(y, "bida_sparse_array")) {

    stopifnot((prod(x$dim) == prod(y$dim)))

    index_yinx  <- match(y$index, x$index, 0L)
    index_xiny  <- match(x$index, y$index, 0L)

    value  <- c(f(x$value[index_xiny == 0], y$default),
                f(y$value[index_yinx == 0], x$default),
                f(x$value[index_yinx], y$value[index_xiny]))
    default <- f(x$default, y$default)
    index  <- c(x$index[index_xiny == 0],
                y$index[index_yinx == 0],
                x$index[index_yinx])
    tmp <- sort.int(index, index.return = TRUE)
    new_bida_sparse_array(value[tmp$ix], tmp$x, x$dim, x$dimnames, default)

  } else {
    stop("y must be a numeric constant or a bida_sparse_array")
  }
}


colIndex <- function(x, dims = 1){
  x$index%/%prod(x$dim[seq_len(dims)])
}
rowIndex <- function(x, dims = 1) {
  x$index%%prod(x$dim[seq_len(dims)])
}


colSums_bida_sparse_array <- function(x, dims = 1) {
  tmp <- seq_len(dims)
  colIndex  <- x$index%/%prod(x$dim[tmp])
  newIndex <- unique(group)
  new_bida_sparse_array(rowsum_fast(x$value, colIndex, newIndex),
                        newIndex,
                        x$dim[-tmp])
}

#' @export
levels.bida_sparse_array <- function(x) {
  if (is.null(x$dimnames)) {
    setNames(lapply(x$dim-1, seq.int, from = 0), paste0("X", seq_along(x$dim)))
  } else {
    x$dimnames
  }
}


get_stride <- function(x, MARGIN = seq_along(x$dim)) {
  c(1, cumprod(x$dim[MARGIN[-length(MARGIN)]]))
}

get_coordinates <- function(x, MARGIN, stride = get_stride(x)) {
  vapply(MARGIN,
         function(i) x$index%/%stride[i]%%x$dim[i],
         numeric(length(x$index)))
}
get_index <- function(x, MARGIN) {
  coord <- get_coordinates(x, MARGIN, get_stride(x))
  coord%*%get_stride(x, MARGIN)
}


marginalize <- function(x, MARGIN) {

  new_dim <- x$dim[MARGIN]
  new_dimnames <- x$dimnames[MARGIN]

  # compute non-missing indicies of MARGIN
  index <- get_index(x, MARGIN)
  new_index <- c(unique(index))

  if (x$default == 0) {
    new_value <- rowsum_fast(x$value, index, new_index)
    new_default <- 0
  } else {
    # count missing elements for each observed level of MARGIN
    # - add to grouped sums below
    k     <- prod(x$dim)/prod(new_dim)
    nmiss <- k - tabulate(match(index, new_index))

    # compute sums over MARGIN
    new_value <- rowsum_fast(x$value, index, new_index) + nmiss*x$default

    # compute new default value
    new_default <- k*x$default
  }

  # return sparse array
  new_bida_sparse_array(new_value, new_index, new_dim, new_dimnames, new_default)
}


backdoor_mean <- function(counts, x, y, z, ess) {

  nxz <- marginalize(counts, c(x, z))
  nz  <- marginalize(nxz, length(x)+seq_along(z))


}
