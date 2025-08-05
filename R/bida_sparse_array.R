


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
#' arr <- as.array(x)
#'
#' # addition
#' x+1
#' y <- new_bida_sparse_array(0:2, 0:2, c(3, 3, 3), default = 1)
#' x+y
#' stopifnot(all(as.array(x+y) == as.array(x) + as.array(y)))
#'
#' colSums(x)
#' # rep
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
new_bida_sparse_array <- function(value, index, dim, dimnames, default) {
  structure(list(value = value,
                 index = index,
                 dim = dim,
                 dimnames = dimnames,
                 default = default),
            class = "bida_sparse_array")
}

#' @rdname bida_sparse_array
#' @export
bida_sparse_array <- function(value, ..., default = 0) {
  UseMethod("bida_sparse_array")
}

#' @rdname bida_sparse_array
#' @export
bida_sparse_array.array <- function(value, default = 0) {
  index = seq_along(value)[!value == default]
  new_bida_sparse_array(value = value[index],
                        index = index-1,
                        dim = dim(value),
                        dimnames = dimnames(value),
                        default = default)
}

#' @rdname bida_sparse_array
#' @export
bida_sparse_array.default <- function(value, index, dim, dimnames = NULL, default = 0) {
  stopifnot(length(index) == length(value))
  stopifnot(all(index%%1 == 0) && min(index) >= 0 && max(index) < prod(dim))
  stopifnot(is.null(dimnames) || all(lengths(dimnames) == 0) || all(lengths(dimnames) == dim))
  new_bida_sparse_array(value, index, dim, dimnames, default)
}


#' @rdname bida_sparse_array
#' @export
as.array.bida_sparse_array <- function(x) {
  y <- array(x$default, x$dim, x$dimnames)
  y[x$index+1] <- x$value
  y
}


# Methods ----
#' @rdname bida_sparse_array
#' @export
dim.bida_sparse_array <- function(x) {
  x$dim
}

#' @rdname bida_sparse_array
#' @export
`dim<-.bida_sparse_array` <- function(x, value) {
  stopifnot(prod(value) > max(x$index))
  x$dim <- value
  x
}


#' @rdname bida_sparse_array
#' @export
dimnames.bida_sparse_array <- function(x) {
  x$dimnames
}
#' @rdname bida_sparse_array
#' @export
`dimnames<-.bida_sparse_array` <- function(x, value) {
  stopifnot(all(lengths(value) == dim(x)) || all(lengths(value) == 0))
  x$dimnames <- value
  x
}

#' @rdname bida_sparse_array
#' @export
rep.bida_sparse_array <- function(x, times = 1, each = 1) {
  if (times == 1 && each == 1) return(x)

  len <- prod(x$dim)
  rep_index   <- rep(x$index, each = each, times = times)
  each_index  <- rep(seq_len(each)-1, times = length(x$index)*times)
  times_index <- rep(each*len*(seq_len(times)-1), each = length(x$index)*each)

  new_value   <- rep(x$value, each = each, times = times)
  new_index   <- each_index + each*rep_index + times_index
  new_dim     <- times*each*prod(x$dim)
  new_bida_sparse_array(new_value, new_index, new_dim, NULL, x$default)
}


#' @rdname bida_sparse_array
#' @export
aperm.bida_sparse_array <- function(x, perm) {

  dims <- dim(x)
  n <- length(dims)

  stopifnot(length(perm) == n)
  if (n == 1) return(x)

  # map to new index
  coord <- arrayInd(x$index+1, dims)-1
  stride <- c(1, cumprod(dims[perm[-n]]))
  x$index <- coord[, perm]%*%stride

  dim(x) <- dims[perm]
  dimnames(x) <- dimnames(x)[perm]
  return(x)
}


## Arithmetics ----
#' @rdname bida_sparse_array
#' @export
sum.bida_sparse_array <- function(x, ...) {
  if (x$default == 0) {
    sum(x$value)
  } else {
    sum(x$value) + (prod(dim(x))-length(x$index))*x$default
  }
}
#' @rdname bida_sparse_array
#' @export
`+.bida_sparse_array` <- function(x, y){
  aritmethics_bida_sparse_array("+", x, y)
}


#' @rdname bida_sparse_array
#' @export
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
    x$value   <- f(x$value, y)
    x$default <- f(x$default, y)
    return(x)
  } else if (inherits(y, "bida_sparse_array")) {

    stopifnot((prod(x$dim) == prod(y$dim)))

    yinx  <- match(y$index, x$index, 0L)
    if (all(yinx == 0)) {
      value <- c(f(x$value, y$default), f(x$default, y$value))
      index <- c(x$index, y$index)
      default <- f(x$default, y$default)
    } else {
      x$index[yinx]
      x$index[-yinx]
      y$index[yinx == 0]
      y$index[yinx > 0]

      value  <- c(f(x$value[-yinx], y$default),
                  f(x$default, y$value[yinx == 0]),
                  f(x$value[yinx], y$value[match(x$index[yinx], y$index, 0L)]))
      index  <- c(x$index[-yinx],
                  y$index[yinx == 0],
                  x$index[yinx])
      default <- f(x$default, y$default)
    }

    new_bida_sparse_array(value, index, x$dim, x$dimnames, default)

  } else {
    stop("y must be a numeric constant or a bida_sparse_array")
  }
}

#' @rdname bida_sparse_array
#' @export
rowSums.bida_sparse_array <- function(x, na.rm = FALSE, dims = 1L){

  dimx <- dim(x)
  keep <- seq.int(1, dims)  # margins to be kept
  new_dim   <- dim(x)[keep]
  new_dimnames <- dimnames(x)[keep]

  # compute row-indicies
  index <- x$index%%prod(dimx[keep])
  new_index <- unique(index)

  # aggregate values
  if (x$default == 0) {
    new_value  <- rowsum_fast(x$value, index, new_index)
    new_default <- 0
  } else {
    k <- prod(dimx[-keep])    # cardinality of margins summed out
    nmiss <- k-tabulate(match(index, new_index))
    new_value <- rowsum_fast(x$value, index, new_index) + nmiss*x$default
    new_default <- k*x$default
  }
  new_bida_sparse_array(new_value, new_index, new_dim, new_dimnames, new_default)
}

#' @rdname bida_sparse_array
#' @export
colSums.bida_sparse_array <- function(x, na.rm = FALSE, dims = 1L){
  dimx <- dim(x)
  keep <- seq.int(dims+1, length(dimx))

  new_dim   <- dim(x)[keep]
  new_dimnames <- dimnames(x)[keep]

  index <- x$index%/%prod(dimx[-keep])
  new_index <- unique(index)

  # aggregate values
  if (x$default == 0) {
    new_value  <- rowsum_fast(x$value, index, new_index)
    new_default <- 0
  } else {
    k <- prod(dimx[-keep])    # cardnality of margins summed out
    nmiss <- k-tabulate(match(index, new_index))
    new_value <- rowsum_fast(x$value, index, new_index) + nmiss*x$default
    new_default <- k*x$default
  }
  new_bida_sparse_array(new_value, new_index, new_dim, new_dimnames, new_default)
}
