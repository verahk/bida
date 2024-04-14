#include <Rcpp.h>
using namespace Rcpp;


//' @rdname factor_product_c
// [[Rcpp::export]]
CharacterVector combine_char_vectors_c(CharacterVector x, CharacterVector y){
  int nx = x.size();
  int ny = y.size();
  int n  = nx + ny;
  int k =  nx;

  CharacterVector out(n);

  for (int i = 0; i < nx; i++) {
    out[i] = x[i];
  }

  for (int j = 0; j < ny; j++) {
    out[k] = y[j];
    k++;
  }

  return out;
}

//' @rdname factor_product_c
// [[Rcpp::export]]
NumericVector named_array_c(NumericVector values, NumericVector dim, CharacterVector names){
  List dimnames(names.size());
  dimnames.attr("names") = names;

  values.attr("dim") = dim;
  values.attr("dimnames") = dimnames;
  return values;
}

//' @rdname factor_product_c
// [[Rcpp::export]]
CharacterVector scope_c(NumericVector x) {
  List tmp = x.attr("dimnames");
  return tmp.names();
}

// [[Rcpp::export]]
NumericVector stride_c(NumericVector x) {
  NumericVector nlev = x.attr("dim");
  int n = nlev.length();

  NumericVector out(n);
  out[0] = 1;
  int cump  = 1;
  for(int i = 1; i < n; i++) {
    cump *= nlev[i-1];
    out[i] = cump;
  }
  return out;
}


//' Compute product of two factors
//'
//' @param x,y numeric vectors with:
//' - dimensions/cardinalities given by `x.attr("dim")` and similarly for `y`
//' - scope given by names of `x.attr("dimnames")` and similarly for `y`
//' @seealso cpqueries_from_cpt_arrays
//' @keywords internal
//'
// [[Rcpp::export]]
NumericVector factor_product_c(NumericVector x, NumericVector y) {

  // scope ----
  CharacterVector scope_x = scope_c(x);
  CharacterVector scope_y = scope_c(y);

  int nx = scope_x.length();
  int ny = scope_y.length();

  // joint scope
  CharacterVector scope = combine_char_vectors_c(scope_x, scope_y);
  scope = scope[duplicated(scope) == 0];
  int n = scope.size();

  // cardinality  and stride
  NumericVector nlev(n);
  NumericVector nlev_x = x.attr("dim");
  NumericVector nlev_y = y.attr("dim");

  NumericVector stride_x(n), stride_y(n);
  NumericVector stride_xx = stride_c(x);
  NumericVector stride_yy = stride_c(y);


  // set stride of x and y to zero for variables not in respective scope

  for (int i = 0; i < nx; i++) {
    nlev[i] = nlev_x[i];
    stride_x[i] = stride_xx[i];
  }

  for (int i = 0; i < ny; i++) {
    for (int j = 0; j < n; j++) {
      if (scope[j] == scope_y[i]) {
        nlev[j] = nlev_y[i];
        stride_y[j] = stride_yy[i];
        break;
      }
    }
  }

  int len = 1;
  for(NumericVector::iterator i = nlev.begin(); i != nlev.end(); ++i) {
    len *= *i;
  }

  if (log(len) > 25) {
    stop("Vector exceeds 2**25 elements");
  }

  // compute factor product
  NumericVector out(len);
  NumericVector ass(n);
  int j = 0;
  int k = 0;
  for (int i = 0; i < len; i++){

    if (i % 1000 == 0)
    Rcpp::checkUserInterrupt();

    out[i] = x[j]*y[k];

    for (int v = 0; v < n; v++){
      if (ass[v] == nlev[v]-1){
        ass[v] = 0;
        j = j - (nlev[v]-1)*stride_x[v];
        k = k - (nlev[v]-1)*stride_y[v];
      } else {
        ass[v] = ass[v]+1;
        j = j + stride_x[v];
        k = k + stride_y[v];
        break;
      }
    }
  }
  return named_array_c(out, nlev, scope);
}

//' @rdname factor_product_c
// [[Rcpp::export]]
NumericVector factors_product_c(List factors) {
  NumericVector f = factors[0];
  for (int i = 1; i < factors.length(); i++)
    f = factor_product_c(f, factors[i]);
  return f;
}

/*** R

new_factor <- function(val, nlev, scope) array(val, nlev, setNames(vector("list", length(nlev)), scope))
x <- new_factor(1:4, c(2, 2), c("a", "b"))
y <- new_factor(1:4, c(2, 2), c("c", "d"))
f <- factor_product_c(x, y)
get_scope(f)
f

all.equal(factor_product_c(x, y),
          factors_product_c(list(x, y)))

k <- rep(2, 3)
x <- new_factor(seq_len(prod(k)), k, letters[seq_along(k)])
y <- new_factor(seq_len(prod(k)), k, LETTERS[seq_along(k)])
microbenchmark::microbenchmark(
  factor_product(x, y),
  factor_product_c(x, y)
)

k <- rep(2, 5)
x <- new_factor(seq_len(prod(k)), k, letters[seq_along(k)])
y <- new_factor(seq_len(prod(k)), k, letters[3+seq_along(k)])
microbenchmark::microbenchmark(
  factor_product(x, y),
  factor_product_c(x, y)
)


*/

