#include <Rcpp.h>
using namespace Rcpp;

//' Identify valid adjustment sets
//'
//' Given the set of ancestors A of X and Y and a valid adjustment set
//' \eqn{Z \subseteq A \setminus Desc(X)}, find the subset of Z nearest X (or Y)
//' that is valid for adjustment.
//' @name find_nearest_adjset
//' @details
//' Based on Zander and Liśkiewicz (2020), with a modified Algorithm 3.1 in Koller (2009, p. 75)
//' for finding reachable nodes.
//' Can be used to find:
//' - O-set: the subset closest to Y.
//' - Parent sets: the subset closets to X.
//' - Minimal o-set: the subset of the O-set closest to Y,
//'   i.e. the subset of nodes that can be reached from both X and Y.
//' - Minimal parent set: the subset of parents of X closets to Y,
//'   i.e. the subset of parents that can be reached from both X and Y.
//' @param G a n-by-n adjacency matrix of a backdoor graph. \code{G[i, j] = 1} if there is a path from node i to node j, zero otherwise.
//' @param dmat a  n-by-n adjacency matrix of ancestor relations.
//'   If  \code{dmat[i, j] = 1} `i' is an ancestor of `j' in the original graph.
//'   Used to check if traversed nodes are ancestors of `Z` looking for v-structures.
//' @param X (integer)
//'   column position of source variable (using C++ indexing, from 0 to n-1)
//' @param A (logical vector)
//'   indicates the set of nodes trails must be in.
//'   The search along a path ends if visiting a node i such that `A[i] = FALSE`.
//' @param Z (integer vector)
//'    column positions of conditioning variables (using C++ indexing, from 0 to n-1)
//' @return
//' - `find_nodes_reachable_via_ancestor` implements Algorithm 3.1 in Koller (2009, p. 75)
//' for finding nodes reachable from `X` given `Z` via active trails in `G, with the restriction
//' that only trails where all nodes are in A is considered. On such paths there
//' can not be any colliders, so this implementation includes no check for colliders
//' in Z.
//' - `find_nearest_adjset` return the subset of Z that is reachable from X.
//' @example man/examples/find_nearest_adjset.R
// [[Rcpp::export]]
LogicalVector find_nodes_reachable_via_ancestors(IntegerMatrix G,
                                                 IntegerMatrix dmat,
                                                 int X,
                                                 LogicalVector A,
                                                 IntegerVector Z = 0L) {

 int n = G.ncol();

 // Starting from X, traverse trails with all nodes in A  ----
 IntegerVector LY(1);          // nodes ..
 IntegerVector Ld(1);          // .. and directions to be visited
 LogicalMatrix V(2, n);        // visisted nodes and dir
 LogicalVector R(n);           // reachable nodes

 // init
 LY[0] = X;
 Ld[0] = 1;

 while (LY.length() > 0){
   Rcpp::checkUserInterrupt();

   int Y = LY[0];
   int d = Ld[0];

   LY.erase(0);
   Ld.erase(0);

   if (A[Y]) {      // if Y in set of ancestors A
     if (!V(d, Y)){ // if not already visited

       R[Y] = true;    // mark Y as reachable
       V(d, Y) = true; // track visited nodes

       // if Y is in Z, it blocks the path and no more nodes need to be visited
       // otherwise paths in the direction d is active, and the next nodes on
       // this path has to be checked
       if (is_false(any(Z == Y))) {

         for (int i = 0; i<n; i++){
           // indep of dir, continue downwards by visiting children of Y from top
           if (G(Y, i) == 1 && !V(0, i)){
             LY.push_back(i);
             Ld.push_back(0);

             // continue upwards by visiting parents of Y from bottom
           } else if ((d == 1) && (G(i, Y) == 1) && !V(1, i)) {
             LY.push_back(i);
             Ld.push_back(1);
           }
         }
       } else if (d == 0) {
         // check if Y is a ancestor of any node j in Z
         for (IntegerVector::iterator j = Z.begin(); j != Z.end(); j++) {
           if (dmat(Y, *j) == 1) {

             // continue by visiting parents of Y from bottom
             for (int i = 0; i<n; i++){
               if (G(i, Y) == 1 && !V(1, i)) {
                 LY.push_back(i);
                 Ld.push_back(1);
               }
             }
             break;
           }
         }
       }
     }
   }
 }

 return R;
}

//' @rdname find_nearest_adjset
//' @export
// [[Rcpp::export]]
IntegerVector find_nearest_adjset(IntegerMatrix G, IntegerMatrix dmat, int X, LogicalVector A, IntegerVector Z) {
   if (Z.length() == 0 || na_omit(Z).length() == 0) {
      return Z;
   } else {
      LogicalVector R = find_nodes_reachable_via_ancestors(G, dmat, X, A, Z);
      LogicalVector indx = R[Z];
      return Z[indx];
   }
}


