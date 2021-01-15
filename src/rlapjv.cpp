#include <Rcpp.h>
#include "lapjv.h"
using namespace Rcpp;




// int_t lapjv_internal(const uint_t n, cost_t *cost[],
//                                int_t *x, int_t *y);

// [[Rcpp::export]]
IntegerVector cpp_lapjv(NumericMatrix cost, bool maximize = false) {
  const uint_t nc = cost.ncol(), nr = cost.nrow();
  IntegerVector x(nc);
  IntegerVector y(nr);
  double max_cost = max(cost);
  if(maximize)
  {
    cost = max_cost - cost;
  }
  // Convert cost to *cost_t[], an array of pointers
  cost_t **cost_p;

  NEW(cost_p, cost_t*, nc);
  for (int i = 0; i < nr; ++i)
  {
    cost_p[i] = cost.column(i).begin();
  }
  int c = lapjv_internal(nr, cost_p, x.begin(), y.begin());
  if(maximize)
  {
    cost = max_cost - cost;
  }
  return y + 1;
}

// int lapmod_internal(
//     const uint_t n, cost_t *cc, uint_t *ii, uint_t *kk,
//     int_t *x, int_t *y, fp_t fp_version)


// [[Rcpp::export]]
IntegerVector cpp_lapmod(int n, NumericVector cc, IntegerVector ii, IntegerVector kk, bool maximize = false) {
  int nr = n;
  // assert(kk.length() == cc.length())
  // assert(ii.length() == nr+1)
  // assert(max(kk) < n)
  IntegerVector x(nr);
  IntegerVector y(nr);
  fp_t fp_version = FP_1;

  double max_cost = max(abs(cc));
  if( maximize )
  {
    cc = max_cost - cc;
  }

  // std::cout << cc << "\n";
  // std::cout << ii << "\n";
  // std::cout << kk << "\n";

  int c = lapmod_internal(nr, cc.begin(), ii.begin(), kk.begin(),
    x.begin(), y.begin(), fp_version, max_cost);
  // std::cout << "Done\n";
  // std::cout << c << "\n";
  // std::cout << x << "\n";
  // std::cout << y << "\n";
  if( maximize )
  {
    cc = max_cost - cc;
  }
  return y + 1;
}