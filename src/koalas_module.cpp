#include <Rcpp.h>

#include "KoalaGroup.h"

struct Settings
{
  bool debug = true;
  double tol = 1e-7;
};
constexpr Settings debug { .debug = true };
// constexpr Settings run { .debug = false };

/* TODO: link to blofeld for this! */
template <class RcppModuleClassName>
RcppModuleClassName* invalidate_default_constructor() {
  Rcpp::stop("Default constructor is disabled for this class");
  return 0;
} //blofeld


using KoalaGroupD3 = KoalaGroup<debug, 3U, 3U, 3U, 3U, 3U>;
RCPP_EXPOSED_AS(KoalaGroupD3)
RCPP_EXPOSED_WRAP(KoalaGroupD3)

using KoalaGroupD1 = KoalaGroup<debug, 1U, 1U, 1U, 1U, 1U>;
RCPP_EXPOSED_AS(KoalaGroupD1)
RCPP_EXPOSED_WRAP(KoalaGroupD1)

using KoalaGroupD0 = KoalaGroup<debug, 0U, 0U, 0U, 0U, 0U>;
RCPP_EXPOSED_AS(KoalaGroupD0)
RCPP_EXPOSED_WRAP(KoalaGroupD0)


RCPP_MODULE(koalas){

	using namespace Rcpp;

  class_<KoalaGroupD3>("KoalaGroupD3")
    .factory(invalidate_default_constructor)
    .constructor<IntegerVector const, NumericVector const, NumericVector const>("Constructor with ncomps vector, parameter vector and state vector")
    .constructor<KoalaGroupD3&>("Copy constructor")
    .method("clone", &KoalaGroupD3::clone, "Clone method (probably better to use the copy c'tor)")

    .method("update", &KoalaGroupD3::update, "Update method")
    .method("active_intervention", &KoalaGroupD3::active_intervention, "Apply an active sampling intervention")

    .property("vitals", &KoalaGroupD3::get_vitals, "Get vitals (yeay/day/alive)")
    .property("state", &KoalaGroupD3::get_state, &KoalaGroupD3::set_state, "Get and set state")
    .property("parameters", &KoalaGroupD3::get_pars, &KoalaGroupD3::set_pars, "Get and set parameters")
  ;

  class_<KoalaGroupD1>("KoalaGroupD1")
    .factory(invalidate_default_constructor)
    .constructor<IntegerVector const, NumericVector const, NumericVector const>("Constructor with ncomps vector, parameter vector and state vector")
    .constructor<KoalaGroupD1&>("Copy constructor")
    .method("clone", &KoalaGroupD1::clone, "Clone method (probably better to use the copy c'tor)")

    .method("update", &KoalaGroupD1::update, "Update method")
    .method("active_intervention", &KoalaGroupD1::active_intervention, "Apply an active sampling intervention")

    .property("vitals", &KoalaGroupD1::get_vitals, "Get vitals (yeay/day/alive)")
    .property("state", &KoalaGroupD1::get_state, &KoalaGroupD1::set_state, "Get and set state")
    .property("parameters", &KoalaGroupD1::get_pars, &KoalaGroupD1::set_pars, "Get and set parameters")
  ;

  class_<KoalaGroupD0>("KoalaGroupD0")
    .factory(invalidate_default_constructor)
    .constructor<IntegerVector const, NumericVector const, NumericVector const>("Constructor with ncomps vector, parameter vector and state vector")
    .constructor<KoalaGroupD0&>("Copy constructor")
    .method("clone", &KoalaGroupD0::clone, "Clone method (probably better to use the copy c'tor)")

    .method("update", &KoalaGroupD0::update, "Update method")
    .method("active_intervention", &KoalaGroupD0::active_intervention, "Apply an active sampling intervention")

    .property("vitals", &KoalaGroupD0::get_vitals, "Get vitals (yeay/day/alive)")
    .property("state", &KoalaGroupD0::get_state, &KoalaGroupD0::set_state, "Get and set state")
    .property("parameters", &KoalaGroupD0::get_pars, &KoalaGroupD0::set_pars, "Get and set parameters")
  ;

}
