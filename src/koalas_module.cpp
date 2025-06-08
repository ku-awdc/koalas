#include <Rcpp.h>

#include "KoalaGroup.h"

struct Settings
{
  bool debug = true;
};
constexpr Settings debug { .debug = true };
// constexpr Settings run { .debug = false };

/* TODO: link to blofeld for this! */
template <class RcppModuleClassName>
RcppModuleClassName* invalidate_default_constructor() {
  Rcpp::stop("Default constructor is disabled for this class");
  return 0;
} //blofeld


RCPP_MODULE(koalas){

	using namespace Rcpp;

  using KoalaGroupD3 = KoalaGroup<debug, 3, 3, 3, 3, 3>;
  class_<KoalaGroupD3>("KoalaGroupD3")
    .factory(invalidate_default_constructor)
    .constructor<IntegerVector const, NumericVector const, NumericVector const>("Constructor with ncomps vector, parameter vector and state vector")

    .method("update", &KoalaGroupD3::update, "Update method")
    .method("active_intervention", &KoalaGroupD3::active_intervention, "Apply an active sampling intervention")

    .property("vitals", &KoalaGroupD3::get_vitals, "Get vitals (yeay/day/alive)")
    .property("state", &KoalaGroupD3::get_state, &KoalaGroupD3::set_state, "Get and set state")
    .property("parameters", &KoalaGroupD3::get_pars, &KoalaGroupD3::set_pars, "Get and set parameters")
  ;

  using KoalaGroupD1 = KoalaGroup<debug, 1, 1, 1, 1, 1>;
  class_<KoalaGroupD1>("KoalaGroupD1")
    .factory(invalidate_default_constructor)
    .constructor<IntegerVector const, NumericVector const, NumericVector const>("Constructor with ncomps vector, parameter vector and state vector")

    .method("update", &KoalaGroupD1::update, "Update method")
    .method("active_intervention", &KoalaGroupD1::active_intervention, "Apply an active sampling intervention")

    .property("vitals", &KoalaGroupD1::get_vitals, "Get vitals (yeay/day/alive)")
    .property("state", &KoalaGroupD1::get_state, &KoalaGroupD1::set_state, "Get and set state")
    .property("parameters", &KoalaGroupD1::get_pars, &KoalaGroupD1::set_pars, "Get and set parameters")
  ;

  using KoalaGroupD0 = KoalaGroup<debug, 0, 0, 0, 0, 0>;
  class_<KoalaGroupD0>("KoalaGroupD0")
    .factory(invalidate_default_constructor)
    .constructor<IntegerVector const, NumericVector const, NumericVector const>("Constructor with ncomps vector, parameter vector and state vector")

    .method("update", &KoalaGroupD0::update, "Update method")
    .method("active_intervention", &KoalaGroupD0::active_intervention, "Apply an active sampling intervention")

    .property("vitals", &KoalaGroupD0::get_vitals, "Get vitals (yeay/day/alive)")
    .property("state", &KoalaGroupD0::get_state, &KoalaGroupD0::set_state, "Get and set state")
    .property("parameters", &KoalaGroupD0::get_pars, &KoalaGroupD0::set_pars, "Get and set parameters")
  ;

}
