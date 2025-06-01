#include <Rcpp.h>

#include "KoalaGroup.h"

struct Settings
{
  bool debug = true;
};
constexpr Settings debug { .debug = true };
constexpr Settings run { .debug = false };

/* TODO: link to blofeld for this! */
template <class RcppModuleClassName>
RcppModuleClassName* invalidate_default_constructor() {
  Rcpp::stop("Default constructor is disabled for this class");
  return 0;
} //blofeld


RCPP_MODULE(koalas){

	using namespace Rcpp;

  using KoalaGroupD3 = KoalaGroup<debug, 3>;
  class_<KoalaGroupD3>("KoalaGroupD3")
    .factory(invalidate_default_constructor)
    .constructor<NumericVector const>("Constructor with parameter vector")

    .property("parameters", &KoalaGroupD3::get_pars_natural, &KoalaGroupD3::set_pars_natural, "Get and set parameters N")
  ;

}
