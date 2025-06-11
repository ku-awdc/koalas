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


template <bool s_t>
class Simple
{
public:
  int i = 0;

  Simple(int ii, bool fake) : i(ii)
  {

  }

  Simple(Simple const& obj)
  {
    Rcpp::Rcout << "Copying...\n";
    i = obj.i+1;
  }

  // TODO: cloning can cause R to abort if the original pointer is invalid
  Simple clone() const
  {
    return Simple(*this);
  }

};
using SimpleTrue = Simple<true>;
RCPP_EXPOSED_AS(SimpleTrue)
RCPP_EXPOSED_WRAP(SimpleTrue)

int testfun()
{
  SimpleTrue tt = SimpleTrue(0, true);
  SimpleTrue copy = SimpleTrue(tt);
  return copy.i;
}


using KoalaGroupD3 = KoalaGroup<debug, 3U, 3U, 3U, 3U, 3U>;
RCPP_EXPOSED_AS(KoalaGroupD3)
RCPP_EXPOSED_WRAP(KoalaGroupD3)
// Note: can't use RCPP_EXPOSED_CLASS with aliases


RCPP_MODULE(koalas){

	using namespace Rcpp;

  function("testfun", &testfun);

  class_<Simple<true>>("Simple")
    .factory(invalidate_default_constructor)
    .constructor<int, bool>("C'tor")
    .constructor<SimpleTrue&>("Copy c'tor")
    .field("i", &Simple<true>::i)
    .method("clone", &Simple<true>::clone)
  ;

  class_<KoalaGroupD3>("KoalaGroupD3")
    .factory(invalidate_default_constructor)
    .constructor<IntegerVector const, NumericVector const, NumericVector const>("Constructor with ncomps vector, parameter vector and state vector")
    //.constructor<KoalaGroupD3&>("Copy constructor")
    .method("clone", &KoalaGroupD3::clone, "Clone method")

    .method("update", &KoalaGroupD3::update, "Update method")
    .method("active_intervention", &KoalaGroupD3::active_intervention, "Apply an active sampling intervention")

    .property("vitals", &KoalaGroupD3::get_vitals, "Get vitals (yeay/day/alive)")
    .property("state", &KoalaGroupD3::get_state, &KoalaGroupD3::set_state, "Get and set state")
    .property("parameters", &KoalaGroupD3::get_pars, &KoalaGroupD3::set_pars, "Get and set parameters")
  ;

  using KoalaGroupD1 = KoalaGroup<debug, 1U, 1U, 1U, 1U, 1U>;
  class_<KoalaGroupD1>("KoalaGroupD1")
    .factory(invalidate_default_constructor)
    .constructor<IntegerVector const, NumericVector const, NumericVector const>("Constructor with ncomps vector, parameter vector and state vector")

    .method("update", &KoalaGroupD1::update, "Update method")
    .method("active_intervention", &KoalaGroupD1::active_intervention, "Apply an active sampling intervention")

    .property("vitals", &KoalaGroupD1::get_vitals, "Get vitals (yeay/day/alive)")
    .property("state", &KoalaGroupD1::get_state, &KoalaGroupD1::set_state, "Get and set state")
    .property("parameters", &KoalaGroupD1::get_pars, &KoalaGroupD1::set_pars, "Get and set parameters")
  ;

  using KoalaGroupD0 = KoalaGroup<debug, 0U, 0U, 0U, 0U, 0U>;
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
