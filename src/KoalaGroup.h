#ifndef KOALA_GROUP_H
#define KOALA_GROUP_H

#include <array>
#include <string_view>
#include <Rcpp.h>

#include "Compartment.h"

template <auto CTS, size_t s_ncomps>
class KoalaGroup{
private:
  
  Compartment<CTS, 1, 1> m_S {};
  Compartment<CTS, s_ncomps, 1> m_V {};
  Compartment<CTS, s_ncomps, 2> m_I;
  Compartment<CTS, s_ncomps, 2> m_N;
  Compartment<CTS, s_ncomps, 1> m_R;
  
  Compartment<CTS, s_ncomps, 2> m_Af;
  Compartment<CTS, 1, 1> m_Cf;

  Compartment<CTS, 1, 1> m_Sf;
  Compartment<CTS, s_ncomps, 1> m_Vf;
  Compartment<CTS, s_ncomps, 2> m_If;
  Compartment<CTS, s_ncomps, 2> m_Nf;
  Compartment<CTS, s_ncomps, 1> m_Rf;
  
  
  double m_vs_rate = 0.0;       // #1
  double m_ni_rate = 0.0;       // #2
  double m_rs_rate = 0.0;       // #3
  double m_beta = 0.0;          // #4
  double m_ia_rate = 0.0;       // #5
  double m_screc_prop = 0.0;    // #6
  double m_disrec_prop = 0.0;   // #7
  double m_birthrate = 0.0;     // #8
  double m_ac_rate = 0.0;       // #9
  double m_mort_nat = 0.0;      // #10
  double m_mort_dis = 0.0;      // #11
  double m_rel_fecundity = 0.0; // #12
  
  
  KoalaGroup() = delete;
  
  [[nodiscard]] auto to_duration(double const rate) const noexcept(!CTS.debug) -> double
  {
    double const duration = 365.0 / rate;
    return duration;
  }

  [[nodiscard]] auto to_rate(double const duration) const noexcept(!CTS.debug) -> double
  {
    double const rate = 365.0 / duration;
    return rate;    
  }
  

public:
  KoalaGroup(Rcpp::NumericVector const parameters) noexcept(!CTS.debug)
  {
    set_pars_natural(parameters);    
  }
  
  auto set_pars_natural(Rcpp::NumericVector const parameters) noexcept(!CTS.debug) -> void
  {
        
  }
  
  [[nodiscard]] auto get_pars_natural() const noexcept(!CTS.debug) -> Rcpp::NumericVector
  {
    using namespace Rcpp;    
    NumericVector pars = NumericVector::create(
      _["vacc_immune_duration"] = to_duration(m_vs_rate),     // #1
      _["vacc_redshed_duration"] = to_duration(m_ni_rate),    // #2
      _["natural_immune_duration"] = to_duration(m_rs_rate),  // #3
      _["beta"] = m_beta,                                     // #4
      _["subcinical_duration"] = to_duration(m_ia_rate),      // #5
      _["subclinical_recover_proportion"] = m_screc_prop,     // #6
      _["diseased_recover_proportion"] = m_disrec_prop,       // #7
      _["birthrate"] = m_birthrate,                           // #8
      _["acute_duration"] = to_duration(m_ac_rate),           // #9
      _["lifespan_natural"] = to_duration(m_mort_nat),        // #10
      _["lifespan_diseased"] = to_duration(m_mort_dis),       // #11
      _["relative_fecundity"] = m_rel_fecundity               // #12
    );
    
    return pars; 
  }
  
  [[nodiscard]] auto get_sum() const noexcept(!CTS.debug) -> double
  {
    double const rv = 2.0;
    return rv;
  }
  
};

#endif // KOALA_GROUP_H
