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

  double m_Z = 0.0;

  int m_year = 0;
  int m_day = 0;
  double m_time = 0.0;

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

  [[nodiscard]] auto to_duration(double const rate) const noexcept(!CTS.debug)
    -> double
  {
    double const duration = 1.0 / (rate * 365.0);
    return duration;
  }

  [[nodiscard]] auto to_rate(double const duration) const noexcept(!CTS.debug)
    -> double
  {
    double const rate = 1.0 / (duration * 365.0);
    return rate;
  }

  auto update_death(double const d_time) noexcept(!CTS.debug)
    -> void
  {
    m_Z += m_S.take_rate(m_mort_nat, d_time);
    m_Z += m_V.take_rate(m_mort_nat, d_time);
    m_Z += m_I.take_rate(m_mort_nat, d_time);
    m_Z += m_N.take_rate(m_mort_nat, d_time);
    m_Z += m_R.take_rate(m_mort_nat, d_time);
    m_Z += m_Cf.take_rate(m_mort_nat, d_time);
    m_Z += m_Sf.take_rate(m_mort_nat, d_time);
    m_Z += m_Vf.take_rate(m_mort_nat, d_time);
    m_Z += m_If.take_rate(m_mort_nat, d_time);
    m_Z += m_Nf.take_rate(m_mort_nat, d_time);
    m_Z += m_Rf.take_rate(m_mort_nat, d_time);

    // Mortality for Af is different:
    m_Z += m_Af.take_rate(m_mort_dis, d_time);

    update_apply();
  }

  auto update_birth(double const d_time) noexcept(!CTS.debug)
    -> void
  {
    double const births = m_birthrate * (get_fertile() + ((1.0 - m_rel_fecundity) * get_infertile())) * d_time;
    // Note: this is correct to be value not rate!
    m_S.insert_value_start(births);
    m_Z -= births;

    update_apply();
  }

  auto update_disease(double const d_time) noexcept(!CTS.debug)
    -> void
  {
    // V(f) -> S(f)
    m_S.insert_value_start( m_V.carry_rate(m_vs_rate, d_time) );
    m_Sf.insert_value_start( m_Vf.carry_rate(m_vs_rate, d_time) );

    // N(f) -> I(f)
    m_I.insert_value_start( m_N.carry_rate(m_ni_rate, d_time) );
    m_If.insert_value_start( m_Nf.carry_rate(m_ni_rate, d_time) );

    // R(f) -> S(f)
    m_S.insert_value_start( m_R.carry_rate(m_rs_rate, d_time) );
    m_Sf.insert_value_start( m_Rf.carry_rate(m_rs_rate, d_time) );

    // S(f) -> I(f)
    {
      double const infrate = m_beta * get_infected() / (get_fertile() + get_infertile());
      {
        double const leaveS = m_S.take_rate(infrate, d_time);
        m_I.insert_value_start(leaveS);
      }
      {
        double const leaveSf = m_Sf.take_rate(infrate, d_time);
        m_If.insert_value_start(leaveSf);
      }
    }

    // I(f) -> Af/R(f)
    {
      double const leaveI = m_I.carry_rate(m_ia_rate, d_time);
      double const toR = leaveI * m_screc_prop;
      m_R.insert_value_start(toR);
      m_Af.insert_value_start(leaveI - toR);
    }
    {
      double const leaveIf = m_If.carry_rate(m_ia_rate, d_time);
      double const toRf = leaveIf * m_screc_prop;
      m_Rf.insert_value_start(toRf);
      m_Af.insert_value_start(leaveIf - toRf);
    }

    // Af -> Cf
    // TODO: checl m_disrec_prop is 0!
    m_Cf.insert_value_start( m_Af.carry_rate(m_ac_rate, d_time) );

    update_apply();
  }

  auto update_passive() noexcept(!CTS.debug)
    -> void
  {

    update_apply();
  }

  auto update_apply() noexcept(!CTS.debug)
    -> void
  {
    m_S.apply_changes();
    m_V.apply_changes();
    m_N.apply_changes();
    m_I.apply_changes();
    m_R.apply_changes();

    m_Af.apply_changes();
    m_Cf.apply_changes();

    m_Sf.apply_changes();
    m_Vf.apply_changes();
    m_Nf.apply_changes();
    m_If.apply_changes();
    m_Rf.apply_changes();
  }

  [[nodiscard]] auto get_infected()
    const noexcept(!CTS.debug)
    -> double
  {
    double const infected = m_I + m_If + m_Af + m_Cf;
    return infected;
  }

  [[nodiscard]] auto get_fertile()
    const noexcept(!CTS.debug)
    -> double
  {
    double const fertile = m_S + m_V + m_I + m_N + m_R;
    return fertile;
  }

  [[nodiscard]] auto get_infertile()
    const noexcept(!CTS.debug)
    -> double
  {
    double const infertile = m_Af + m_Cf + m_Sf + m_Vf + m_If + m_Nf + m_Rf;
    return infertile;
  }


public:
  KoalaGroup(Rcpp::NumericVector const parameters) noexcept(!CTS.debug)
  {
    set_pars_natural(parameters);
    m_S.set_sum(255.0);
    m_Cf.set_sum(45.0);
    m_Z = -300.0;
  }

  auto set_pars_natural(Rcpp::NumericVector const parameters) noexcept(!CTS.debug)
    -> void
  {
    m_vs_rate = to_rate(parameters["vacc_immune_duration"]);     // #1
    m_ni_rate = to_rate(parameters["vacc_redshed_duration"]);    // #2
    m_rs_rate = to_rate(parameters["natural_immune_duration"]);  // #3
    m_beta = parameters["beta"]/365.0;                           // #4
    m_ia_rate = to_rate(parameters["subcinical_duration"]);      // #5
    m_screc_prop = parameters["subclinical_recover_proportion"]; // #6
    m_disrec_prop = parameters["diseased_recover_proportion"];   // #7
    m_birthrate = parameters["birthrate"]/365.0;                 // #8
    m_ac_rate = to_rate(parameters["acute_duration"]);           // #9
    m_mort_nat = to_rate(parameters["lifespan_natural"]);        // #10
    m_mort_dis = to_rate(parameters["lifespan_diseased"]);       // #11
    m_rel_fecundity = parameters["relative_fecundity"];          // #12
  }

  [[nodiscard]] auto get_pars_natural() const noexcept(!CTS.debug)
    -> Rcpp::NumericVector
  {
    using namespace Rcpp;
    NumericVector pars = NumericVector::create(
      _["vacc_immune_duration"] = to_duration(m_vs_rate),     // #1
      _["vacc_redshed_duration"] = to_duration(m_ni_rate),    // #2
      _["natural_immune_duration"] = to_duration(m_rs_rate),  // #3
      _["beta"] = m_beta*365.0,                               // #4
      _["subcinical_duration"] = to_duration(m_ia_rate),      // #5
      _["subclinical_recover_proportion"] = m_screc_prop,     // #6
      _["diseased_recover_proportion"] = m_disrec_prop,       // #7
      _["birthrate"] = m_birthrate*365.0,                     // #8
      _["acute_duration"] = to_duration(m_ac_rate),           // #9
      _["lifespan_natural"] = to_duration(m_mort_nat),        // #10
      _["lifespan_diseased"] = to_duration(m_mort_dis),       // #11
      _["relative_fecundity"] = m_rel_fecundity               // #12
    );

    return pars;
  }

  auto update(int const days, double const d_time) noexcept(!CTS.debug)
    -> void
  {
    int newdays = 0;
    while (newdays < days)
    {
      update_death(d_time);
      update_birth(d_time);
      update_disease(d_time);

      m_time += d_time;
      if(m_time >= 1.0)
      {
        update_passive();

        m_time = 0.0;
        m_day++;
        newdays++;
        if(m_day >= 366)
        {
          m_day = 1;
          m_year++;
        }
      }

      Rcpp::checkUserInterrupt();
    };

  }

  [[nodiscard]] auto get_state() const noexcept(!CTS.debug)
    -> Rcpp::DataFrame
  {
    using namespace Rcpp;
    DataFrame rv = DataFrame::create(
      _["Year"] = m_year,
      _["Day"] = m_day,
      _["S"] = m_S.get_sum(),
      _["V"] = m_V.get_sum(),
      _["I"] = m_I.get_sum(),
      _["N"] = m_N.get_sum(),
      _["R"] = m_R.get_sum(),
      _["Af"] = m_Af.get_sum(),
      _["Cf"] = m_Cf.get_sum(),
      _["Sf"] = m_Sf.get_sum(),
      _["Vf"] = m_Vf.get_sum(),
      _["If"] = m_If.get_sum(),
      _["Nf"] = m_Nf.get_sum(),
      _["Rf"] = m_Rf.get_sum(),
      _["Z"] = m_Z
    );
    return rv;
  }

  auto set_state(Rcpp::DataFrame const state) noexcept(!CTS.debug)
    -> void
  {
    // TODO
  }

};

#endif // KOALA_GROUP_H
