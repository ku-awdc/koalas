#ifndef KOALA_GROUP_H
#define KOALA_GROUP_H

#include <array>
#include <string_view>
#include <Rcpp.h>

#include "Compartment.h"

template <auto CTS, unsigned int s_nV, unsigned int s_nI, unsigned int s_nN, unsigned int s_nR, unsigned int s_nA>
class KoalaGroup{
private:

  static constexpr unsigned int s_nS = 1U;
  static constexpr unsigned int s_nC = 1U;

  Compartment<CTS, s_nS> m_S;
  Compartment<CTS, s_nV> m_V;
  Compartment<CTS, s_nI> m_I;
  Compartment<CTS, s_nN> m_N;
  Compartment<CTS, s_nR> m_R;

  Compartment<CTS, s_nA> m_Af;
  Compartment<CTS, s_nC> m_Cf;

  Compartment<CTS, s_nS> m_Sf;
  Compartment<CTS, s_nV> m_Vf;
  Compartment<CTS, s_nI> m_If;
  Compartment<CTS, s_nN> m_Nf;
  Compartment<CTS, s_nR> m_Rf;

  double m_Z = 0.0;
  double m_sumTx = 0.0;
  double m_sumVx = 0.0;
  double m_sumRx = 0.0;
  double m_sumMx = 0.0;

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

  double m_se = 1.0;
  double m_sp = 1.0;
  // Destinations: (1) R, (2) N, (3) A/C/I, (4) Z
  double m_tx_dest1 = 0.33;
  double m_tx_dest2 = 0.33;
  double m_tx_dest3 = 0.01;
  double m_tx_dest4 = 0.33;
  double m_vx_eff = 1.0;
  double m_vx_bst = 1.0;
  double m_passive_rate = 0.0;


  KoalaGroup() = delete;

  [[nodiscard]] auto to_duration(double const rate)
    const noexcept(!CTS.debug)
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
    double const dmort = m_Af.take_rate(m_mort_dis, d_time);
    m_sumMx += dmort;
    m_Z += dmort;

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

  auto update_passive(double const d_time) noexcept(!CTS.debug)
    -> void
  {
    if (m_passive_rate > 0.0 )
    {
      double const prop = 1.0 - std::exp(-m_passive_rate * d_time);
      treat_vacc_all(prop);

      update_apply();
    }
  }

  template <typename SrcT, typename DstT>
  auto treat_vacc_noninf(SrcT& src, DstT& dst, double const prop) noexcept(!CTS.debug)
    -> void
  {
    double const test = src.get_sum() * prop;

    // Test positives are also treated but there is no other difference:
    m_sumTx += (1.0 - m_sp) * test;
    m_sumVx += test;

    // Move to V(f):
    dst.insert_value_start( src.take_prop(m_vx_eff * prop) );
  }

  enum class Shedding { negative, positive };

  template <Shedding s_shed, typename SrcT, typename DstT1, typename DstT2>
  auto treat_vacc_inf(SrcT& src, DstT1& dst1, DstT2& dst2, double const prop) noexcept(!CTS.debug)
    -> void
  {
    double const test = src.get_sum() * prop;

    // Test positives are treated and have a chance to be cured:
    double const testpos = test * ((s_shed==Shedding::negative) ? (1.0-m_sp) : m_se);
    double const testneg = test - testpos;
    m_sumTx += testpos;

    double const cure = m_tx_dest1 * testpos;
    double const nshd = m_tx_dest2 * testpos;
    double const back = m_tx_dest3 * testpos;
    m_sumVx += (cure + nshd + back);
    double const remv = test - (cure + nshd + back);
    m_sumRx += remv;


    // Everything except (back+testneg) * (1-m_vx_bst) is removed from src:
    double const restart = (back+testneg) * m_vx_bst;
    src.remove_number( cure+nshd+remv+restart );

    // Go back to the src category N(f) or R(f) but with fresh vaccination:
    src.insert_value_start(restart);
    // Destinations: (1) R, (2) N, (3) A/C/I, (4) Z
    dst1.insert_value_start(cure);
    dst2.insert_value_start(nshd);
    // Dest 3 is not taken out to start with
    m_Z += remv;
  }

  auto treat_vacc_all(double const prop) noexcept(!CTS.debug)
    -> void
  {
    // Susceptible:
    treat_vacc_noninf(m_S, m_V, prop);
    treat_vacc_noninf(m_Sf, m_Vf, prop);

    // Already vaccinated:
    treat_vacc_noninf(m_V, m_V, prop);
    treat_vacc_noninf(m_Vf, m_Vf, prop);

    // Non-shedding:
    treat_vacc_inf<Shedding::negative>(m_N, m_R, m_N, prop);
    treat_vacc_inf<Shedding::negative>(m_Nf, m_Rf, m_Nf, prop);

    // Infectious:
    treat_vacc_inf<Shedding::positive>(m_I, m_R, m_N, prop);
    treat_vacc_inf<Shedding::positive>(m_If, m_Rf, m_Nf, prop);

    // Diseased:
    treat_vacc_inf<Shedding::positive>(m_Af, m_Rf, m_Nf, prop);
    treat_vacc_inf<Shedding::positive>(m_Cf, m_Rf, m_Nf, prop);

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
  KoalaGroup(Rcpp::IntegerVector ncomps, Rcpp::NumericVector const parameters, Rcpp::NumericVector const state) noexcept(!CTS.debug)
    : m_S(s_nS), m_V(ncomps["V"]), m_I(ncomps["I"]), m_N(ncomps["N"]), m_R(ncomps["R"]),
      m_Af(ncomps["A"]), m_Cf(s_nC),
      m_Sf(s_nS), m_Vf(ncomps["V"]), m_If(ncomps["I"]), m_Nf(ncomps["N"]), m_Rf(ncomps["R"])
  {
    set_pars(parameters);
    set_state(state);
  }

  auto set_pars(Rcpp::NumericVector const parameters) noexcept(!CTS.debug)
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

    m_se = parameters["sensitivity"];
    m_sp = parameters["specificity"];
    m_tx_dest1 = parameters["treatment_dest_R"];
    m_tx_dest2 = parameters["treatment_dest_N"];
    m_tx_dest3 = parameters["treatment_dest_IAC"];
    m_tx_dest4 = parameters["treatment_dest_remove"];
    m_vx_eff = parameters["vaccine_efficacy"];
    m_vx_bst = parameters["vaccine_booster"];
    m_passive_rate = parameters["passive_intervention_rate"] / 365.0;
        // -std::log(1.0 - parameters["passive_proportion"]) / 365.0

  }

  [[nodiscard]] auto get_pars() const noexcept(!CTS.debug)
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
      _["relative_fecundity"] = m_rel_fecundity,              // #12

      _["sensitivity"] = m_se,
      _["specificity"] = m_sp,
      _["treatment_dest_R"] = m_tx_dest1,
      _["treatment_dest_N"] = m_tx_dest2,
      _["treatment_dest_IAC"] = m_tx_dest3,
      _["treatment_dest_remove"] = m_tx_dest4,
      _["vaccine_efficacy"] = m_vx_eff,
      _["vaccine_booster"] = m_vx_bst
    ); // Max number of elements is 20

    pars.push_back(
      m_passive_rate*365.0,  // 1.0 - std::exp(-m_passive_rate * 365.0),
      "passive_intervention_rate"
    );

    return pars;
  }

  auto set_state(Rcpp::NumericVector const state) noexcept(!CTS.debug)
    -> void
  {
    m_S.set_sum(state["S"]);
    m_V.set_sum(state["V"]);
    m_I.set_sum(state["I"]);
    m_N.set_sum(state["N"]);
    m_R.set_sum(state["R"]);
    m_Af.set_sum(state["Af"]);
    m_Cf.set_sum(state["Cf"]);
    m_Sf.set_sum(state["Sf"]);
    m_Vf.set_sum(state["Vf"]);
    m_If.set_sum(state["If"]);
    m_Nf.set_sum(state["Nf"]);
    m_Rf.set_sum(state["Rf"]);

    m_year = static_cast<int>(state["Year"]);
    m_day = static_cast<int>(state["Day"]);
    m_sumTx = state["SumTx"];
    m_sumVx = state["SumVx"];
    m_sumRx = state["SumRx"];
    m_sumMx = state["SumMx"];

    m_Z = -(m_S+m_V+m_I+m_N+m_R+m_Af+m_Cf+m_Vf+m_If+m_Nf+m_Rf);
  }

  [[nodiscard]] auto get_state() const noexcept(!CTS.debug)
    -> Rcpp::NumericVector
  {
    using namespace Rcpp;
    NumericVector rv = NumericVector::create(
      _["Year"] = static_cast<double>(m_year),
      _["Day"] = static_cast<double>(m_day),
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
      _["Z"] = m_Z,
      _["SumTx"] = m_sumTx,
      _["SumVx"] = m_sumVx,
      _["SumRx"] = m_sumRx,
      _["SumMx"] = m_sumMx
    );
    return rv;
  }

  auto active_intervention(double const prop) noexcept(!CTS.debug)
    -> void
  {
    treat_vacc_all(prop);
    update_apply();
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
        // d_time is 1.0 as this happens once per day:
        update_passive(1.0);

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

  auto get_vitals() const noexcept(!CTS.debug)
    -> Rcpp::NumericVector
  {
    using namespace Rcpp;
    NumericVector rv = NumericVector::create(
      _["Year"] = static_cast<double>(m_year),
      _["Day"] = static_cast<double>(m_day),
      _["Alive"] = get_fertile() + get_infertile()
    );

    return rv;
  }

};

#endif // KOALA_GROUP_H
