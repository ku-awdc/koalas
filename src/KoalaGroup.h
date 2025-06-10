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
  double m_mort_acute = 0.0;    // #11
  double m_mort_chronic = 0.0;  // #11
  double m_rel_fecundity = 0.0; // #12

  double m_se = 1.0;
  double m_sp = 1.0;
  // Destinations: (1) R, (2) N, (3) A/C/I, (4) Z
  double m_cure_N = 0.9;
  double m_cure_I = 0.9;
  double m_cure_A = 0.75;
  double m_cure_C = 0.6;
  double m_vx_eff = 1.0;
  double m_vx_bst = 1.0;
  double m_passive_rate = 0.0;

  // Assume that passive cull rates are fixed and hard-coded for ease:
  static constexpr double const s_passive_cull_positive = 0.0;
  static constexpr double const s_passive_cull_acute = 0.2;
  static constexpr double const s_passive_cull_chronic = 0.3;

  bool m_recording = false;

  KoalaGroup() = delete;

  [[nodiscard]] auto to_duration(double const rate)
    const noexcept(!CTS.debug)
    -> double
  {
    if (rate == 0.0) return R_PosInf;

    double const duration = 1.0 / (rate * 365.0);
    return duration;
  }

  [[nodiscard]] auto to_rate(double const duration) const noexcept(!CTS.debug)
    -> double
  {
    if (Rcpp::traits::is_infinite<REALSXP>(duration)) return 0.0;

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
    m_Z += m_Sf.take_rate(m_mort_nat, d_time);
    m_Z += m_Vf.take_rate(m_mort_nat, d_time);
    m_Z += m_If.take_rate(m_mort_nat, d_time);
    m_Z += m_Nf.take_rate(m_mort_nat, d_time);
    m_Z += m_Rf.take_rate(m_mort_nat, d_time);

    // Mortality for Af and Cf is different:
    double const dmort = m_Af.take_rate(m_mort_acute, d_time) + m_Cf.take_rate(m_mort_chronic, d_time);
    m_sumMx += dmort;
    m_Z += dmort;

    update_apply();
  }

  auto update_birth(double const d_time) noexcept(!CTS.debug)
    -> void
  {
    double const births = m_birthrate * (get_fertile() + (m_rel_fecundity * get_infertile())) * d_time;
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

  auto update_passive(double const d_time)
    noexcept(!CTS.debug)
    -> void
  {
    if (m_passive_rate > 0.0 )
    {
      double const prop = 1.0 - std::exp(-m_passive_rate * d_time);
      treat_vacc_all(prop, s_passive_cull_positive, s_passive_cull_acute, s_passive_cull_chronic);

      update_apply();
    }
  }

  template <typename SrcT, typename DstT>
  auto treat_vacc_noninf(SrcT& src, DstT& dst, double const prop, double const efficacy) noexcept(!CTS.debug)
    -> void
  {
    double const test = src.get_sum() * prop;

    // Test positives are treated+vaccinated rather than just vaccinated, but there is no other difference:
    m_sumTx += (1.0 - m_sp) * test;
    m_sumVx += m_sp * test;

    // Move to or restart in V(f) or R(f):
    dst.insert_value_start( src.take_prop(efficacy * prop) );
  }

  enum class Status { nonshedding, subclinical, diseased };

  template <Status s_status, typename Src, typename DstF, typename DstN, typename DstR>
  auto treat_vacc_inf(Src& src, DstF& dstF, DstN& dstN, DstR& dstR,
                      double const prop, double const cull_prob,
                      double const cure_prob) noexcept(!CTS.debug)
    -> void
  {
    double const n_sampled = src.get_sum() * prop;

    double progressed = prop;

    // First round of testing - perfect for Af/Cf, se for I/If, 1-sp for N/Nf:
    double const testpos = (s_status==Status::diseased ? 1.0 : (s_status==Status::subclinical ? m_se : (1.0-m_sp)));
    // All animals that are test negative are vaccinated, but this only affects non-shedding:
    double const to_vacc = progressed * (1.0-testpos);
    m_sumVx += (n_sampled * to_vacc);
    if constexpr (s_status==Status::nonshedding)
    {
      dstN.insert_value_start(src.take_prop(m_vx_bst * to_vacc));
    }

    progressed *= testpos;

    // A proportion of test-positive animals are culled:
    double const to_cull = progressed * cull_prob;
    {
      double const n_culled = src.take_prop(to_cull);
      m_Z += n_culled;
      m_sumRx += n_culled;
    }

    progressed *= (1.0 - cull_prob);

    // All animals remaining end up being treated:
    m_sumTx += (n_sampled * progressed);

    // Most have a complete treatment course but some are released early due to a false negative:
    double const treat_complete = (s_status==Status::nonshedding ? 1.0 : m_se); // Note: deliberately not 1-sp here for N!!!
    double const to_frel = progressed * (1.0-treat_complete);
    if constexpr (s_status==Status::diseased)
    {
      // Move diseased to start of I (failed treatment, but no longer clinically diseased)
      dstF.insert_value_start(src.take_prop(to_frel));
    }
    else if constexpr (s_status==Status::subclinical)
    {
      // Do nothing - don't restart animals in I
    }
    else if constexpr (s_status==Status::nonshedding)
    {
      // Restart non-shedding due to vaccine (this will have no effect as test sensitivity above is 1)
      dstN.insert_value_start(src.take_prop(m_vx_bst * to_frel));
    }

    progressed *= treat_complete;

    // Of the complete treatment number, a proportion recover and the rest are non-sheding:
    double const to_nshed = progressed * (1.0 - cure_prob);
    if constexpr (s_status==Status::nonshedding) {
      // For N, restart due to vaccine:
      dstN.insert_value_start(src.take_prop(m_vx_bst * to_nshed));
    }else{
      // Otherwise all move to start of N:
      dstN.insert_value_start(src.take_prop(to_nshed));
    }

    // All cured animals go to start of R:
    double const to_cure = progressed * cure_prob;
    dstR.insert_value_start(src.take_prop(to_cure));

    // Final checks:
    if constexpr (CTS.debug)
    {
      if(std::abs(prop - (to_vacc+to_cull+to_frel+to_nshed+to_cure)) > CTS.tol){
        Rcpp::Rcout << prop << " != " << to_vacc << " + " << to_cull << " + " << to_frel << " + " << to_nshed << " + " << to_cure << "\n";
        Rcpp::stop("Logic error in treat_vacc_inf");
      }
    }

  }

  auto treat_vacc_all(double const prop, double const cull_positive,
                      double const cull_acute, double const cull_chronic)
    noexcept(!CTS.debug)
    -> void
  {
    // Susceptible:
    treat_vacc_noninf(m_S, m_V, prop, m_vx_eff);
    treat_vacc_noninf(m_Sf, m_Vf, prop, m_vx_eff);

    // Already vaccinated:
    treat_vacc_noninf(m_V, m_V, prop, m_vx_bst);
    treat_vacc_noninf(m_Vf, m_Vf, prop, m_vx_bst);

    // Already recovered:
    treat_vacc_noninf(m_R, m_R, prop, m_vx_bst);
    treat_vacc_noninf(m_Rf, m_Rf, prop, m_vx_bst);

    // Non-shedding:
    treat_vacc_inf<Status::nonshedding>(m_N, m_N, m_N, m_R, prop, cull_positive, m_cure_N);
    treat_vacc_inf<Status::nonshedding>(m_Nf, m_Nf, m_Nf, m_Rf, prop, cull_positive, m_cure_N);

    // Infectious:
    treat_vacc_inf<Status::subclinical>(m_I, m_I, m_N, m_R, prop, cull_positive, m_cure_I);
    treat_vacc_inf<Status::subclinical>(m_If, m_If, m_Nf, m_Rf, prop, cull_positive, m_cure_I);

    // Diseased:
    treat_vacc_inf<Status::diseased>(m_Af, m_If, m_Nf, m_Rf, prop, cull_acute, m_cure_A);
    treat_vacc_inf<Status::diseased>(m_Cf, m_If, m_Nf, m_Rf, prop, cull_chronic, m_cure_C);

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

    check_state();
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

  auto check_state()
    const noexcept(!CTS.debug)
    -> void
  {
    if constexpr (CTS.debug)
    {
      if(m_S.get_sum() < 0.0) Rcpp::stop("Negative value in S");
      if(m_V.get_sum() < 0.0) Rcpp::stop("Negative value in V");
      if(m_I.get_sum() < 0.0) Rcpp::stop("Negative value in I");
      if(m_N.get_sum() < 0.0) Rcpp::stop("Negative value in N");
      if(m_R.get_sum() < 0.0) Rcpp::stop("Negative value in R");

      if(m_Af.get_sum() < 0.0) Rcpp::stop("Negative value in Af");
      if(m_Cf.get_sum() < 0.0) Rcpp::stop("Negative value in Cf");

      if(m_Sf.get_sum() < 0.0) Rcpp::stop("Negative value in Sf");
      if(m_Vf.get_sum() < 0.0) Rcpp::stop("Negative value in Vf");
      if(m_If.get_sum() < 0.0) Rcpp::stop("Negative value in If");
      if(m_Nf.get_sum() < 0.0) Rcpp::stop("Negative value in Nf");
      if(m_Rf.get_sum() < 0.0) Rcpp::stop("Negative value in Rf");

      double const N = get_fertile() + get_infertile();
      if(std::abs(N + m_Z) > CTS.tol){
        Rcpp::Rcout << -m_Z << " != " << get_fertile() << " + " << get_infertile() << "\n";
        Rcpp::stop("-Z != Fertile + Infertile");
      }
    }
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
    m_mort_acute = to_rate(parameters["lifespan_acute"]);        // #11
    m_mort_chronic = to_rate(parameters["lifespan_chronic"]);    // #11
    m_rel_fecundity = parameters["relative_fecundity"];          // #12

    m_se = parameters["sensitivity"];
    m_sp = parameters["specificity"];
    m_cure_N = parameters["cure_prob_N"];
    m_cure_I = parameters["cure_prob_I"];
    m_cure_A = parameters["cure_prob_A"];
    m_cure_C = parameters["cure_prob_C"];
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
      _["lifespan_acute"] = to_duration(m_mort_acute),        // #11
      _["lifespan_chronic"] = to_duration(m_mort_chronic),    // #11
      _["relative_fecundity"] = m_rel_fecundity,              // #12

      _["sensitivity"] = m_se,
      _["specificity"] = m_sp,
      _["cure_prob_N"] = m_cure_N,
      _["cure_prob_I"] = m_cure_I,
      _["cure_prob_A"] = m_cure_A,
      _["cure_prob_C"] = m_cure_C
    ); // Max number of elements is 20

    pars.push_back(
      m_vx_eff,
      "vaccine_efficacy"
    );

    pars.push_back(
      m_vx_bst,
      "vaccine_booster"
    );

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

    m_day = static_cast<int>(state["Day"]);
    m_sumTx = state["SumTx"];
    m_sumVx = state["SumVx"];
    m_sumRx = state["SumRx"];
    m_sumMx = state["SumMx"];

    m_Z = -(m_S+m_V+m_I+m_N+m_R+m_Af+m_Cf+m_Vf+m_If+m_Nf+m_Rf);
    check_state();

  }

  [[nodiscard]] auto get_state() const noexcept(!CTS.debug)
    -> Rcpp::NumericVector
  {
    using namespace Rcpp;
    NumericVector rv = NumericVector::create(
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

  auto active_intervention(double const prop, double const cull_positive,
                            double const cull_acute, double const cull_chronic)
    noexcept(!CTS.debug)
    -> void
  {
    treat_vacc_all(prop, cull_positive, cull_acute, cull_chronic);
    update_apply();
  }

  auto update(int const days, double const d_time, bool const record) noexcept(!CTS.debug)
    -> Rcpp::List
  {
    if constexpr (CTS.debug)
    {
      if(days <= 0) Rcpp::stop("Invalid days argument");
      if(m_recording && !record) Rcpp::stop("Can't stop recording when already started!");
    }

    int const len = (m_recording ? days : days+1);
    // Rcpp::Rcout << "Updating for " << days << " days (" << len << " rows)\n";
    Rcpp::List rv(len);

    int ii = 0;
    if(!m_recording){
      Rcpp::NumericVector state = get_state();
      rv[ii] = state;
      ii++;
      
      // Reset cumulative counters:
      m_sumTx = 0.0;
      m_sumVx = 0.0;
      m_sumRx = 0.0;
      m_sumMx = 0.0;
    }
    m_recording = record;

    check_state();

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

        m_time -= 1.0;
        m_day++;
        newdays++;

        Rcpp::NumericVector state = get_state();

        if constexpr (CTS.debug){
          if(ii >= rv.size()) Rcpp::stop("Logic error in update");
        }

        rv[ii] = state;
        ii++;
      }

      check_state();
      Rcpp::checkUserInterrupt();
    };

    check_state();

    return rv;
  }

  auto get_vitals() const noexcept(!CTS.debug)
    -> Rcpp::NumericVector
  {
    check_state();

    using namespace Rcpp;
    NumericVector rv = NumericVector::create(
      _["Day"] = static_cast<double>(m_day),
      _["Alive"] = get_fertile() + get_infertile()
    );

    return rv;
  }

};

#endif // KOALA_GROUP_H
