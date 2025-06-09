#ifndef COMPARTMENT_H
#define COMPARTMENT_H

#include <array>

// TODO: should take/carry be immediate?  So only insert is done on delay?
// TODO: allow s_ncomps = 0 in which case it is set at run time

template <auto CTS, unsigned int s_ncomps>//, unsigned int s_ndests>
class Compartment{
private:
  // TODO: these need to be array or vector, depending on s_ncomps==0
  std::array<double, s_ncomps> m_values { };
  std::array<double, s_ncomps> m_changes { };
  
  const int m_ncomps;

  auto check_compartments() const noexcept(!CTS.debug)
    -> void
  {
    if constexpr (!CTS.debug && s_ncomps!=0U)
    {
      if (m_ncomps != s_ncomps) Rcpp::stop("Non-matching number of sub-compartments");
    }
    if constexpr (s_ncomps <= 0U)
    {
      Rcpp::stop("Negative s_ncomps are not allowed, and values of 0 are not yet supported");
    }    
  }

  Compartment() = delete;
  
public:
  explicit Compartment(int const ncomps) noexcept(!CTS.debug)
    : m_ncomps(ncomps)
  {
    check_compartments();
  }

  Compartment(int const ncomps, double const value) noexcept(!CTS.debug)
    : m_ncomps(ncomps)
  {
    set_sum(value);
    check_compartments();
  }

  auto apply_changes()
    noexcept(!CTS.debug)
    -> void
  {
    for(int i=0; i<s_ncomps; ++i)
    {
      m_values[i] += m_changes[i];
      m_changes[i] = 0.0;
    }
  }

  [[nodiscard]] auto get_sum()
    const noexcept(!CTS.debug)
    -> double
  {
    double const rv = std::accumulate(m_values.begin(), m_values.end(), 0.0);
    return rv;
  }

  auto set_sum(double const value) noexcept(!CTS.debug)
    -> void
  {
    // TOODO: switch to distribute balanced or all in first box
    for(auto& val : m_values){
      // val = value / static_cast<double>(s_ncomps);
      val = 0.0;
    }
    m_values[0] = value;
  }

  /*
  auto set_rate(double const rate, std::array<double, s_ndests> const& proportions,
                std::array<double, s_ndests> const& proportions_final)
    noexcept(!CTS.debug)
    -> void
  {

  }

  auto set_rate(double const rate, std::array<double, s_ndests> const& proportions) noexcept(!CTS.debug)
    -> void
  {
    set_rate(rate, proportions, proportions);
  }
  */

  [[nodiscard]] auto take_rate(double const rate, double const d_time) noexcept(!CTS.debug)
    -> double
  {
    double const prop = 1.0 - std::exp(-rate * s_ncomps * d_time);
    double const rv = take_prop(prop);
    return rv;
  }

  auto remove_number(double const number) noexcept(!CTS.debug)
    -> void
  {
    if (number > 0.0)
    {
      double const propr = 1.0 - (number / get_sum());
      if constexpr (CTS.debug)
      {
        if (propr > 1.0) Rcpp::stop("Negative number supplied to remove from compartment");
        if (propr < 0.0) Rcpp::stop("Attempt to remove too high a number from compartment");
      }
      for(auto& val : m_values){
        val *= propr;
      }
    }
  }

  [[nodiscard]] auto take_prop(double const prop) noexcept(!CTS.debug)
    -> double
  {
    double rv = 0.0;
    for(int i=0; i<s_ncomps; ++i){
      double const tt = m_values[i] * prop;
      rv += tt;
      m_changes[i] -= tt;
    }
    return rv;
  }

  [[nodiscard]] auto carry_rate(double const rate, double const d_time)
    noexcept(!CTS.debug)
    -> double
  {
    double const prop = 1.0 - std::exp(-rate * s_ncomps * d_time);
    double carry = 0.0;
    for(int i=0; i<s_ncomps; ++i){
      m_changes[i] += carry;
      carry = m_values[i] * prop;
      m_changes[i] -= carry;
    }
    return carry;
  }

  auto insert_value_start(double const value)
    noexcept(!CTS.debug)
    -> void
  {
    m_changes[0] += value;
  }

  auto ptr()
    noexcept(!CTS.debug)
    -> std::array<double, s_ncomps>&
  {
    return m_values;
  }

  auto ptr()
    const noexcept(!CTS.debug)
    -> std::array<double, s_ncomps> const&
  {
    return m_values;
  }

  // Note: unusual + overloads return double
  [[nodiscard]] auto operator+(double const sum)
    const noexcept(!CTS.debug)
      -> double
  {
    return sum + get_sum();
  }
  template <auto sc, unsigned int ss>//, size_t sn>
  [[nodiscard]] auto operator+(Compartment<sc, ss> const& obj)
    const noexcept(!CTS.debug)
      -> double
  {
    return obj.get_sum() + get_sum();
  }

};

template <auto sc, unsigned int ss>//, size_t sn>
[[nodiscard]] auto operator+(double const sum, Compartment<sc, ss> const& obj)
  -> double
  {
    return sum + obj.get_sum();
  }

template <auto sc, unsigned int ss>//, size_t sn>
[[nodiscard]] auto operator+(Compartment<sc, ss> const& obj, double const sum)
  -> double
  {
    return sum + obj.get_sum();
  }

#endif // COMPARTMENT_H
