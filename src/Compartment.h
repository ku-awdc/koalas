#ifndef COMPARTMENT_H
#define COMPARTMENT_H

#include <array>

template <auto CTS, size_t s_ncomps, size_t s_ndests>
class Compartment{
private:
  std::array<double, s_ncomps> m_values { };

public:
  Compartment() noexcept(!CTS.debug) = default;

  Compartment(double const value) noexcept
  {
    set_sum(value);
  }

  [[nodiscard]] auto get_sum() const noexcept(!CTS.debug) -> double
  {
    double const rv = std::accumulate(m_values.begin(), m_values.end(), 0.0);
    return rv;
  }

  auto set_sum(double const value) noexcept(!CTS.debug) -> void
  {
    for(auto& val : m_values){
      val = value / static_cast<double>(s_ncomps);
    }
  }

  auto set_rate(double const rate, std::array<double, s_ndests> const& proportions,
                std::array<double, s_ndests> const& proportions_final) noexcept(!CTS.debug) -> void
  {

  }

  auto set_rate(double const rate, std::array<double, s_ndests> const& proportions) noexcept(!CTS.debug) -> void
  {
    set_rate(rate, proportions, proportions);
  }

  [[nodiscard]] auto take_out(double const prop)
    -> double
  {
    double rv = 0.0;
    for(auto& val : m_values){
      double const tt = val * prop;
      rv += tt;
      val -= tt;
    }
    return rv;
  }

  auto ptr() noexcept(!CTS.debug) ->
    std::array<double, s_ncomps>&
  {
    return m_values;
  }


};

#endif // COMPARTMENT_H
