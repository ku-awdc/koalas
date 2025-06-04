#ifndef COMPARTMENT_H
#define COMPARTMENT_H

#include <array>

template <auto CTS, size_t s_ncomps, size_t s_ndests>
class Compartment{
private:
  std::array<double, s_ncomps> m_values { };
  std::array<double, s_ncomps> m_changes { };
  
public:
  Compartment() noexcept(!CTS.debug) = default;

  Compartment(double const value) noexcept(!CTS.debug)
  {
    set_sum(value);
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
    for(auto& val : m_values){
      val = value / static_cast<double>(s_ncomps);
    }
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
      if(i>0){
        m_changes[i] += carry;
      }
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
  template <auto sc, size_t ss, size_t sn>
  [[nodiscard]] auto operator+(Compartment<sc, ss, sn> const& obj)
    const noexcept(!CTS.debug)
      -> double
  {
    return obj.get_sum() + get_sum();
  }

};

template <auto sc, size_t ss, size_t sn>
[[nodiscard]] auto operator+(double const sum, Compartment<sc, ss, sn> const& obj)
  -> double
  {
    return sum + obj.get_sum();
  }

template <auto sc, size_t ss, size_t sn>
[[nodiscard]] auto operator+(Compartment<sc, ss, sn> const& obj, double const sum)
  -> double
  {
    return sum + obj.get_sum();
  }

#endif // COMPARTMENT_H
