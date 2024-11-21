#pragma once

#include "t8_element.hxx"

namespace t8_mra {

template <>
typename t8_multiscale<T8_ECLASS_TRIANGLE>::lmi_t
t8_multiscale<T8_ECLASS_TRIANGLE>::t8_element_id_to_lmi(
    const t8_element_t* element, t8_eclass_scheme_c* eclass_scheme,
    t8_locidx_t element_id, size_t offset) {
  levelmultiindex<DIM> lmi;
  lmi.level() = eclass_scheme->t8_element_level(element);
  lmi.mi(0) = element_id / std::pow(lmi[0], 4);

  auto tmp = lmi.mi(0);  /// Path starts with basecell

  for (auto l = 0; l < static_cast<int>(lmi.level()); ++l)
    tmp +=
        std::pow(l + 1 + offset, 10) *
        (static_cast<size_t>(
             element_id / std::pow(static_cast<int>(lmi.level()) - l - 1, 4)) %
         4);

  lmi.mi(1) = tmp;

  return lmi;
}

}  // namespace t8_mra
