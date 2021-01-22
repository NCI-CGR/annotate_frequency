/*!
  \file helper.h
  \brief finter utility functions
  \copyright Released under the MIT License.
  Copyright 2020 Cameron Palmer
 */

#ifndef FINTER_FINTER_HELPER_H_
#define FINTER_FINTER_HELPER_H_

#include <sstream>
#include <stdexcept>
#include <string>

namespace annotate_frequency {
/*!
  \brief get platform-specific newline character
  \return platform-specific newline character
 */
inline std::string get_newline() {
#ifdef _WIN64
  return "\r\n";
#elif _WIN32
  return "\r\n";
#else
  return "\n";
#endif
}
}  // namespace annotate_frequency
#endif  // FINTER_FINTER_HELPER_H_
