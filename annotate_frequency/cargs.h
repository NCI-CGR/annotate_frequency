/*!
 \file cargs.h
 \brief command line argument handling
 \author Cameron Palmer
 \note requires boost::program_options library + headers
 \copyright Released under the MIT License.
 Copyright 2021 Cameron Palmer

 Thanks to
 https://www.boost.org/doc/libs/1_70_0/doc/html/program_options/tutorial.html#id-1.3.32.4.3
 */

#ifndef ANNOTATE_FREQUENCY_CARGS_H_
#define ANNOTATE_FREQUENCY_CARGS_H_

#include <stdexcept>
#include <string>
#include <vector>
#include "boost/program_options.hpp"

namespace annotate_frequency {
/*!
  \class cargs
  \brief command line argument parser using boost::program_options
 */
class cargs {
 public:
  /*!
    \brief constructor with program arguments
    @param argc number of arguments including program name
    @param argv string array containing actual arguments
   */
  cargs(int argc, char **argv) : _desc("Recognized options") {
    initialize_options();
    boost::program_options::store(
        boost::program_options::parse_command_line(argc, argv, _desc), _vm);
    boost::program_options::notify(_vm);
  }
  /*!
    \brief copy constructor
    @param obj existing cargs object
   */
  cargs(const cargs &obj) : _desc(obj._desc), _vm(obj._vm) {}
  /*!
    \brief destructor
   */
  ~cargs() throw() {}

  /*!
    \brief set user help documentation and default values for parameters as
    needed

    Note the weird syntax with overloaded () operators. The lists are not
    separated by commas.
   */
  void initialize_options();

  /*!
    \brief access first imputed data info file
    \return name of first imputed data info file, if specified

    Note that this option automatically supports .gz and .bz2 file extensions
   */
  std::string get_set1_imputed_info() const {
    return compute_parameter<std::string>("set1-imputed-info");
  }

  std::string get_input_filename() const {
    return compute_parameter<std::string>("input-filename");
  }

  std::string get_output_filename() const {
    return compute_parameter<std::string>("output-filename");
  }

  std::string get_supercontinent() const {
    return compute_parameter<std::string>("supercontinent");
  }

  std::string get_frequency_filename() const {
    return compute_parameter<std::string>("frequency-filename");
  }

  std::string get_frequency_metadata_filename() const {
    return compute_parameter<std::string>("frequency-metadata-filename");
  }
  /*!
    \brief determine whether the user has requested help documentation
    \return whether the user has requested help documentation

    This test is separate from whether the user has run the software with no
    flags
   */
  bool help() const { return compute_flag("help"); }

  /*!
    \brief find status of arbitrary flag
    @param tag name of flag
    \return whether the flag is set

    This is the underlying accessor function used by the custom flag accessors,
    and can be used for arbitrary flag additions if you don't want to type out
    the custom access functions or want to allow dynamic specification from a
    config file.
   */
  bool compute_flag(const std::string &tag) const { return _vm.count(tag); }
  /*!
    \brief find value of arbitrary parameter
    @tparam value_type class to which the value should be cast
    @param tag name of parameter
    \return value of parameter if specified

    \warning throws exception if parameter was not set and had no default
   */
  template <class value_type>
  value_type compute_parameter(const std::string &tag) const {
    if (_vm.count(tag)) {
      return _vm[tag].as<value_type>();
    }
    throw std::domain_error("cargs: requested parameter \"" + tag + "\" unset");
  }

  /*!
    \brief report help documentation to arbitrary output stream
    @param out stream to which to write help documentation

    Parameter should probably be std::cout or std::cerr at your preference.
   */
  void print_help(std::ostream &out) {
    if (!(out << _desc))
      throw std::domain_error("cargs::print_help: unable to write to stream");
  }

 private:
  /*!
    \brief default constructor
    \warning disabled
   */
  cargs() { throw std::domain_error("cargs: do not use default constructor"); }
  boost::program_options::options_description
      _desc;  //!< help documentation string
  boost::program_options::variables_map
      _vm;  //!< storage of parsed command line settings
};
}  // namespace annotate_frequency

#endif  // ANNOTATE_FREQUENCY_CARGS_H_
