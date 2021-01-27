/*!
  \file cargs.cc
  \brief method implementation for command line argument parser class
  \copyright Released under the MIT License.
  Copyright 2021 Cameron Palmer
*/

#include "annotate_frequency/cargs.h"

void annotate_frequency::cargs::initialize_options() {
  _desc.add_options()("help,h", "emit this help message")(
      "input-filename,i", boost::program_options::value<std::string>(),
      "name of input association file")(
      "supercontinent,s", boost::program_options::value<std::string>(),
      "name of supercontinent (if using combined frequency file)")(
      "frequency-filename,f", boost::program_options::value<std::string>(),
      "name of file with reference frequency data")(
      "frequency-metadata-filename,m",
      boost::program_options::value<std::string>()->default_value(""),
      "name of file with reference frequency SNP annotations (optional)")(
      "output-filename,o", boost::program_options::value<std::string>(),
      "name of output results file");
}
