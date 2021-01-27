/*!
  \file finter.h
  \brief main entry point for finter classes. Pulls in other headers
  and provides primary IO utility functions.
  \copyright Released under the MIT License.
  Copyright 2020 Cameron Palmer
 */

#ifndef FINTER_FINTER_FINTER_H_
#define FINTER_FINTER_FINTER_H_

#include <iostream>
#include <string>

#include "finter/finter_reader.h"
#include "finter/finter_writer.h"

namespace annotate_frequency {
/*!
  \brief utility function to facilitate file input while respecting compression
  suffixes
  @param filename name of file to be opened
  \return pointer to allocated object
  \warning compression type is determined from suffix of provided file
  \warning memory must be managed by caller
 */
inline finter_reader *reconcile_reader(const std::string &filename) {
  finter_reader *ptr = 0;
  if (filename.find(".gz") == filename.size() - 3) {
    ptr = new finter_reader_gzip;
    ptr->open(filename);
  } else if (filename.find(".bz2") == filename.size() - 4) {
    ptr = new finter_reader_bzip2;
    ptr->open(filename);
  } else {
    ptr = new finter_reader_flat;
    ptr->open(filename);
  }
  return ptr;
}
/*!
  \brief utility function to facilitate file output while respecting compression
  suffixes
  @param filename name of file to be opened
  \return pointer to allocated object
  \warning compression type is determined from suffix of provided file
  \warning memory must be managed by caller
 */
inline finter_writer *reconcile_writer(const std::string &filename) {
  finter_writer *ptr = 0;
  if (filename.find(".gz") == filename.size() - 3) {
    ptr = new finter_writer_gzip;
    ptr->open(filename);
  } else if (filename.find(".bz2") == filename.size() - 4) {
    ptr = new finter_writer_bzip2;
    ptr->open(filename);
  } else {
    ptr = new finter_writer_flat;
    ptr->open(filename);
  }
  return ptr;
}
}  // namespace annotate_frequency

#endif  // FINTER_FINTER_FINTER_H_
