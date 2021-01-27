/*!
  \file finter_reader_flat.cc
  \brief method implementations for finter flat file parsing

  This is really just an unnecessarily complicated wrapper for ifstream.
  But I wanted to have finter be capable of handling plaintext files
  without external library support, just for coding convenience down the line.

  \copyright Released under the MIT License.
  Copyright 2020 Cameron Palmer
 */

#include "finter/finter_reader_flat.h"

void annotate_frequency::finter_reader_flat::open(const char *filename) {
  _input.open(filename);
  if (!_input.is_open())
    throw std::domain_error(
        "annotate_frequency::finter_reader_flat::open: cannot open file "
        "\"" +
        std::string(filename) + "\"");
}

bool annotate_frequency::finter_reader_flat::getline(std::string *line) {
  if (!line) throw std::runtime_error("flat_reader::getline: null pointer");
  *line = "";
  if (_input.peek() == EOF) return false;
  std::getline(_input, *line);
  return true;
}

void annotate_frequency::finter_reader_flat::read(char *buf,
                                                  std::streamsize n) {
  _input.read(buf, n);
}
