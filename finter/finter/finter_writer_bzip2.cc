/*!
  \file finter_writer_bzip2.cc
  \brief method implementations for finter bzip2 writer
  \copyright Released under the MIT License.
  Copyright 2020 Cameron Palmer
 */

#include "finter/finter_writer_bzip2.h"

void annotate_frequency::finter_writer_bzip2::open(const char *filename) {
  if (_raw_output)
    throw std::domain_error(
        "annotate_frequency::finter_writer_bzip2: attempted to reopen "
        "in-use handle");
  _raw_output = fopen(filename, "w");
  if (!_raw_output)
    throw std::domain_error(
        "annotate_frequency::finter_writer_bzip2: cannot open file \"" +
        std::string(filename) + "\"");
  int error = 0;
  _bz_output = BZ2_bzWriteOpen(&error, _raw_output, 9, 0, 0);
  if (error == BZ_CONFIG_ERROR) {
    throw std::domain_error(
        "annotate_frequency::finter_writer_bzip2::open: bzip2 writing "
        "library reports it was compiled improperly");
  } else if (error == BZ_PARAM_ERROR) {
    _fail = true;
  } else if (error != BZ_OK) {
    _bad = true;
  }
}

void annotate_frequency::finter_writer_bzip2::close() {
  int error = 0;
  if (_bz_output) {
    BZ2_bzWriteClose(&error, _bz_output, 0, 0, 0);
    if (error == BZ_SEQUENCE_ERROR) {
      throw std::domain_error(
          "annotate_frequency::finter_writer_bzip2::close: bzip reports "
          "write/close operation called on read handle");
    }
    _bz_output = 0;
  }
  if (_raw_output) {
    fclose(_raw_output);
    _raw_output = 0;
  }
}

void annotate_frequency::finter_writer_bzip2::clear() {
  _good = true;
  _bad = _fail = false;
}

bool annotate_frequency::finter_writer_bzip2::is_open() const {
  return (_raw_output && _bz_output);
}

void annotate_frequency::finter_writer_bzip2::put(char c) {
  int error = 0;
  BZ2_bzWrite(&error, _bz_output, reinterpret_cast<void *>(&c), 1);
  if (error == BZ_PARAM_ERROR || error == BZ_IO_ERROR) {
    _fail = true;
  } else if (error == BZ_SEQUENCE_ERROR) {
    throw std::domain_error(
        "annotate_frequency::finter_writer_bzip2::put: bzip reports "
        "write operation called on read handle");
  }
}

void annotate_frequency::finter_writer_bzip2::writeline(
    const std::string &orig_line) {
  int error = 0;
  std::string line = orig_line + get_newline();
  BZ2_bzWrite(&error, _bz_output,
              const_cast<void *>(reinterpret_cast<const void *>(line.c_str())),
              static_cast<int>(line.size()));
  if (error == BZ_PARAM_ERROR || error == BZ_IO_ERROR) {
    _fail = true;
  } else if (error == BZ_SEQUENCE_ERROR) {
    throw std::domain_error(
        "annotate_frequency::finter_writer_bzip2::writeline: bzip "
        "reports "
        "write operation called on read handle");
  }
}

void annotate_frequency::finter_writer_bzip2::write(char *buf,
                                                    std::streamsize n) {
  int error = 0;
  BZ2_bzWrite(&error, _bz_output, reinterpret_cast<void *>(buf),
              static_cast<int>(n));
  if (error == BZ_PARAM_ERROR || error == BZ_IO_ERROR) {
    _fail = true;
  } else if (error == BZ_SEQUENCE_ERROR) {
    throw std::domain_error(
        "annotate_frequency::finter_writer_bzip2::write: bzip reports "
        "write operation called on read handle");
  }
}
