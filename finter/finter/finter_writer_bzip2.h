/*!
  \file finter_writer_bzip2.h
  \brief bzip2-specific writer class definitions
  \copyright Released under the MIT License.
  Copyright 2020 Cameron Palmer
 */

#ifndef FINTER_FINTER_FINTER_WRITER_BZIP2_H_
#define FINTER_FINTER_FINTER_WRITER_BZIP2_H_

#include <bzlib.h>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include "finter/config.h"
#include "finter/finter_writer_parent.h"

namespace annotate_frequency {
/*!
  \class finter_writer_bzip2
  \brief interface for bzip2 file output that doesn't break my brain
 */
class finter_writer_bzip2 : public finter_writer {
 public:
  /*!
    \brief constructor
   */
  finter_writer_bzip2() : finter_writer(), _raw_output(0), _bz_output(0) {}
  /*!
    \brief destructor
   */
  ~finter_writer_bzip2() throw() { close(); }
  /*!
    \brief open a file
    @param filename name of file to open
   */
  void open(const char *filename);
  /*!
    \brief close the file
   */
  void close();
  /*!
    \brief clear internal state flags
   */
  void clear();
  /*!
    \brief test whether the file is currently open
    \return whether the file is currently open
   */
  bool is_open() const;
  /*!
    \brief write a single character to file
    @param c character to write to file
   */
  void put(char c);
  /*!
    \brief write a line to file, appending a system-specific newline
    @param s line to write to file
   */
  void writeline(const std::string &s);
  /*!
    \brief write a specified number of characters to file
    @param buf character buffer containing contents to be written
    @param n number of characters from buffer to write to file
   */
  void write(char *buf, std::streamsize n);
  /*!
    \brief test whether end of file has been reached
    \return whether end of file has been reached
   */
  bool eof() const { return false; }
  /*!
    \brief test whether connection is currently valid
    \return whether connection is currently valid
   */
  bool good() const { return _good; }
  /*!
    \brief test whether a write operation has failed for this stream
    \return whether a write operation has failed for this stream
   */
  bool fail() const { return _fail; }
  /*!
    \brief test whether connection is currently invalid
    \return whether connection is currently invalid
   */
  bool bad() const { return _bad; }

 private:
  FILE *_raw_output;   //!< C-style file handle used by bz2 library
  BZFILE *_bz_output;  //!< bz2 library interface to underlying file pointer
};
}  // namespace annotate_frequency

#endif  // FINTER_FINTER_FINTER_WRITER_BZIP2_H_
