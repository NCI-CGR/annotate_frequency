/*!
  \file finter_writer_parent.h
  \brief virtual shared writer class declarations
  \copyright Released under the MIT License.
  Copyright 2020 Cameron Palmer
 */

#ifndef FINTER_FINTER_FINTER_WRITER_PARENT_H_
#define FINTER_FINTER_FINTER_WRITER_PARENT_H_

#include <ios>
#include <string>

#include "finter/helper.h"

namespace annotate_frequency {
/*!
  \class finter_writer
  \brief parent class for all finter writers
 */
class finter_writer {
 public:
  /*!
    \brief constructor
   */
  finter_writer() : _good(true), _bad(false), _fail(false) {}
  /*!
    \brief virtual destructor
   */
  virtual ~finter_writer() throw() {}
  /*!
    \brief open a file, converting name from string to char*
    @param filename name of file to open
   */
  void open(const std::string &filename) { open(filename.c_str()); }
  /*!
    \brief virtual file open
    @param filename of file to open
   */
  virtual void open(const char *filename) = 0;
  /*!
    \brief virtual file close
   */
  virtual void close() = 0;
  /*!
    \brief virtual clear state flags
   */
  virtual void clear() = 0;
  /*!
    \brief virtual test if file is currently open
    \return whether file is currently open
   */
  virtual bool is_open() const = 0;
  /*!
    \brief virtual write a single character to file
    @param c character to write to file
   */
  virtual void put(char c) = 0;
  /*!
    \brief virtual write a line to file, appending newline
    @param s line to write to file
   */
  virtual void writeline(const std::string &s) = 0;
  /*!
    \brief virtual write a specified number of characters to file
    @param buf array containing characters to write to file
    @param n how many characters in buffer to write to file
   */
  virtual void write(char *buf, std::streamsize n) = 0;
  /*!
    \brief virtual whether end of file has been reached
    \return whether end of file has been reached
   */
  virtual bool eof() const = 0;
  /*!
    \brief virtual test of whether stream is in valid state
    \return whether stream is in valid state
   */
  virtual bool good() const = 0;
  /*!
    \brief virtual test of whether write operation failed
    \return whether write operation failed
   */
  virtual bool fail() const = 0;
  /*!
    \brief virtual test of whether stream is in invalid state
    \return whether stream is in invalid state
   */
  virtual bool bad() const = 0;

 protected:
  bool _good;  //!< internal state flag: whether stream is in valid state
  bool _bad;   //!< internal state flag: whether stream is in invalid state
  bool _fail;  //!< internal state flag: whether past write operation has failed
};
}  // namespace annotate_frequency

#endif  // FINTER_FINTER_FINTER_WRITER_PARENT_H_
