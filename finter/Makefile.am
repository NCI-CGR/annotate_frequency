# Automake file

ACLOCAL_AMFLAGS = -I m4
# Build rules for libraries.
lib_LTLIBRARIES = lib/libfinter.la

lib_libfinter_la_SOURCES = finter/finter_reader_bzip2.cc \
  finter/finter_reader_flat.cc \
  finter/finter_writer_bzip2.cc \
  finter/finter_writer_flat.cc \
  finter/finter_reader_gzip.cc \
  finter/finter_writer_gzip.cc

library_includedir = $(includedir)/libfinter/include

library_include_HEADERS = \
  finter/finter_reader.h \
  finter/finter_reader_parent.h \
  finter/finter_reader_bzip2.h \
  finter/finter_reader_flat.h \
  finter/finter_reader_gzip.h \
  finter/finter_writer.h \
  finter/finter_writer_parent.h \
  finter/finter_writer_bzip2.h \
  finter/finter_writer_flat.h \
  finter/finter_writer_gzip.h \
  finter/helper.h \
  finter/finter.h \
  finter/config.h

lib_libfinter_la_libincludedir = $(libdir)/libfinter/include/finter

#nodist_lib_libfinter_la_libinclude_HEADERS = finter/config.h

lib_libfinter_la_LDFLAGS = -version-info 1:0:0

dist_doc_DATA = README
#ACLOCAL_AMFLAGS = -I m4

#lib_libfinter_la_libADD = lib/libfinter.la
