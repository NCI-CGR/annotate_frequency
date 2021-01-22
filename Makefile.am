
SUBDIRS = finter

bin_PROGRAMS = annotate_frequency.out

AM_CXXFLAGS = -I$(srcdir)/finter -I$(srcdir)/finter/include -ggdb -Wall -std=c++17

annotate_frequency_out_SOURCES = annotate_frequency/cargs.cc annotate_frequency/cargs.h annotate_frequency/main.cc
annotate_frequency_out_LDADD = finter/lib/libfinter.la -lz -lbz2
dist_doc_DATA = README
ACLOCAL_AMFLAGS = -I m4
## TAP support
TEST_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(top_srcdir)/tap-driver.sh
TESTS = tests/placeholder.test
EXTRA_DIST = $(TESTS)
