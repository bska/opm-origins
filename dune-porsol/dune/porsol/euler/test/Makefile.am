# $Date$
# $Revision: duneproject 5489 2009-03-25 11:19:24Z sander $

noinst_PROGRAMS = euler_upstream_test

AM_CPPFLAGS += $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS) $(SUPERLU_CPPFLAGS)
AM_LDFLAGS  += $(DUNEMPILDFLAGS) $(BOOST_LDFLAGS) $(SUPERLU_LDFLAGS)

LDADD = $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
        $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)   \
        $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)          \
        $(DUNEMPILIBS) $(SUPERLU_LIBS) 

euler_upstream_test_SOURCES = euler_upstream_test.cpp

#TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
