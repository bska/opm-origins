# $Date$
# $Revision$

check_PROGRAMS = 
noinst_PROGRAMS = \
        aniso_implicit_steadystate_test \
        aniso_steadystate_test \
        implicit_steadystate_test \
        steadystate_test \
        upscaling_test

AM_CPPFLAGS += $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS) $(SUPERLU_CPPFLAGS)
AM_LDFLAGS  += $(DUNEMPILDFLAGS) $(BOOST_LDFLAGS) $(SUPERLU_LDFLAGS)
LDADD = $(DUNEMPILIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
        $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)   \
        $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SUPERLU_LIBS)


upscaling_test_SOURCES = upscaling_test.cpp

steadystate_test_SOURCES = steadystate_test.cpp

aniso_steadystate_test_SOURCES = aniso_steadystate_test.cpp

implicit_steadystate_test_SOURCES = implicit_steadystate_test.cpp

aniso_implicit_steadystate_test_SOURCES = aniso_implicit_steadystate_test.cpp

#upscale_perm_SOURCES = upscale_perm.C

#upscale_relperm_SOURCES = upscale_relperm.C

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
