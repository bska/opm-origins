# $Date$
# $Revision$

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

SUBDIRS = m4 common grid lib doc

ACLOCAL_AMFLAGS = -I m4

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5
# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) \
                            --with-dune-grid=$(DUNE_GRID_ROOT)  \
                            CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules
