# $Id: duneproject 5842 2010-01-20 18:48:34Z joe $

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

SUBDIRS = m4 dune doc

if BUILD_DOCS
# TODO: set up documentation tree automatically
#SUBDIRS += doc
endif

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) \
                            --with-dune-grid=$(DUNE_GRID_ROOT) \
                            --with-dune-istl=$(DUNE_ISTL_ROOT) \
                            --with-dune-cornerpoint=$(DUNE_CORNERPOINT_ROOT) \
                            --with-dune-porsol=$(DUNE_PORSOL_ROOT)  \
                            CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules
