#/* -*- Mode: Makefile -*- */
SHELL = /bin/sh

TOP=..

include $(TOP)/Path.incl 
include $(TOP)/$(SCIDIR)/Makefile.incl

NSP=$(TOP)/$(SCIDIR)/bin/nsp
#NAME=
LIB=libnspqhull

# if Makelib exists we run it since it can use parralel compilation
# then we always run the builder to be sure that loader.sce files are
# also generated

all::
	@if test -f Makelib; then \
		$(MAKE) -s -f Makelib $(MFLAGS); \
	fi
	@echo "running builder"; \
		$(NSP) -nw -e "exec('builder.sce');quit" -errcatch > /dev/null ;

tests	:: all

clean	::
	@if test -f Makelib; then make -s -f Makelib clean;fi

distclean::
	@if test -f Makelib; then make -s -f Makelib distclean;fi
	@$(RM) -f Makelib




