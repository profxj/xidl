###############################################################################
# Sloan Digital Sky Survey (SDSS)
# IDL support code for products: idlmapper, idlspec2d
#
# S. Burles & D. Schlegel
###############################################################################

# For Macs with Tiger, copy libcc_dynamic.a into /usr/lib
# Then, do:  sudo ranlib libcc_dynamic.a
#
# IDL support utilities for spectro2d and the fibermapper
#
SHELL = /bin/sh
#
.c.o :
	$(CC) -c $(CCCHK) $(CFLAGS) $*.c
#
CFLAGS  = $(SDSS_CFLAGS) -DCHECK_LEAKS -I../include -DDARWIN

SUBDIRS = src

all :
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) all ); \
	done

#
# Install things in their proper places in $(XIDL_DIR)
#
install :
	@echo "You should be sure to have updated before doing this."
	@if [ "$(XIDL_DIR)" = "" ]; then \
		echo You have not specified a destination directory >&2; \
		exit 1; \
	fi 
	@echo "You will be installing in \$$XIDL_DIR=$$XIDL_DIR"
	@echo ""
	@echo "I'll give you 5 seconds to think about it"
	@sleep 5
	@echo
	@ rm -rf $(XIDL_DIR)
	@ mkdir $(XIDL_DIR)
	@ for f in $(SUBDIRS); do \
		(mkdir $(XIDL_DIR)/$$f; cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) install ); \
	done
	/bin/cp Makefile $(XIDL_DIR)


clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) clean ); \
	done
