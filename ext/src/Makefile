###############################################################################
# Sloan Digital Sky Survey (SDSS)
# IDL support code for products: idlmapper, idlspec2d
#
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh
#
.c.o :
	$(CC) -c $(CCCHK) $(CFLAGS) $(X_CFLAGS) $*.c
#
CFLAGS  = -I$(IDL_DIR)/external $(SDSS_CFLAGS) -DCHECK_LEAKS

LIB = $(DJS_DIR)/lib

RECIPES = nrecipes.o nrutil.o
OBJECTS = \
	arravsigclip.o \
	arrmedsigclip.o \
	arravsigmask.o \
	qtst.o \
	arrmedsigmask.o \
	arrmedian.o \
	ccorrelate.o \
	$(RECIPES)

all : $(LIB)/libmath.$(SO_EXT)


$(LIB)/libmath.$(SO_EXT): $(OBJECTS)
	$(LD) $(X_LD_FLAGS) -o $(LIB)/libmath.$(SO_EXT) $(OBJECTS)
#	nm -s ../lib/libmath.$(SO_EXT)

#
# Install things in their proper places in $(DJS_DIR)
#
install : 
	-@ cp *.c $(DJS_DIR)/src/math
	-@ cp *.h $(DJS_DIR)/src/math
	-@ cp Makefile $(DJS_DIR)/src/math

clean :
	- /bin/rm -f *~ core *.o so_locations
