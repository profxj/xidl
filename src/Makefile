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

#XIDL_DIR = $(XIDL_DIR)
LIB = $(XIDL_DIR)/lib

RECIPES = nrecipes.o nrutil.o
OBJECTS = \
	arclincorr.o \
	arrmedsigclip.o \
	arrmedsclipwgt.o \
	arrmedsigmask.o \
	rectify.o \
	vector_median.o \
	rebin2dspec.o \
	invert_arc.o \
	x_gsmooth.o \
	$(RECIPES)

all : $(LIB)/libxmath.$(SO_EXT)


$(LIB)/libxmath.$(SO_EXT): $(OBJECTS)
	$(LD) $(X_LD_FLAGS) -o $(LIB)/libxmath.$(SO_EXT) $(OBJECTS)
#	nm -s ../lib/libxmath.$(SO_EXT)

#
# Install things in their proper places in $(XIDL_DIR)
#
install : 
	-@ cp *.c $(XIDL_DIR)/src
	-@ cp *.h $(XIDL_DIR)/src
	-@ cp Makefile $(XIDL_DIR)/src

clean :
	- /bin/rm -f *~ core *.o so_locations
